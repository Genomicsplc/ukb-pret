# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Disease Pheconstructor
#
# This notebook contains code to reproduce the phenotype for binary traits used in the UK Biobank PRS release.
#
# To reproduce the phenotype:
#
# 1. Convert this file to an Jupyter notebook using jupytext
# 2. Upload the notebook to your DNANexus project with `dx upload disease_pheconstructor.ipynb`
# 3. Upload the YAML of valid traits to your DNANexus project with `dx upload diseases.yaml`
# 4. Start a JupyterLab instance from the DNANexus web interface, ensuring that it is running on a Spark cluster
# 5. Open the pheconstructor notebook and replace the `TRAIT_CODE` variable in the first Python cell of this notebook with the trait of interest
# 6. Replace the `'<dataset_id>.dataset'` string in the first Python cell with the ID of your UKB dataset
# 7. Run the notebook - this should output the phenotype data to a CSV file named `<TRAIT_CODE>_phenotype.csv`, which is compatible with the UKB-PRET tool.
# If the Spark queries in the notebook fail, try increasing the number of cores and/or processing power of your instance

TRAIT_CODE = 'HT'
DATASET_ID = '<dataset_id>.dataset'
CONFIG_FILENAME = 'diseases.yaml'

# +
import pyspark
import dxdata

from typing import List
import pandas
import numpy
import yaml
# -

# Load the database into pyspark and dxdata
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
dataset = dxdata.load_dataset(DATASET_ID)

# +
# Load the phenotype configuration YAML

try:
    with open(f'/mnt/project/{CONFIG_FILENAME}', 'r') as f:
        pheno_config = yaml.safe_load(f)
except FileNotFoundError as e:
    raise FileNotFoundError(f'File {CONFIG_FILENAME} was not found. '
                            f'Please upload the file with: dx upload {CONFIG_FILENAME}')
# -

[pheno_config] = [conf for conf in pheno_config if conf['trait_code'] == TRAIT_CODE]


# ### Loading the self-reported data
#
# These data can be found in the 20002 field (for non-cancer cases), and the corresponding dates in the 20008 field
#
# Both fields are split into a number of instances (i0, i1, i2, i3), which denote the specific visit the patient made for assessment
#
# For 20002 instances, each row contains a list of categorical codes which indicate the specific disease they self reported
#
# Each 20008 instance consists of a number of arrays(a0, a1, ...), where each array element contains a single year for date of diagnosis

def load_all_self_reported_data(self_reported_config: dict):
        
    self_reported_code_prefix = f'p{self_reported_config["phenotype_datafield"]}'
    self_reported_date_code_prefix = f'p{self_reported_config["date_datafield"]}'
    disease_codes = self_reported_config['coding_values']

    [participant_data] = [e for e in dataset.entities if e.name == 'participant']
    sr_codes = [f.name for f in participant_data.fields if f.name.startswith(self_reported_code_prefix)]
    sr_names = ['eid'] + sr_codes

    # Load the self-reported data
    self_reported_data = participant_data.retrieve_fields(names=sr_names, engine=dxdata.connect())
    self_reported_df = self_reported_data.toPandas()
    self_reported_df = self_reported_df.set_index('eid')

    # Get the date of diagnosis for self reported data
    all_date_fields = [f for f in participant_data.fields if f.name.startswith(self_reported_date_code_prefix)]
    date_fields = [f.name for f in all_date_fields]
    all_fields = date_fields + ['eid']
    self_reported_data = participant_data.retrieve_fields(names=all_fields, engine=dxdata.connect())
    self_reported_data_df = self_reported_data.toPandas()
    self_reported_data_df = self_reported_data_df.set_index('eid')

    for suffix in ['_i0', '_i1', '_i2', '_i3']:
        instance_date_fields = [f for f in date_fields if suffix in f]
        self_reported_data_df['date_of_diagnosis' + suffix] = [[y for y in x if pandas.notna(y)] for x in self_reported_data_df[instance_date_fields].values.tolist()]
        self_reported_data_df['date_of_diagnosis' + suffix] = self_reported_data_df['date_of_diagnosis' + suffix].apply(lambda y: numpy.nan if len(y)==0 else y)

        # Add the date of diagnosis to self_reported_df
        self_reported_df['date_of_diagnosis' + suffix] = self_reported_data_df['date_of_diagnosis' + suffix]

    self_reported_df = self_reported_df.dropna(how='all')
    
    # Convert all remaining Nans / Nones to empty lists
    for c in self_reported_df.columns:
        self_reported_df[c] = self_reported_df[c].apply(lambda d: d if isinstance(d, list) else [])
    
    exploded_df = pandas.DataFrame(columns=['date_of_diagnosis', 'disease_code'])
    exploded_df.index.name = 'eid'

    for suffix in ['_i0', '_i1', '_i2', '_i3']:
        sub_df = self_reported_df[['date_of_diagnosis' + suffix, self_reported_code_prefix + suffix]].copy().reset_index()
        sub_df = sub_df[sub_df['date_of_diagnosis' + suffix].str.len() == sub_df[self_reported_code_prefix + suffix].str.len()]

        exploded_sub_df = sub_df.apply(pandas.Series.explode).reset_index().set_index('eid')
        exploded_sub_df['suffix'] = int(suffix[-1])
        exploded_df = pandas.concat([exploded_df,exploded_sub_df.rename(columns={'date_of_diagnosis' + suffix:'date_of_diagnosis', self_reported_code_prefix + suffix:'disease_code'})])
        
    # Filtering to the trait of interest
    disease_df = exploded_df[exploded_df['disease_code'].isin(disease_codes)]
    # Remove instances with unknown date (-1) or preferred not to answer (-3)
    disease_df['date_of_diagnosis'] = disease_df['date_of_diagnosis'].replace(-1, numpy.NaN).replace(-3, numpy.NaN)
    # Convert the date_of_diagnosis to a date
    disease_df['date_of_diagnosis'] = pandas.to_datetime(disease_df['date_of_diagnosis'].astype(float), format='%Y')
    
    # For each case, choose the earliest date of diagnosis
    disease_df = disease_df.reset_index()
    dates_series = disease_df.groupby('eid').apply(lambda x: x.loc[x['suffix'].idxmin()])['date_of_diagnosis']
    out_df = pandas.DataFrame({'date_of_diagnosis': dates_series, 'disease_code': 1})
    return out_df


# ### Loading the Hospital Inpatient data
# The `hesin` field contains the dates of diagnosis, and `hesin_diag` contains the ICD10 codes
#
# The dates in `hesin` are named `epistart` (which corresponds to the episode start date), and `admidate` (the admission date to the hospital) 
#
# Where available, we choose to use epistart as the date_of_diagnosis, and fall back onto admidate

def load_hesin_data(icd10_codes: List[str]):

    # Load the hesin pandas dataframe
    [hesin] = [e for e in dataset.entities if e.name == 'hesin']
    hes_data = hesin.retrieve_fields(names=['eid', 'epistart', 'admidate', 'ins_index'], engine=dxdata.connect())
    hes_df = hes_data.toPandas()
    hes_df = hes_df.set_index('eid')
    hes_df.index = hes_df.index.astype(str)
    
    # Convert dates to datetime format
    hes_df["epistart"] = pandas.to_datetime(hes_df.epistart, format="%Y-%m-%d")
    hes_df["admidate"] = pandas.to_datetime(hes_df.admidate, format="%Y-%m-%d")
    hes_df = hes_df.set_index("ins_index", append=True)
    
    # Collapse dates to single epistart column
    record_dates = hes_df.loc[:, ["epistart"]]
    record_dates.epistart = numpy.where(hes_df.epistart.isnull(), hes_df.admidate, hes_df.epistart)
    
    # Load the hesin_diag pandas dataframe
    [hesin_diag] = [e for e in dataset.entities if e.name == 'hesin_diag']
    hes_diag_data = hesin_diag.retrieve_fields(names=['eid', 'ins_index', 'diag_icd10', 'level'], engine=dxdata.connect())
    hes_diag_df = hes_diag_data.toPandas()
    hes_diag_df = hes_diag_df.set_index('eid')
    hes_diag_df.index = hes_diag_df.index.astype(str)
    hes_diag_df = hes_diag_df.set_index("ins_index", append=True)
    hes_diag_df = hes_diag_df.dropna()
    
    # Merging data from the two tables
    df = hes_diag_df.merge(record_dates, how="inner", left_index=True, right_index=True)
    df = df.reset_index()
    df = df.rename(columns={"level": "source", "epistart": "date_of_diagnosis", "diag_icd10": "code"})
    df = df.set_index("eid")
    
    # Get the date of last encounter of each individual
    df['date_of_last_encounter'] = df.groupby(df.index)["date_of_diagnosis"].max()

    df = df[df['code'].isin(icd10_codes)]
    date_of_diag_series = df.groupby(df.index)["date_of_diagnosis"].min()

    return date_of_diag_series.to_frame()


# ### Loading patient data
#
# Load the generic patient data. 
#
# _Note that the date_of_birth is estimated from the month and year of birth, set to be the 15th in all cases_

def load_binary_patient_data():
    
    [participants] = [e for e in dataset.entities if e.name == 'participant']

    date_of_death_code = 'p40000_i0' # Note that there is also a p40000_i1 column, but it's mostly empty (maybe for a few James Bonds'?)
    year_of_birth_code = 'p34'
    month_of_birth_code = 'p52'
    date_of_first_assessment_code = 'p53_i0'

    patient_data = participants.retrieve_fields(names=['eid', date_of_first_assessment_code, year_of_birth_code, month_of_birth_code, date_of_death_code], engine=dxdata.connect())
    patient_df = patient_data.toPandas()
    patient_df = patient_df[~patient_df[date_of_first_assessment_code].isna()]
    patient_df = patient_df.set_index('eid')
    patient_df = patient_df.rename(columns={date_of_first_assessment_code: 'date_of_first_assessment', date_of_death_code: 'date_of_death'})
    
    # Estimate date of birth from month and year of birth
    patient_df['date_of_birth'] = pandas.to_datetime(patient_df[year_of_birth_code].map(int).map(str) + '-' + patient_df[month_of_birth_code].map(int).map(str) + '-15', format="%Y-%m-%d")
    patient_df = patient_df.drop([year_of_birth_code, month_of_birth_code], axis=1)
    return patient_df


# # Generate Disease Phenotype
#
# Now a generic disease phenotype is created by combining self-reported data with hospital inpatient data

[self_reported_config] = [d for d in pheno_config['data_fields'] if d['data_type'] == 'self-report']
self_reported_df = load_all_self_reported_data(self_reported_config)

# +
[icd10_config] = [d for d in pheno_config['data_fields'] if d['data_type'] == 'icd10']

icd10_codes = icd10_config['coding_values']
hesin_df = load_hesin_data(icd10_codes)
# -

# Load the generic patient data
patient_df = load_binary_patient_data()

# Rename columns
self_reported_df = self_reported_df.rename(columns={'date_of_diagnosis':'date_of_diag_sr', 'disease_code':f'{pheno_config["trait_code"]}_sr'})
hesin_df = hesin_df.rename(columns={'date_of_diagnosis':'date_of_diag_icd10'})

# +
# Merging the patient data with hospital and self-reported data
patient_df = patient_df.merge(self_reported_df, left_index=True, right_index=True, how="left")
patient_df[f"{pheno_config['trait_code']}_sr"] = (~patient_df[f"{pheno_config['trait_code']}_sr"].isnull()).astype(int)

patient_df = patient_df.merge(hesin_df, left_index=True, right_index=True, how="left")
patient_df[f"{pheno_config['trait_code']}_main_icd10"] = (~patient_df["date_of_diag_icd10"].isnull()).astype(int)


# +
# Logical helper functions
def df_or(df):
    cols = [df[c] for c in df.columns]
    return vec_or(*cols)

def vec_or(*args):
    res = numpy.nan_to_num(args[0])
    for c in args[1:]:
        res = numpy.logical_or(res, numpy.nan_to_num(c))
    return numpy.array(res, dtype=int)


# -

# Creating the overall hypertension and date_of_diagnosis columns
# This is done with a logical OR on both columns, and choosing the minimum date of diagnosis
patient_df[pheno_config['trait_code']] = df_or(patient_df[[f'{pheno_config["trait_code"]}_sr', f'{pheno_config["trait_code"]}_main_icd10']])
patient_df['date_of_diagnosis'] = patient_df[["date_of_diag_sr", "date_of_diag_icd10"]].min(axis=1)

# Adding the incident column
patient_df['incident'] = (patient_df['date_of_diagnosis'] > patient_df["date_of_first_assessment"]).astype(int)

# Set the index 
patient_df.index.names = ['eid']
patient_df.index = patient_df.index.map(str)


# +
# Add the age_at_first_assessment

def delta_date(datename1, datename2):
    return numpy.abs((datename1 - datename2).dt.days)

def infer_age(birthdate, reference, newname="age_at_first_assessment"):
    return numpy.around(delta_date(birthdate, reference) / 365., 2)



# -

patient_df['age_at_first_assessment'] = infer_age(patient_df['date_of_birth'], pandas.to_datetime(patient_df['date_of_first_assessment']))

EXPECTED_BINARY_HEADERS = [pheno_config['trait_code'],'incident','age_at_first_assessment','date_of_first_assessment','date_of_diagnosis','date_of_death']
patient_df = patient_df[EXPECTED_BINARY_HEADERS]

# Output to file
patient_df.to_csv(f'{pheno_config["trait_code"]}_phenotype.csv')


