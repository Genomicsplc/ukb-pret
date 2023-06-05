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

# # Quantitative Pheconstructor
#
# This notebook contains code to reproduce the phenotype for quantitative traits used in the UK Biobank PRS release.
#
# To reproduce the phenotype:
#
# 1. Convert this file to an Jupyter notebook using jupytext
# 2. Upload the notebook to your DNANexus project with `dx upload qt_pheconstructor.ipynb`
# 3. Upload the YAML of valid traits to your DNANexus project with `dx upload qt.yaml`
# 4. Start a JupyterLab instance from the DNANexus web interface, ensuring that it is running on a Spark cluster
# 5. Open the pheconstructor notebook and replace the `TRAIT_CODE` variable in the first Python cell of this notebook with the trait of interest
# 6. Replace the `'<dataset_id>.dataset'` string in the first Python cell with the ID of your UKB dataset
# 7. Run the notebook - this should output the phenotype data to a CSV file named `<TRAIT_CODE>_phenotype.csv`, which is compatible with the UKB-PRET tool.
# If the Spark queries in the notebook fail, try increasing the number of cores and/or processing power of your instance

TRAIT_CODE = 'PDCL'
DATASET_ID = '<dataset_id>.dataset'
CONFIG_FILENAME = 'qt.yaml'

# +
import pyspark
import dxdata

import yaml
# -

# Load the database into pyspark and dxdata
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
dataset = dxdata.load_dataset(DATASET_ID)

# Load the phenotype configuration YAML
try:
    with open(f'/mnt/project/{CONFIG_FILENAME}', 'r') as f:
        pheno_config = yaml.safe_load(f)
except FileNotFoundError as e:
    raise FileNotFoundError(f'File {CONFIG_FILENAME} was not found. '
                            f'Please upload the file with: dx upload {CONFIG_FILENAME}')

[pheno_config] = [conf for conf in pheno_config if conf['trait_code'] == TRAIT_CODE]
[ukb_code_config] = [d for d in pheno_config['data_fields'] if d['data_type'] == 'quantitative']

# +
# Load the UKB data fields corresponding to the trait

ukb_code = ukb_code_config['phenotype_datafield']
[participant_data] = [e for e in dataset.entities if e.name == 'participant']
codes = [f.name for f in participant_data.fields if f.name.startswith(f'p{ukb_code}')]
names = ['eid'] + codes

data = participant_data.retrieve_fields(names=names, engine=dxdata.connect())
df = data.toPandas()
df = df.set_index('eid')

# -

# Remove any instances where all instances for an individual are NaN
df = df.dropna(how='all')

# Take the value of the first instance where available, and the second instance otherwise
df[TRAIT_CODE] = df[f'p{ukb_code}_i0'].fillna(df[f'p{ukb_code}_i1'])

# Select the trait column only and write the phenotype to CSV
out_df = df[TRAIT_CODE].to_frame()
out_df.to_csv(f'{pheno_config["trait_code"]}_phenotype.csv')




