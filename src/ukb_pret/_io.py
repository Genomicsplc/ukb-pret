"""
Scripts for reading and writing PRS scores for evaluation
"""

from datetime import timedelta
import logging
import numpy
import os
import pandas
import sys
import yaml

from .constants import SUPPORTED_PHENOTYPE_HEADERS, SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS,\
    BINARY_METRIC_FIELDS, BINARY_METRIC_NO_SURVIVAL_FIELDS, QUANTITATIVE_METRIC_FIELDS, SUPPORTED_PC_HEADERS,\
    PHENOTYPE_TRAIT_DEFS_PATH, UKB_ENHANCED_PRS_CODES_PATH, UKB_RELEASE_PRS_FILE_HEADERS


def _load_csv_file(path: str, delimiter: str = ','):
    return pandas.read_csv(path, delimiter=delimiter)


def _add_survival_vars(df: pandas.DataFrame, outcome_phenotype: str, follow_up_months: float = 120,
                       censoring_date: str = '2020-03-01') -> pandas.DataFrame:
    df = _add_prevalent_data(df, outcome_phenotype)
    recode = ["date_of_death", "date_of_first_assessment", "date_of_diagnosis"]
    for i in recode:
        df[i] = pandas.to_datetime(df[i])

    follow_up_years = follow_up_months / 12
    days_to_censor = round((365 / 12) * follow_up_months)
    df = _add_censoring_data(df, censoring_date, days_to_censor)

    df['age_at_diagnosis'] = df['age_at_first_assessment'] + ((df.date_of_diagnosis -
                                                               df.date_of_first_assessment) / timedelta(days=365))
    df['date_of_last_encounter'] = pandas.to_datetime(df[['date_of_death', 'censoring_date']].min(axis=1))
    df['age_at_censoring'] = df['age_at_first_assessment'] + ((df.date_of_last_encounter -
                                                               df.date_of_first_assessment) / timedelta(days=365))
    df['age_event'] = df[['age_at_diagnosis', 'age_at_censoring']].min(axis=1)
    df.loc[(df[outcome_phenotype] == 1) & numpy.isnat(
        df['date_of_diagnosis']), 'age_event'] = df['age_at_first_assessment']

    df['time_survival'] = (df['censoring_date'] -
                           df['date_of_first_assessment']) / timedelta(days=365)

    df.loc[df['date_of_diagnosis'].notna(), 'time_survival'] = (
            (df['date_of_diagnosis'] - df['date_of_first_assessment']) / timedelta(days=365))

    df.loc[(df['date_of_death'].notna()) &
           (df['date_of_death'] < df['date_of_diagnosis']), 'time_survival'] = (
            (df['date_of_death'] - df['date_of_first_assessment']) / timedelta(days=365))

    df[outcome_phenotype] = numpy.where(
        (df['time_survival'] > follow_up_years) & (df['incident'] == 1), 0, df[outcome_phenotype])

    df['time_survival'] = numpy.where(df['prevalent'] == 1, numpy.nan, df['time_survival'])

    df['incident'] = numpy.where(
        (df['time_survival'] > follow_up_years) & (df['incident'] == 1), 0, df['incident']
    )
    df['time_survival'] = numpy.where(
        (df['time_survival'] > follow_up_years), follow_up_years, df['time_survival']
    )
    return df.dropna(subset=['age_event', outcome_phenotype])


def _add_censoring_data(df: pandas.DataFrame, censoring_date: str, days_to_censor: float):
    df['censoring_date'] = pandas.to_datetime(censoring_date)
    df['censoring_date_follow_up'] = df.date_of_first_assessment + timedelta(days=days_to_censor)
    df['censoring_date'] = pandas.to_datetime(
        df[['censoring_date_follow_up', 'censoring_date']].min(axis=1))
    return df


def _add_prevalent_data(df, trait_code):
    try:
        df.loc[(df[trait_code] == 1) & (df['incident'] == 1), 'prevalent'] = 0
        df.loc[(df[trait_code] == 1) & (df['incident'] == 0), 'prevalent'] = 1
        df.loc[(df[trait_code] == 0) & (df['incident'] == 0), 'prevalent'] = 0
        df.loc[(df[trait_code] == 0) & (df['incident'] == 1), 'prevalent'] = numpy.nan
        df.loc[(df[trait_code] == 0) & (df['incident'] == 1), 'incident'] = numpy.nan
        df.loc[(df['prevalent'] == 1) & (df['incident'] == 0), 'incident'] = numpy.nan

        prev_miss = df['prevalent'].isnull().sum()
        age_miss = df['age_at_first_assessment'].isnull().sum()
        if prev_miss > 0:
            logging.info("Number of individuals dropped due to incorrect disease coding: {}".format(
                prev_miss))
        if age_miss > 0:
            logging.info("Number of individuals dropped due to missing age: {}".format(
                format(age_miss)))
        df.dropna(axis=0, subset=['age_at_first_assessment', 'prevalent'], inplace=True)
    except BaseException as e:
        logging.info(f"Prevalent and Incident phenotype data not generated...")
        raise e
    return df


def load_phenotype_dictionary():
    """Load supported phenotype definitions from file."""
    return read_traits_file(PHENOTYPE_TRAIT_DEFS_PATH)


def load_ukb_enhanced_prs_code(input_path: str):
    codes_df = read_ukb_enhanced_prs_codes_file(UKB_ENHANCED_PRS_CODES_PATH)
    ukb_code = int(resolve_column_header(input_path, ignored_headers=UKB_RELEASE_PRS_FILE_HEADERS).lstrip('p'))
    return codes_df.loc[ukb_code]['gplc_trait_code'], ukb_code


def read_ukb_enhanced_prs_codes_file(input_path: str, id_field: str = 'ukb_enhanced_prs_code', delimiter: str = '\t'):
    df = _load_csv_file(input_path, delimiter)
    expected_headers = {'gplc_trait_code', 'ukb_field_description', id_field}
    if set(df.columns) != expected_headers:
        print(f'The input file {input_path} should only contain columns with headers {expected_headers}')
        sys.exit(1)
    return df.set_index(id_field)


def read_prs_file(input_path: str, prs_field: str, id_field: str = 'eid', delimiter: str = ','):
    """Loads the PRS from file. Resolve the data_tag from the input column header"""

    df = _load_csv_file(input_path, delimiter)
    expected_headers = {prs_field, id_field}
    if set(df.columns) != expected_headers:
        print(f'The input file {input_path} should only contain columns with headers {expected_headers}')
        sys.exit(1)
    return df.set_index(id_field)


def read_ukb_prs_file(input_path: str,  ukb_code: int, id_field: str = 'eid', delimiter: str = ','):

    df = _load_csv_file(input_path, delimiter)
    expected_headers = set(UKB_RELEASE_PRS_FILE_HEADERS + [id_field, f'p{str(ukb_code)}'])
    if set(df.columns) != expected_headers:
        print(f'The input file {input_path} should only contain columns with headers {expected_headers}')
        sys.exit(1)
    return df.set_index(id_field)


def read_pheno_file(input_path: str, pheno_field: str, phenotype_dict: dict,
                    id_field: str = 'eid', delimiter: str = ',') -> pandas.DataFrame:
    """Loads the phenotypes from file"""

    df = _load_csv_file(input_path, delimiter)
    standard_headers = {pheno_field, id_field}
    standard_headers_with_sex = standard_headers | {'sex'}
    extended_headers = standard_headers | set(SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS)
    extended_headers_with_sex = extended_headers | {'sex'}

    supported_headers = [standard_headers, standard_headers_with_sex, extended_headers, extended_headers_with_sex]
    if not any(hdrs == set(df.columns) for hdrs in supported_headers):
        print(f'The input file {input_path} should only contain columns with headers in one of the '
              f'following supported formats:')
        for headers in supported_headers:
            print(headers)
        sys.exit(1)

    df = df.set_index(id_field)
    if set(df.columns) | {id_field} == extended_headers or set(df.columns) | {id_field} == extended_headers_with_sex:
        df = _add_survival_vars(df, pheno_field, float(phenotype_dict['follow_up_months']))
    return df


def read_pop_cluster_centers(input_path: str, id_field: str = 'population', delimiter: str = '\t'):
    df = _load_csv_file(input_path, delimiter)
    expected_headers = set(SUPPORTED_PC_HEADERS) | {id_field}
    if set(df.columns) != expected_headers:
        print(f'The input file {input_path} should only contain columns with headers {expected_headers}')
        sys.exit(1)
    return df.set_index(id_field)


def read_principal_components_file(input_path: str, id_field: str = 'eid', delimiter: str = ','):
    """Read the principal components file containing the first 4 genetically inferred pcs [eid,pc1,pc2,pc3,pc4]"""
    if input_path is None:
        return pandas.DataFrame()
    df = _load_csv_file(input_path, delimiter)
    expected_headers = set(SUPPORTED_PC_HEADERS) | {id_field}
    if set(df.columns) != expected_headers:
        print(f'The input file {input_path} should only contain columns with headers {expected_headers}')
        sys.exit(1)
    return df.set_index(id_field)


def read_traits_file(input_path: str):
    with open(input_path, 'r') as f:
        traits_data = yaml.load(f, Loader=yaml.BaseLoader)
    return traits_data


def resolve_column_header(input_path: str, id_field: str = 'eid', delimiter: str = ',',
                          ignored_headers: list = SUPPORTED_PHENOTYPE_HEADERS):
    # TODO: Re-organise such that only one read from the file is performed

    df = pandas.read_csv(input_path, delimiter=delimiter)
    df = df[[x for x in df.columns if x not in ignored_headers]]
    if len(df.columns) != 2:
        print(f'The file {input_path} contains some unsupported headers. Please compare '
              f'{list(df.columns)} to the expected inputs')
        sys.exit(1)

    if id_field not in df.columns:
        print(f'No field matching the expected ID column {id_field} found in {input_path}')
        sys.exit(1)

    [header] = [f for f in df.columns if f != id_field]
    return header


def resolve_rap_inputs(ukb_release_prs_file_path: str):
    """
    This function contains all logic for resolving dataframes from Research Analysis Platform output
    """

    gplc_trait_code, ukb_code = load_ukb_enhanced_prs_code(ukb_release_prs_file_path)
    ukb_df = read_ukb_prs_file(ukb_release_prs_file_path, ukb_code)
    ukb_prs_df = ukb_df[[f'p{ukb_code}']]
    ukb_prs_label = f'UKB_Enhanced_PRS_{gplc_trait_code}'
    ukb_prs_df = ukb_prs_df.rename(columns={f'p{ukb_code}': ukb_prs_label})
    pcs_df = ukb_df[SUPPORTED_PC_HEADERS]
    sex_df = ukb_df['sex']
    return ukb_prs_df, pcs_df, sex_df, ukb_prs_label, gplc_trait_code


def write_metrics_csv(results: dict, output_path: str):
    """
    Write out a CSV file of evaluation metrics
    """

    if results['flavour'] == 'binary':
        valid_fields = BINARY_METRIC_FIELDS if results['survival_data'] else BINARY_METRIC_NO_SURVIVAL_FIELDS
    else:
        valid_fields = QUANTITATIVE_METRIC_FIELDS

    csv_dict = {metric: [] for metric in valid_fields}
    for data_source, by_sex in results['evaluation'].items():
        for sex, by_anc in by_sex.items():
            for ancestry, data in by_anc.items():
                for metric_key, metric in data.get('metrics', {}).items():
                    # if metric_key in valid_fields:
                    csv_dict[metric_key].append(metric['value'])
                    csv_dict[metric_key + '_l95ci'].append(metric['lci'])
                    csv_dict[metric_key + '_u95ci'].append(metric['uci'])

                csv_dict['ancestry'].append(ancestry)
                csv_dict['sex'].append(sex)
                csv_dict['data_source'].append(data_source)
                csv_dict['n_samples'].append(data['sample_data']['n_samples'])
                if results['survival_data']:
                    csv_dict['n_cases'].append(data['sample_data']['n_cases'])
                    csv_dict['n_controls'].append(data['sample_data']['n_controls'])
    df = pandas.DataFrame.from_dict(csv_dict).set_index('data_source')
    df.to_csv(output_path)
