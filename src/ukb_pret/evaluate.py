"""
Evaluations and comparisons of PRS
"""

import logging
import os
import pandas
from typing import Tuple

from ._io import read_prs_file, read_pheno_file, read_ancestry_file, read_principal_components_file, \
    resolve_column_header
from .calculate import calculate_binary_metrics, calculate_quantitative_metrics, calculate_prs_sumstats, \
    calculate_sample_information
from .constants import ANCESTRY_MAPPINGS, SEX_MAPPING, SUPPORTED_PC_HEADERS, MIN_EVALUATION_SAMPLES_PER_GROUP,\
    BINARY_METRIC_KEYS, QUANTITATIVE_METRIC_KEYS
from .error import UkbPretImportError
from .plot import generate_binary_plots, generate_quantitative_plots, generate_cross_ancestry_plots


def read_and_prepare_data(prs_file_path1: str, prs_file_path2: str, pheno_file_path: str,
                          trait_code: str, phenotype_dict: dict, ancestry_file_path: str,
                          pcs_file_path: str) -> Tuple[pandas.DataFrame, pandas.DataFrame, dict, dict]:
    """Loads and prepares the input DataFrames from csv files. Returns the dataframes and dictionary of metadata"""

    df1, meta_dict1 = independent_sample_filtering(prs_file_path1, pheno_file_path, trait_code,
                                                   phenotype_dict, ancestry_file_path, pcs_file_path)
    df2, meta_dict2 = independent_sample_filtering(prs_file_path2, pheno_file_path, trait_code,
                                                   phenotype_dict, ancestry_file_path, pcs_file_path)

    df1, df2, n_removed_from_sample_intersection1, n_removed_from_sample_intersection2 = \
        _filter_to_input_data_intersection(df1, df2)

    meta_dict1['n_removed_from_input_intersection'] = {'value': n_removed_from_sample_intersection1,
                                                       'description': f'Number removed after intersection with the '
                                                                      f'other input PRS: '
                                                                      f'{n_removed_from_sample_intersection1}'}
    meta_dict1['n_included_samples'] = {'value': len(df1),
                                        'description': f'Total number of included samples: {len(df1)}'}

    meta_dict2['n_removed_from_input_intersection'] = {'value': n_removed_from_sample_intersection2,
                                                       'description': f'Number removed after intersection with the '
                                                                      f'other PRS data: '
                                                                      f'{n_removed_from_sample_intersection2}'}
    meta_dict2['n_included_samples'] = {'value': len(df2),
                                        'description': f'Total number of included samples: {len(df2)}'}

    return df1, df2, meta_dict1, meta_dict2


def independent_sample_filtering(prs_file: str, pheno_file_path: str, trait_code: str,
                                 phenotype_dict: dict, ancestry_file_path: str, pcs_file_path: str):
    prs_column_name = resolve_column_header(prs_file)
    prs_df = read_prs_file(prs_file, prs_column_name)
    n_input_samples = len(prs_df)
    pheno_df = read_pheno_file(pheno_file_path, trait_code, phenotype_dict)

    prs_df, pheno_df = _get_overlapping_results(prs_df, pheno_df)
    prs_df, pheno_df, n_prs_missing, n_pheno_missing = _filter_null_values(prs_df, pheno_df, trait_code)

    df = prs_df.join(pheno_df, on=prs_df.index)
    df = _prepare_ancestry_data(df, ancestry_file_path)
    df = _prepare_pcs_data(df, pcs_file_path)

    df, n_removed_from_ukb_wbu_filtering = _filter_to_ukb_wbu_testing(df)
    meta_dict = {'n_removed_prs': {'value': n_prs_missing,
                                   'description': f'Number of excluded PRS samples: {n_prs_missing}'},
                 'n_removed_pheno': {'value': n_pheno_missing,
                                     'description': f'Number of excluded phenotype samples: {n_pheno_missing}'},
                 'n_input_samples': {'value': n_input_samples,
                                     'description': f'Total number of input PRS samples: {n_input_samples}'},
                 'n_removed_from_ukb_wbu_filtering': {'value': n_removed_from_ukb_wbu_filtering,
                                                      'description': 'Number removed after filtering to the '
                                                                     'UKB PRS Release testing set: '
                                                                     f'{n_removed_from_ukb_wbu_filtering}'},

                 }

    return df, meta_dict


def evaluate_prs(df: pandas.DataFrame, trait_code: str, prs_column_name: str, phenotype_dict: dict, output_dir: str):
    """ Load PRS and compare performance"""

    eval_func = _determine_evaluation_method(phenotype_dict['scale'])

    eval_dict, cross_ancestry_eval_dict = dict(), dict()
    sexes_to_evaluate_in = ['both', 'female', 'male'] if phenotype_dict['gender'] == 'both' \
        else [phenotype_dict['gender']]
    for sex in sexes_to_evaluate_in:
        eval_dict[sex] = {}
        if sex == 'both':
            input_df = df
        else:
            input_df = df[df['sex'] == SEX_MAPPING[sex]].copy()

        eval_dict = eval_func(input_df, eval_dict, prs_column_name, trait_code, sex, phenotype_dict, output_dir)

        if all(p in df.columns for p in ['ancestry'] + SUPPORTED_PC_HEADERS):
            cross_ancestry_eval_dict[sex] = cross_ancestry_evaluation(df, prs_column_name, sex, phenotype_dict,
                                                                      output_dir)
        elif 'ancestry' in df.columns:
            cross_ancestry_eval_dict[sex] = cross_ancestry_evaluation(df, prs_column_name, sex,
                                                                      phenotype_dict, output_dir)

    return eval_dict, cross_ancestry_eval_dict


def _filter_to_ukb_wbu_testing(df: pandas.DataFrame):
    n_before = len(df)
    df = df[df['in_ukb_wbu_testing'] == 1]
    return df, n_before - len(df)


def _filter_to_input_data_intersection(df1: pandas.DataFrame, df2: pandas.DataFrame):
    """Filters the input samples to the intersection of values"""

    n_before_prs1, n_before_prs2 = len(df1), len(df2)
    df1_out = df1[df1.index.isin(df2.index)]
    df2_out = df2[df2.index.isin(df1.index)]
    return df1_out, df2_out, n_before_prs1 - len(df1_out), n_before_prs2 - len(df2_out)


def _determine_evaluation_method(trait_type: str):
    if trait_type == 'binary':
        eval_func = binary_evaluation
    elif trait_type == 'quantitative':
        eval_func = quantitative_evaluation
    else:
        raise Exception(f'Unsupported trait type: {trait_type}')

    return eval_func


def binary_evaluation(df: pandas.DataFrame, eval_dict: dict,
                      prs_column_name: str, trait_code: str, sex: str, phenotype_dict: dict, output_dir: str):

    if 'ancestry' in df.columns:
        ancestry_data = df.groupby(['ancestry']).apply(
            lambda x: generate_binary_results(x, prs_column_name, trait_code, sex, x['ancestry'].iloc[0],
                                              phenotype_dict, output_dir)
        )
        eval_dict[sex].update(ancestry_data.to_dict())

    else:
        eval_dict[sex]['All'] = generate_binary_results(df, prs_column_name, trait_code, sex, 'All', phenotype_dict,
                                                        output_dir)

    return eval_dict


def quantitative_evaluation(df: pandas.DataFrame, eval_dict: dict,
                            prs_column_name: str, trait_code: str, sex: str, phenotype_dict: dict, output_dir: str):
    if 'ancestry' in df.columns:
        ancestry_data = df.groupby(['ancestry']).apply(
            lambda x: generate_quantitative_results(x, prs_column_name, trait_code, sex, x['ancestry'].iloc[0],
                                                    phenotype_dict, output_dir)
        )
        eval_dict[sex].update(ancestry_data.to_dict())
    else:
        eval_dict[sex]['All'] = generate_quantitative_results(df, prs_column_name, trait_code, sex, 'All',
                                                              phenotype_dict, output_dir)

    return eval_dict


def _prepare_ancestry_data(df: pandas.DataFrame, ancestry_file_path: str):
    if ancestry_file_path is None:
        return df
    else:
        ancestry_df = read_ancestry_file(ancestry_file_path).dropna()
        for key, value in ANCESTRY_MAPPINGS.items():
            ancestry_df = ancestry_df.replace(key, value)
        # Remove AMR_NAT (too few samples)
        ancestry_df = ancestry_df[ancestry_df['ancestry'] != 'AMR']
        df_with_ancestry = df.join(ancestry_df)
        return df_with_ancestry


def _prepare_pcs_data(df: pandas.DataFrame, pcs_file_path: str):
    if pcs_file_path is None:
        return df
    else:
        # Prepare ancestry data
        pcs_df = read_principal_components_file(pcs_file_path).dropna()
        df = df.join(pcs_df)
        return df


def generate_binary_results(df: pandas.DataFrame, prs_column_name: str, trait_code: str,
                            sex: str, ancestry: str, phenotype_dict: dict, output_dir: str):
    sample_data = calculate_sample_information(df)
    prs_sumstats_dict = calculate_prs_sumstats(df, prs_column_name)
    if sample_data['n_cases'] > MIN_EVALUATION_SAMPLES_PER_GROUP and \
            sample_data['n_controls'] > MIN_EVALUATION_SAMPLES_PER_GROUP:
        metrics_dict = calculate_binary_metrics(df, prs_column_name, trait_code)
        plots_dict = generate_binary_plots(
            df, trait_code, prs_column_name, sex,
            ancestry, phenotype_dict, os.path.join(output_dir, 'plots', f'{prs_column_name}_{ancestry}'))
    else:
        logging.info(f'Cannot performing analysis in {ancestry} ancestry and {sex} sex - number of cases '
                     f'and/or controls is less than the minimum allowed value {MIN_EVALUATION_SAMPLES_PER_GROUP}')
        metrics_dict, plots_dict = {k: {'value': '', 'lci': '', 'uci': ''} for k in BINARY_METRIC_KEYS}, dict()
    return {'metrics': metrics_dict, 'plots': plots_dict, 'prs_sumstats': prs_sumstats_dict, 'sample_data': sample_data}


def generate_quantitative_results(df: pandas.DataFrame, prs_column_name: str, trait_code: str,
                                  sex: str, ancestry: str, phenotype_dict: dict, output_dir: str):

    prs_sumstats_dict = calculate_prs_sumstats(df, prs_column_name)
    sample_data = {'n_samples': len(df)}
    if sample_data['n_samples'] > MIN_EVALUATION_SAMPLES_PER_GROUP:
        metrics_dict = calculate_quantitative_metrics(df, prs_column_name, trait_code)
        plots_dict = generate_quantitative_plots(df, trait_code, prs_column_name,
                                                 os.path.join(output_dir, 'plots', f'{prs_column_name}_{ancestry}'),
                                                 sex, ancestry)
    else:
        logging.info(f'Cannot performing analysis in {ancestry} ancestry and {sex} sex - number of samples '
                     f'is less than the minimum allowed value {MIN_EVALUATION_SAMPLES_PER_GROUP}')
        metrics_dict, plots_dict = {k: {'value': '', 'lci': '', 'uci': ''} for k in QUANTITATIVE_METRIC_KEYS}, dict()

    return {'metrics': metrics_dict, 'plots': plots_dict, 'prs_sumstats': prs_sumstats_dict, 'sample_data': sample_data}


def cross_ancestry_evaluation(df: pandas.DataFrame, prs_column_name: str, sex: str, phenotype_dict: dict,
                              output_dir: str):
    if sex == phenotype_dict['gender']:
        plots_dict = generate_cross_ancestry_plots(df, prs_column_name,
                                                   os.path.join(output_dir,
                                                                'plots', f'{prs_column_name}_cross_ancestry'),
                                                   sex)
    else:
        plots_dict = dict()
    return {'plots': plots_dict}


def _get_overlapping_results(df_a: pandas.DataFrame, df_b: pandas.DataFrame, index: str = 'eid'):
    col_a, col_b = df_a.columns, df_b.columns
    for a in col_a:
        if a in col_b and a != index:
            raise UkbPretImportError(f'Column name {a} used in an input PRS file and a phenotype file - please use '
                                     f'distinct column names (excluding {index})')
    combined_df = df_a.join(df_b, how='inner')
    df_a_new, df_b_new = combined_df[col_a], combined_df[col_b]

    return df_a_new, df_b_new


def _filter_null_values(prs_df: pandas.DataFrame, pheno_df: pandas.DataFrame, trait_code: str):
    """Logic for processing null values"""

    # TODO: Specifics on dealing with nulls TBC. For now simply removing nulls
    prs_df_clean = prs_df.dropna()
    pheno_df_clean = pheno_df[pheno_df[trait_code].notna()]
    removed_prs_entries = len(prs_df) - len(prs_df_clean)
    removed_pheno_entries = len(pheno_df) - len(pheno_df_clean)

    logging.info(f'Removed {removed_prs_entries} null PRS entries')
    logging.info(f'Removed {removed_pheno_entries} null phenotype entries')

    out_prs_df, out_pheno_df = _get_overlapping_results(prs_df_clean, pheno_df_clean)

    return out_prs_df, out_pheno_df, removed_prs_entries, removed_pheno_entries
