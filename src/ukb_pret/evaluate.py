"""
Evaluations and comparisons of PRS
"""

import logging
import numpy
import os
import pandas
from scipy.special import softmax
from sklearn.metrics.pairwise import cosine_similarity
from typing import Tuple, List

from ._io import read_pop_cluster_centers, write_metrics_csv
from .calculate import calculate_binary_metrics, calculate_quantitative_metrics, calculate_prs_sumstats, \
    calculate_sample_information
from .constants import SEX_MAPPING, SUPPORTED_PC_HEADERS, MIN_EVALUATION_SAMPLES_PER_GROUP,\
    QUANTITATIVE_METRIC_KEYS, POPULATION_CLUSTER_CENTRES_PATH, SOFTMAX_BETA, \
    SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS, BINARY_METRIC_KEYS_NO_SURVIVAL, SURVIVAL_DATA_HEADERS

from .error import UkbPretImportError
from .plot import generate_binary_plots, generate_quantitative_plots, generate_cross_ancestry_plots
from .report import generate_report


def prepare_data(prs1_df: pandas.DataFrame, prs2_df: pandas.DataFrame, pheno_df: pandas.DataFrame,
                 pcs_df: pandas.DataFrame, trait_code: str) -> Tuple[pandas.DataFrame, pandas.DataFrame, dict, dict]:
    """Loads and prepares the input DataFrames from csv files. Returns the dataframes and dictionary of metadata"""

    pop_cluster_centers = read_pop_cluster_centers(POPULATION_CLUSTER_CENTRES_PATH)

    df1, meta_dict1 = independent_sample_filtering(prs1_df, pheno_df, pcs_df, trait_code)
    df2, meta_dict2 = independent_sample_filtering(prs2_df, pheno_df, pcs_df, trait_code)

    df1, df2, n_removed_from_sample_intersection1, n_removed_from_sample_intersection2 = \
        _filter_to_input_data_intersection(df1, df2)

    meta_dict1['n_removed_from_input_intersection'] = {'value': n_removed_from_sample_intersection1,
                                                       'description': f'Number removed after intersection with the '
                                                                      f'other input PRS: '
                                                                      f'{n_removed_from_sample_intersection1}'}
    meta_dict1['n_included_samples'] = {'value': len(df1),
                                        'description': f'Total number of included samples: {len(df1)}'}

    meta_dict1 = _add_total_removed_from_prs(meta_dict1)

    meta_dict2['n_removed_from_input_intersection'] = {'value': n_removed_from_sample_intersection2,
                                                       'description': f'Number removed after intersection with the '
                                                                      f'other PRS data: '
                                                                      f'{n_removed_from_sample_intersection2}'}
    meta_dict2['n_included_samples'] = {'value': len(df2),
                                        'description': f'Total number of included samples: {len(df2)}'}

    meta_dict2 = _add_total_removed_from_prs(meta_dict2)

    df1, df2 = infer_all_ancestries([df1, df2], pop_cluster_centers)
    return df1, df2, meta_dict1, meta_dict2


def _add_total_removed_from_prs(meta_dict: dict) -> dict:
    n_prs_missing_tot = meta_dict['n_removed_prs']['value'] + meta_dict['n_removed_from_input_intersection']['value']
    meta_dict['n_removed_prs_tot'] = {'value': n_prs_missing_tot,
                                      'description': f'Number of excluded PRS samples: {n_prs_missing_tot}'}
    return meta_dict


def infer_all_ancestries(dfs: List[pandas.DataFrame], pop_cluster_centers: pandas.DataFrame):
    out_dfs = []
    for df in dfs:
        out_dfs.append(_infer_ancestry_from_pcs(df, pop_cluster_centers))
    return out_dfs


def _infer_ancestry_from_pcs(df: pandas.DataFrame, pop_cluster_centres: pandas.DataFrame) -> pandas.DataFrame:

    if not all(pc in df.columns for pc in SUPPORTED_PC_HEADERS):
        return df

    # Calculate the cosine similarity between the population PC centers and the input PC data
    ancestries = pop_cluster_centres.index.values
    pc_data = df[SUPPORTED_PC_HEADERS]
    cos_sim = cosine_similarity(pc_data, pop_cluster_centres)
    proj_proportion = softmax(cos_sim * SOFTMAX_BETA, axis=1)

    proj_pred_ancestry = numpy.array([ancestries[i] for i in proj_proportion.argmax(axis=1)], dtype='object')
    proj_pred_ancestry_df = pandas.DataFrame(index=df.index, data=proj_pred_ancestry,
                                             columns=['ancestry']).rename_axis(index='eid')

    # Remove AMR due to insufficient samples in UKB
    proj_pred_ancestry_df = proj_pred_ancestry_df[proj_pred_ancestry_df['ancestry'] != 'AMR']
    df_with_ancestry = df.join(proj_pred_ancestry_df)
    return df_with_ancestry


def independent_sample_filtering(prs_df: pandas.DataFrame, pheno_df: pandas.DataFrame, pcs_df: pandas.DataFrame,
                                 trait_code: str):

    n_input_samples = len(prs_df)

    prs_df, pheno_df = _get_overlapping_results(prs_df, pheno_df)
    prs_df, pheno_df, n_prs_missing, n_pheno_missing = _filter_null_values(prs_df, pheno_df, trait_code)

    df = prs_df.join(pheno_df, on=prs_df.index)
    if pcs_df is None:
        logging.warn('No Principal Component (PC) file provided - analysis will be carried out across all individuals '
                     'and there will be no quality control section.')
    else:
        df = df.join(pcs_df)

    meta_dict = {'n_removed_prs': {'value': n_prs_missing,
                                   'description': f'Number of excluded PRS samples from overlap with :'
                                                  f' {n_prs_missing}'},
                 'n_removed_pheno': {'value': n_pheno_missing,
                                     'description': f'Number of excluded phenotype samples: {n_pheno_missing}'},
                 'n_input_samples': {'value': n_input_samples,
                                     'description': f'Total number of input PRS samples: {n_input_samples}'}
                 }

    return df, meta_dict


def compare_prs(prs1_df: pandas.DataFrame, prs2_df: pandas.DataFrame, pheno_df: pandas.DataFrame, metadata: dict,
                pcs_df=None):

    eval_dict, cross_ancestry_eval_dict = dict(), dict()
    prs1_df, prs2_df, meta_dict1, meta_dict2 = prepare_data(prs1_df, prs2_df, pheno_df, pcs_df, metadata['trait_code'])
    meta_dict = {metadata['prs1_label']: meta_dict1, metadata['prs2_label']: meta_dict2}

    eval_dict[metadata['prs1_label']], cross_ancestry_eval_dict[metadata['prs1_label']] = \
        evaluate_prs(prs1_df, metadata['trait_code'], metadata['prs1_label'],
                     metadata['phenotype_dictionary'], metadata['output_dir'])

    eval_dict[metadata['prs2_label']], cross_ancestry_eval_dict[metadata['prs2_label']] = \
        evaluate_prs(prs2_df, metadata['trait_code'], metadata['prs2_label'],
                     metadata['phenotype_dictionary'], metadata['output_dir'])

    results = construct_results_dict(eval_dict, metadata['phenotype_dictionary'], cross_ancestry_eval_dict, meta_dict,
                                     metadata['trait_code'])
    results['survival_data'] = all(x in pheno_df.columns for x in SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS)
    write_metrics_csv(results, os.path.join(metadata['output_dir'], 'evaluation_metrics.csv'))
    generate_report(results, os.path.join(metadata['output_dir'], 'prs_evaluation_report.pdf'))


def evaluate_prs(df: pandas.DataFrame, trait_code: str, prs_column_name: str, phenotype_dict: dict, output_dir: str):
    """ Load PRS and compare performance"""

    eval_func = _determine_evaluation_method(phenotype_dict['scale'])

    eval_dict, cross_ancestry_eval_dict = dict(), dict()

    if 'sex' not in df.columns:
        df['sex'] = 'unspecified'
        sexes_to_evaluate_in = ['unspecified']
    else:
        sexes_to_evaluate_in = ['both', 'female', 'male'] if phenotype_dict['sex'] == 'both' \
            else [phenotype_dict['sex']]

    for sex in sexes_to_evaluate_in:
        eval_dict[sex] = {}
        if sex in ['both', 'unspecified']:
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


def _are_sufficient_samples(sample_data: dict) -> bool:
    sufficient_samples = False
    if sample_data['n_cases'] is None and sample_data['n_controls'] is None:
        if sample_data['n_samples'] > MIN_EVALUATION_SAMPLES_PER_GROUP:
            sufficient_samples = True
    else:
        if sample_data['n_cases'] > MIN_EVALUATION_SAMPLES_PER_GROUP and \
                sample_data['n_controls'] > MIN_EVALUATION_SAMPLES_PER_GROUP:
            sufficient_samples = True
    return sufficient_samples


def resolve_binary_headers(df: pandas.DataFrame):
    """
    Resolve the headers expected in output metrics from input dataframe
    """

    metrics = BINARY_METRIC_KEYS_NO_SURVIVAL
    if all(x in df.columns for x in ['age_at_first_assessment', 'sex']):
        metrics.append('cond_odds_ratio_1sd')

    if all(h in df.columns for h in SURVIVAL_DATA_HEADERS):
        metrics.append('hazard_ratio_1sd')

        if all(x in df.columns for x in ['age_at_first_assessment', 'sex']):
            metrics.append('cond_hazard_ratio_1sd')

    return metrics


def generate_binary_results(df: pandas.DataFrame, prs_column_name: str, trait_code: str,
                            sex: str, ancestry: str, phenotype_dict: dict, output_dir: str):
    sample_data = calculate_sample_information(df)
    prs_sumstats_dict = calculate_prs_sumstats(df, prs_column_name)
    if _are_sufficient_samples(sample_data):
        metrics_dict = calculate_binary_metrics(df, prs_column_name, trait_code)
        plots_dict = generate_binary_plots(
            df, trait_code, prs_column_name, sex,
            ancestry, phenotype_dict, os.path.join(output_dir, 'plots', f'{prs_column_name}_{ancestry}'))
    else:
        logging.info(f'Cannot performing analysis in {ancestry} ancestry and {sex} sex - number of cases '
                     f'and/or controls is less than the minimum allowed value {MIN_EVALUATION_SAMPLES_PER_GROUP}')
        metrics_dict, plots_dict = {k: {'value': '', 'lci': '', 'uci': ''} for k in resolve_binary_headers(df)}, dict()
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
    if sex == phenotype_dict['sex'] or sex == 'unspecified':
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


def construct_results_dict(eval_dict: dict, phenotype_dict: dict, cross_ancestry_eval_dict: dict, meta_dict: dict,
                           trait_code: str):

    results = {'phenotype_code': trait_code,
               'phenotype_description': phenotype_dict['full_name'],
               'flavour': phenotype_dict['scale'],
               'sex_major': phenotype_dict['sex'],
               'evaluation': eval_dict,
               'cross_ancestry_evaluation': cross_ancestry_eval_dict,
               'metadata': meta_dict
               }
    return results
