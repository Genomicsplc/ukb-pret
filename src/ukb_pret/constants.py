"""
Constants used throughout ukb-pret
"""
import os


def _add_metric_names(metric_names: list, generic_metrics: list):
    output_metrics = generic_metrics[:]
    for metric in metric_names:
        output_metrics += [metric, metric + '_l95ci', metric + '_u95ci']
    return output_metrics


PLOT_COLOURS = ['#E53935', '#9C27B0', '#3F51B5', '#2196F3', '#00BCD4', '#009688', '#43A047', '#FFC107', '#FF7043']
MIN_EVALUATION_SAMPLES_PER_GROUP = 100

####################################
# FILE PATHS
####################################

_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
PHENOTYPE_TRAIT_DEFS_PATH = os.path.join(_data_dir, 'traits.yaml')
POPULATION_CLUSTER_CENTRES_PATH = os.path.join(_data_dir, 'pop_cluster_centres.tsv')

####################################
# ANCESTRY INFERENCE PARAMETERS
####################################

SOFTMAX_BETA = 3

####################################
# DATA FILE HEADERS
####################################

SUPPORTED_PHENOTYPE_HEADERS = ['age_at_first_assessment', 'date_of_diagnosis', 'date_of_first_assessment',
                               'date_of_death', 'incident', 'sex', 'in_ukb_wbu_testing']
SUPPORTED_PC_HEADERS = [f'pc{i}' for i in range(1, 5)]
SURVIVAL_DATA_HEADERS = SUPPORTED_PHENOTYPE_HEADERS + ['age_at_diagnosis', 'age_at_censoring', 'age_event', 'prevalent']

####################################
# METRIC CONSTANTS
####################################

ODDS_RATIO_THRESHOLD_LIST = [3, 5, 10]
ODDS_RATIO_TITLE = 'odds_ratio_INtop{}'
GENERIC_METRIC_FIELDS = ['data_source', 'ancestry', 'sex', 'n_samples']
BINARY_METRIC_KEYS = ['auc', 'odds_ratio_1sd', 'hazard_ratio_1sd', 'cond_odds_ratio_1sd', 'cond_hazard_ratio_1sd'] + \
                     [ODDS_RATIO_TITLE.format(t) for t in ODDS_RATIO_THRESHOLD_LIST]
BINARY_METRIC_FIELDS = _add_metric_names(BINARY_METRIC_KEYS,
                                         GENERIC_METRIC_FIELDS + ['n_cases', 'n_controls'])
QUANTITATIVE_METRIC_KEYS = ['rsq', 'beta_1sd']
QUANTITATIVE_METRIC_FIELDS = _add_metric_names(QUANTITATIVE_METRIC_KEYS, GENERIC_METRIC_FIELDS)

####################################
# LABEL MAPPINGS
####################################

CUM_INC_X_LABEL_MAPPINGS = {'time_survival': 'Years since assessment', 'age_event': 'Age in years'}
SEX_MAPPING = {'male': 1., 'female': 0.}
REPORT_METRIC_MAPPINGS = {'auc': 'AUC', 'odds_ratio_1sd': 'OR per std', 'hazard_ratio_1sd': 'HR per std',
                          'rsq': 'Rsq', 'beta_1sd': 'Beta per std', 'cond_odds_ratio_1sd': 'Conditional OR',
                          'cond_hazard_ratio_1sd': 'Conditional HR'}
SUMMARY_TABLE_TO_RESULTS_MAPPING = {'Mean': 'mean', 'Mean LCI': 'mean_lci', 'Mean UCI': 'mean_uci',
                                    'Standard Deviation': 'std'}
SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_BINARY = {'Number of Samples': 'n_samples', 'Number of Cases': 'n_cases',
                                                 'Number of Controls': 'n_controls',
                                                 'Number of Incident Cases': 'n_incident',
                                                 'Number of Prevalent Cases': 'n_prevalent'}
SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_QUANTITATIVE = {'Number of Samples': 'n_samples'}
for t in ODDS_RATIO_THRESHOLD_LIST:
    REPORT_METRIC_MAPPINGS[ODDS_RATIO_TITLE.format(t)] = 'OR per std in top {} %'.format(t)

####################################
# REPORT CONSTANTS
####################################

REPORT_ANCESTRY_ORDER = ['EUR', 'EAS', 'SAS', 'AFR', 'All']
REPORT_SEX_ORDER = ['both', 'female', 'male']
BINARY_METRICS_TO_INCLUDE_IN_REPORT = ['cond_odds_ratio_1sd', 'cond_hazard_ratio_1sd', 'auc']
QUANTITATIVE_METRICS_TO_INCLUDE_IN_REPORT = ['rsq', 'beta_1sd']
