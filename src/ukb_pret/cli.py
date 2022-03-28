"""
Command line interface for evaluating PRS
"""

import argparse
import os
import sys

from ._io import resolve_column_header, write_metrics_csv, load_phenotype_dictionary
from .constants import SUPPORTED_PHENOTYPE_HEADERS
from .evaluate import evaluate_prs, read_and_prepare_data
from .report import generate_report


class ExitWithList(argparse.Action):

    def __init__(self, option_strings, dest, **kwargs):
        self.values = None
        super().__init__(option_strings, dest, nargs=0, **kwargs)

    def __call__(self, *args, **kwargs):
        fields = '-  ' + '\n-  '.join(self.values)
        print(f"Listing available fields for {self.option_strings[0]}:\n\n"
              f"{fields}\n\n"
              f"Please run without {self.option_strings[0]} to make a query.")
        sys.exit(0)


class ListPhenotypes(ExitWithList):
    """Lists the supported phenotype definitions provided by Genomics plc"""

    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)
        self.values = _parse_supported_phenotypes()


# TODO: Provide interface with Gplc supported phenotype definitions
def _parse_supported_phenotypes():
    # "Config containing supported phenotypes has not yet been implemented")
    return []


def build_parser():
    """Defines the arguments taken by the command line parser."""

    parser = argparse.ArgumentParser(description="A CLI tool for evaluating a set of PRS against a set of phenotypes.")

    parser.add_argument('--prs-files',
                        nargs='+',
                        required=True,
                        help='Paths to two files containing Polygenic Risk Score (PRS) and '
                             'participant eIDs in CSV format. Headers should be [eid,<data_tag>], '
                             'where <data_tag> is the field used to identify the PRS in the output')
    parser.add_argument('--pheno-file',
                        required=True,
                        help='Paths to a file containing phenotype data and participant '
                             'eIDs in CSV format. Headers should contain at least '
                             '[eid,<trait_code>,sex,in_ukb_wbu_testing], '
                             'where <trait_code> matches a Genomics plc phenotype definition '
                             '(type "evaluate-prs --list-phenotypes" for supported '
                             'phenotypes). Additionally, this file can include the following columns to enable '
                             f'survival analysis in binary traits: {SUPPORTED_PHENOTYPE_HEADERS}')
    parser.add_argument('--ancestry-file',
                        required=False,
                        default=None,
                        help='Path to a file containing inferred ancestry labels and participant. Headers should be '
                             '[eid,ancestry]')
    parser.add_argument('--pcs-file', required=False, default=None,
                        help='Path to a file containing the first 4 genetically inferred principal components. '
                             'Headers should be [eid,pc1,pc2,pc3,pc4]')
    parser.add_argument('--output-dir', required=False, default='.', help='Output directory for evaluation '
                                                                          'report and CSV containing metrics')
    parser.add_argument('--list-phenotypes', action=ListPhenotypes, help='List supported Genomics plc '
                                                                         'phenotype definitions')
    return parser


def run_from_command_line(argv: list = None):
    args_dict = preprocess_command_line_inputs(argv)

    if not os.path.exists(os.path.join(args_dict['output_dir'], 'plots')):
        os.mkdir(os.path.join(args_dict['output_dir'], 'plots'))

    trait_code = resolve_column_header(args_dict['pheno_file'])
    phenotype_dictionary = load_phenotype_dictionary(trait_code)
    eval_dict, cross_ancestry_eval_dict, meta_dict = dict(), dict(), dict()
    if len(args_dict['prs_files']) != 2:
        raise AssertionError('The UKB PRS Evaluation Tool only support side-by-side analysis of two input PRS.'
                             ' Please only enter two filepaths to PRS files.')
    prs_file1, prs_file2 = args_dict['prs_files']

    df1, df2, meta_dict1, meta_dict2 = read_and_prepare_data(prs_file1, prs_file2, args_dict['pheno_file'], trait_code,
                                                             phenotype_dictionary, args_dict['ancestry_file'],
                                                             args_dict['pcs_file'])

    # Analyse prs_file1
    data_tag = resolve_column_header(prs_file1)
    eval_dict[data_tag], cross_ancestry_eval_dict[data_tag] = evaluate_prs(df1, trait_code, data_tag,
                                                                           phenotype_dictionary,
                                                                           args_dict['output_dir'])
    meta_dict[data_tag] = meta_dict1

    # Analyse prs_file2
    data_tag = resolve_column_header(prs_file2)
    eval_dict[data_tag], cross_ancestry_eval_dict[data_tag] = evaluate_prs(df2, trait_code, data_tag,
                                                                           phenotype_dictionary,
                                                                           args_dict['output_dir'])
    meta_dict[data_tag] = meta_dict2

    results = construct_results_dict(eval_dict, phenotype_dictionary, cross_ancestry_eval_dict, meta_dict, trait_code)
    write_metrics_csv(results, os.path.join(args_dict['output_dir'], 'evaluation_metrics.csv'))
    generate_report(results, os.path.join(args_dict['output_dir'], 'prs_evaluation_report.pdf'))


def construct_results_dict(eval_dict: dict, phenotype_dict: dict, cross_ancestry_eval_dict: dict, meta_dict: dict,
                           trait_code: str):

    results = {'phenotype_code': trait_code,
               'phenotype_description': phenotype_dict['full_name'],
               'flavour': phenotype_dict['scale'],
               'sex_major': phenotype_dict['gender'],
               'evaluation': eval_dict,
               'cross_ancestry_evaluation': cross_ancestry_eval_dict,
               'metadata': meta_dict
               }
    return results


def main(argv: list = None) -> None:
    """Entrypoint for ukb-pret command line usage

    Parameters
    ----------
    argv: list
        list of options collected from command line"""
    run_from_command_line(argv)


def preprocess_command_line_inputs(argv: list = None):
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()
    args_dict = vars(parser.parse_args(argv))
    return args_dict
