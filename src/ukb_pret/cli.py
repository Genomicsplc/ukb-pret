"""
Command line interface for evaluating PRS
"""

import argparse
import os
import pandas
import sys

from ._io import resolve_column_header, resolve_rap_inputs, load_phenotype_dictionary, read_prs_file, read_pheno_file, \
    read_principal_components_file

from .constants import SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS, DF_COLUMN_HEADERS, UNSUPPORTED_DF_CHARS
from .evaluate import compare_prs


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


def build_parser():
    """Defines the arguments taken by the command line parser."""

    parser = argparse.ArgumentParser(description="A CLI tool for evaluating a set of PRS against a phenotype.")

    parser.add_argument('--prs-files',
                        nargs=2,
                        required=True,
                        help='Paths to two files, each containing Polygenic Risk Score (PRS) and '
                             'participant eIDs in CSV format. Headers should be [eid,<data_tag>], '
                             'where <data_tag> is a field without spaces or special characters that is '
                             'used to identify the PRS in the output (REQUIRED)')
    parser.add_argument('--pheno-file',
                        required=True,
                        help='Paths to a file containing phenotype data and participant '
                             'eIDs in CSV format. Headers should contain at least '
                             '[eid,<trait_code>], '
                             'where <trait_code> is a field without spaces or special characters that '
                             'can either correspond to an existing Gplc phenotype '
                             'definition or be defined by the user. [sex] can also be '
                             'included as a header for stratified analysis using coding {0: female, 1: male}. '
                             'Additionally, this file can include the following columns to enable '
                             f'survival analysis in binary traits: '
                             f'{SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS} (REQUIRED)')
    parser.add_argument('--pcs-file', required=False, default=None,
                        help='Path to a file containing the first 4 genetically inferred principal components. '
                             'Headers should be [eid,pc1,pc2,pc3,pc4] (OPTIONAL) (when omitted, evaluation is '
                             'carried out across all ancestries & the report does not contain a '
                             'quality control section)')
    parser.add_argument('--output-dir', required=False, default='.', help='Output directory for evaluation '
                                                                          'report and CSV containing metrics'
                                                                          ' (default is current working '
                                                                          'directory) (OPTIONAL)')
    return parser


def build_parser_rap():

    parser = argparse.ArgumentParser(description="A CLI tool for evaluating a set of PRS against a phenotype on the "
                                                 "Research Analysis Platform.")

    parser.add_argument('--ukb-release-prs-file',
                        required=True,
                        help='Path to a CSV file containing a Polygenic Risk Score (PRS), participant eIDs, sex and '
                             'the first 4 genetically inferred principal components. '
                             'Headers should be [eid,p<dataset_id>,sex,pc1,pc2,pc3,pc4], '
                             'where <dataset_id> is the field used to identify the phenotype used to generate the PRS '
                             '(REQUIRED)')
    parser.add_argument('--user-prs-file', required=True, default=None,
                        help='Path to a CSV file containing a Polygenic Risk Score (PRS) and participant eIDs '
                             'provided by the user. Headers should be [eid,<data_tag>], where <data_tag> is a field '
                             'without spaces or special characters that is '
                             'used to identify the PRS in the output (REQUIRED)')
    parser.add_argument('--pheno-file',
                        required=True,
                        help='Paths to a file containing phenotype data and participant '
                             'eIDs in CSV format. Headers should contain at least '
                             '[eid,<trait_code>], where <trait_code> is a field without spaces or special characters '
                             'that can correspond to an existing Gplc phenotype '
                             'definition or be defined by the user. Additionally, this file can include the following '
                             f'columns to enable survival analysis in binary traits:'
                             f' {SUPPORTED_PHENOTYPE_HEADERS_SURVIVAL_ANALYSIS} (REQUIRED)')
    parser.add_argument('--output-dir', required=False, default='.',
                        help='Output directory for evaluation report and CSV containing metrics '
                             '(default is current working directory) (OPTIONAL)')
    return parser


def validate_inputs(metadata: dict):

    for k, v in metadata.items():
        if k in DF_COLUMN_HEADERS and any(x in v for x in UNSUPPORTED_DF_CHARS):
            print(f'The input {v} is invalid - please ensure it does not contain whitespace or invalid characters')
            sys.exit(1)

    if not os.path.exists(os.path.join(metadata['output_dir'], 'plots')):
        os.makedirs(os.path.join(metadata['output_dir'], 'plots'))


def run_from_command_line(argv: list = None) -> None:
    """
    Entrypoint for ukb-pret command line usage

    Parameters
    ----------
    argv: list
        list of options collected from command line
    """

    args_dict = preprocess_command_line_inputs(build_parser(), argv)
    trait_code = resolve_column_header(args_dict['pheno_file'])

    phenotype_dictionary = load_phenotype_dictionary()
    if trait_code in phenotype_dictionary.keys():
        phenotype_dictionary = phenotype_dictionary[trait_code]
    else:
        print('ukb-pret does not support custom trait codes via the evaluate-prs entrypoint\n'
              'The input phenotype file must use a valid Gplc trait code to identify the phenotype data\n'
              f'Found {trait_code} in {args_dict["pheno_file"]} instead.')
        sys.exit(1)

    prs1_filepath, prs2_filepath = args_dict['prs_files']

    prs1_label = resolve_column_header(prs1_filepath)
    prs1_df = read_prs_file(prs1_filepath, prs1_label)
    prs2_label = resolve_column_header(prs2_filepath)
    prs2_df = read_prs_file(prs2_filepath, prs2_label)
    pheno_df = read_pheno_file(args_dict['pheno_file'], trait_code, phenotype_dictionary)
    pcs_df = read_principal_components_file(args_dict['pcs_file']).dropna()

    ukb_pret(trait_code, prs1_label, prs2_label, args_dict['output_dir'], phenotype_dictionary,
             prs1_df, prs2_df, pheno_df, pcs_df)


def run_from_command_line_rap(argv: list = None) -> None:
    """Entrypoint for ukb-pret for use with the DNANexus Research Analysis Platform (RAP)

    Parameters
    ----------
    argv: list
        list of options collected from command line
    """

    args_dict = preprocess_command_line_inputs(build_parser_rap(), argv)
    prs1_label = resolve_column_header(args_dict['user_prs_file'])
    prs1_df = read_prs_file(args_dict['user_prs_file'], prs1_label)
    prs2_df, pcs_df, sex_df, prs2_label, gplc_trait_code = resolve_rap_inputs(args_dict['ukb_release_prs_file'])
    trait_code = resolve_column_header(args_dict['pheno_file'])
    phenotype_dictionary = load_phenotype_dictionary()
    if trait_code in phenotype_dictionary.keys():
        phenotype_dictionary = phenotype_dictionary[trait_code]
    else:
        # For a user-defined trait code that doesn't directly correspond to a Gplc definition,
        # inherit metadata from the Gplc phenotype definition
        phenotype_dictionary = phenotype_dictionary[gplc_trait_code]
        phenotype_dictionary['full_name'] = 'User-defined phenotype'

    pheno_df = read_pheno_file(args_dict['pheno_file'], trait_code, phenotype_dictionary)
    if 'sex' in pheno_df.columns:
        print('WARNING: sex data are provided via the UKB PRS Release and the phenotype file.\n'
              'Using phenotype sex data...')
    else:
        pheno_df = pheno_df.join(sex_df)

    ukb_pret(trait_code, prs1_label, prs2_label, args_dict['output_dir'], phenotype_dictionary,
             prs1_df, prs2_df, pheno_df, pcs_df)


def ukb_pret(trait_code: str, prs1_label: str, prs2_label: str, output_dir: str, phenotype_dictionary: dict,
             prs1_data: pandas.DataFrame, prs2_data: pandas.DataFrame, pheno_data: pandas.DataFrame,
             pc_data: pandas.DataFrame):
    metadata = {'prs1_label': prs1_label,
                'prs2_label': prs2_label,
                'trait_code': trait_code,
                'output_dir': output_dir,
                'phenotype_dictionary': phenotype_dictionary}
    validate_inputs(metadata)
    compare_prs(prs1_data, prs2_data, pheno_data, metadata, pc_data)


def preprocess_command_line_inputs(parser, argv: list = None):
    if argv is None:
        argv = sys.argv[1:]
    args_dict = vars(parser.parse_args(argv))
    return args_dict
