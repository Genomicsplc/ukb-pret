import numpy
import pandas
import tempfile
import unittest
from unittest.mock import patch

from ukb_pret.cli import preprocess_command_line_inputs, build_parser, validate_inputs, run_from_command_line, \
    run_from_command_line_rap
from ukb_pret.constants import SURVIVAL_DATA_HEADERS, SUPPORTED_PC_HEADERS


UKB_PRET_CLI_PATCH_PATH = 'ukb_pret.cli'


def generate_basic_metadata():
    tmp_dir = tempfile.TemporaryDirectory()
    metadata = {'prs1_label': 'valid_label',
                'prs2_label': 'other_valid_label',
                'trait_code': 'ATT',
                'output_dir': tmp_dir.name,
                'phenotype_dictionary': {}}
    return metadata


class ExpectedUkbPretInputs:
    """
    A class containing the standard expected inputs to the ukb_pret function for testing
    The outputs of this class may be tweaked inside specific tests to match the use case
    """

    ast_pheno_dict = {'follow_up_months': '120', 'full_name': 'Asthma', 'sex': 'both', 'scale': 'binary'}
    expected_prs1_df = pandas.DataFrame.from_dict({'user_prs': {'1001': 0.5, '1002': -0.3, '1003': 1.02}})
    expected_prs2_df = pandas.DataFrame.from_dict({'UKB_Enhanced_PRS_AST': {'1001': 0.2, '1002': -0.3, '1003': 0.4}})
    expected_pheno_df = pandas.DataFrame.from_dict({'AST': {'1001': 0, '1002': 1, '1003': 0},
                                                    'incident': {'1001': 0.0, '1002': None, '1003': 0.0},
                                                    'age_at_first_assessment': {
                                                        '1001': 47.04, '1002': 65.07, '1003': 50.3
                                                    },
                                                    'date_of_first_assessment': {
                                                        '1001': pandas.Timestamp('2009-08-07 00:00:00'),
                                                        '1002': pandas.Timestamp('2010-02-03 00:00:00'),
                                                        '1003': pandas.Timestamp('2010-02-04 00:00:00')},
                                                    'date_of_diagnosis': {
                                                        '1001': pandas.NaT,
                                                        '1002': pandas.Timestamp('2011-02-03 00:00:00'),
                                                        '1003': pandas.NaT
                                                    },
                                                    'date_of_death': {
                                                        '1001': pandas.NaT, '1002': pandas.NaT, '1003': pandas.NaT
                                                    },
                                                    'prevalent': {'1001': 0.0, '1002': 1.0, '1003': 0.0},
                                                    'censoring_date': {
                                                        '1001': pandas.Timestamp('2019-08-05 00:00:00'),
                                                        '1002': pandas.Timestamp('2020-02-01 00:00:00'),
                                                        '1003': pandas.Timestamp('2020-02-02 00:00:00')
                                                    },
                                                    'censoring_date_follow_up': {
                                                        '1001': pandas.Timestamp('2019-08-05 00:00:00'),
                                                        '1002': pandas.Timestamp('2020-02-01 00:00:00'),
                                                        '1003': pandas.Timestamp('2020-02-02 00:00:00')
                                                    },
                                                    'age_at_diagnosis': {'1001': numpy.nan, '1002': 66.07,
                                                                         '1003': numpy.nan},
                                                    'date_of_last_encounter': {
                                                        '1001': pandas.Timestamp('2019-08-05 00:00:00'),
                                                        '1002': pandas.Timestamp('2020-02-01 00:00:00'),
                                                        '1003': pandas.Timestamp('2020-02-02 00:00:00')
                                                    },
                                                    'age_at_censoring': {'1001': 57.04, '1002': 75.07, '1003': 60.3},
                                                    'age_event': {'1001': 57.04, '1002': 66.07, '1003': 60.3},
                                                    'time_survival': {'1001': 10.0, '1002': numpy.nan, '1003': 10.0}
                                                    })
    expected_pheno_df_with_sex = expected_pheno_df.join(pandas.DataFrame.from_dict({'eid': ['1001', '1002', '1003'],
                                                                                    'sex': [0, 1, 1]}).set_index('eid'))

    expected_pcs_df = pandas.DataFrame.from_dict({
        'pc1': {'1001': 0.002, '1002': 0.3453, '1003': -1.257},
        'pc2': {'1001': 0.32547, '1002': 0.1234, '1003': -1.3456},
        'pc3': {'1001': -0.3567, '1002': 0.2345, '1003': -1.5467},
        'pc4': {'1001': 0.245, '1002': 0.2346, '1003': -1.2346}
    })

    expected_dfs = [expected_prs1_df, expected_prs2_df, expected_pheno_df_with_sex, expected_pcs_df]
    for df in expected_dfs:
        df.index.name = 'eid'


class TestCommandLineInterface(unittest.TestCase):

    def test_no_prs_files(self):
        input = ['--pheno-file', 'a/b/c.csv']
        with self.assertRaises(SystemExit) as se:
            preprocess_command_line_inputs(build_parser(), input)

    def test_no_pheno_file(self):
        input = ['--prs-files', 'a/b/c.csv']
        with self.assertRaises(SystemExit) as se:
            preprocess_command_line_inputs(build_parser(), input)

    def test_multiple_prs_files(self):
        input = ['--prs-files', 'a/b/c.csv', '/d/e/f.csv', '--pheno-file', 'test.pheno']
        preprocess_command_line_inputs(build_parser(), input)

    def test_validate_inputs(self):
        md = generate_basic_metadata()
        validate_inputs(md)

    def test_validate_inputs_fails_with_whitespace(self):
        md = generate_basic_metadata()
        md['prs2_label'] = 'invalid due to whitespace'
        with self.assertRaises(SystemExit):
            validate_inputs(md)

    def test_validate_inputs_fails_with_special_character(self):
        md = generate_basic_metadata()
        md['prs2_label'] = 'a+crazy?trait:'
        with self.assertRaises(SystemExit):
            validate_inputs(md)


def mock_user_prs_df():
    return pandas.DataFrame({'eid': ['1001', '1002', '1003'], 'user_prs': [0.5, -0.3, 1.02]})


def mock_alt_user_prs_df():
    return pandas.DataFrame({'eid': ['1001', '1002', '1003'], 'UKB_Enhanced_PRS_AST': [0.2, -0.3, 0.4]})


def mock_ukb_release_prs_file_df():
    return pandas.DataFrame({
        'eid': ['1001', '1002', '1003'],
        'p26211': [0.2, -0.3, 0.4],
        'pc1': [0.002, 0.3453, -1.257],
        'pc2': [0.32547, 0.1234, -1.3456],
        'pc3': [-0.3567, 0.2345, -1.5467],
        'pc4': [0.245, 0.2346, -1.2346],
        'sex': [0, 1, 1]
    })


def mock_pheno_file(trait_code: str = 'AST'):
    return pandas.DataFrame({
        'eid': ['1001', '1002', '1003'],
        trait_code: [0, 1, 0],
        'incident': [0, 0, 0],
        'age_at_first_assessment': [47.04, 65.07, 50.3],
        'date_of_first_assessment': ['2009-08-07', '2010-02-03', '2010-02-04'],
        'date_of_diagnosis': ['', '2011-02-03', ''],
        'date_of_death': ['', '', '']
    })


def mock_pheno_file_with_sex(trait_code: str = 'AST'):
    pheno_df = mock_pheno_file(trait_code)
    pheno_df['sex'] = [0, 1, 1]
    return pheno_df


def mock_pcs_file():
    return pandas.DataFrame({
        'eid': ['1001', '1002', '1003'],
        'pc1': [0.002, 0.3453, -1.257],
        'pc2': [0.32547, 0.1234, -1.3456],
        'pc3': [-0.3567, 0.2345, -1.5467],
        'pc4': [0.245, 0.2346, -1.2346]
    })


class TestEvaluatePrsEntrypoint(unittest.TestCase):

    args = ['--prs-files', '/fake/path/1', '/fake/path/2', '--pheno-file', '/third/path']
    expected = ExpectedUkbPretInputs()

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_gplc_trait_code(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_file_with_sex(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_file_with_sex(), mock_pcs_file()]
            args = self.args + ['--pcs-file', '/a/pcs/file/path']
            run_from_command_line(args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), self.expected.expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_custom_trait_code(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_file_with_sex('ATT'), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_file_with_sex('ATT'), mock_pcs_file()]
            args = self.args + ['--pcs-file', '/a/pcs/file/path']
            with self.assertRaises(SystemExit):
                run_from_command_line(args)

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_sex(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_file(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_file(), mock_pcs_file()]
            args = self.args + ['--pcs-file', '/a/pcs/file/path']
            run_from_command_line(args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected = [self.expected.expected_prs1_df, self.expected.expected_prs2_df, self.expected.expected_pheno_df,
                    self.expected.expected_pcs_df]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_survival_covariates(self, mock_ukb_pret):
        mock_pheno_df = mock_pheno_file_with_sex()[['eid', 'AST', 'sex']]
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_df.copy(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_df.copy(), mock_pcs_file()]
            args = self.args + ['--pcs-file', '/a/pcs/file/path']
            run_from_command_line(args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected = [self.expected.expected_prs1_df, self.expected.expected_prs2_df,
                    self.expected.expected_pheno_df_with_sex[['AST', 'sex']], self.expected.expected_pcs_df]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_survival_covariates_and_sex(self, mock_ukb_pret):
        mock_pheno_df = mock_pheno_file()[['eid', 'AST']]
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_df.copy(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_df.copy(), mock_pcs_file()]
            args = self.args + ['--pcs-file', '/a/pcs/file/path']
            run_from_command_line(args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected = [self.expected.expected_prs1_df, self.expected.expected_prs2_df,
                    self.expected.expected_pheno_df_with_sex[['AST']], self.expected.expected_pcs_df]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_principal_components(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_file_with_sex(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_file_with_sex(), mock_pcs_file()]
            run_from_command_line(self.args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected_dfs = self.expected.expected_dfs[:-1] + [pandas.DataFrame()]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_principal_components_and_sex(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_file(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_file(), mock_pcs_file()]
            run_from_command_line(self.args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected = [self.expected.expected_prs1_df, self.expected.expected_prs2_df, self.expected.expected_pheno_df,
                    pandas.DataFrame()]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_principal_components_and_sex_and_survival_variables(self, mock_ukb_pret):
        mock_pheno_df = mock_pheno_file()[['eid', 'AST']]
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_pheno_df.copy(), mock_user_prs_df(), mock_user_prs_df(),
                                       mock_alt_user_prs_df(), mock_alt_user_prs_df(),
                                       mock_pheno_df.copy(), mock_pcs_file()]
            run_from_command_line(self.args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        expected = [self.expected.expected_prs1_df, self.expected.expected_prs2_df,
                    self.expected.expected_pheno_df_with_sex[['AST']], pandas.DataFrame()]
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected[i].sort_index(axis=1))


def mock_ukb_enhanced_trait_codes():
    return pandas.DataFrame({'gplc_trait_code': ['AST'], 'ukb_enhanced_prs_code': [26211],
                             'ukb_field_description': 'Description of Asthma'})


class TestEvaluatePrsRapEntrypoint(unittest.TestCase):

    args = ['--ukb-release-prs-file', '/fake/path', '--user-prs-file', '/other/path',
            '--pheno-file', '/third/path']

    expected = ExpectedUkbPretInputs()

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_custom_trait_code(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), mock_ukb_release_prs_file_df(),
                                       mock_ukb_release_prs_file_df(), mock_pheno_file(), mock_pheno_file()]
            run_from_command_line_rap(self.args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), self.expected.expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_gplc_trait_code(self, mock_ukb_pret):
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), mock_ukb_release_prs_file_df(),
                                       mock_ukb_release_prs_file_df(), mock_pheno_file(trait_code='ATT'),
                                       mock_pheno_file(trait_code='ATT')]
            run_from_command_line_rap(self.args)
        pheno_dict = self.expected.ast_pheno_dict.copy()
        pheno_dict['full_name'] = 'User-defined phenotype'
        simple_args = ('ATT', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', pheno_dict)
        expected_pheno_df = self.expected.expected_pheno_df_with_sex.copy()
        expected_pheno_df = expected_pheno_df.rename(columns={'AST': 'ATT'})
        expected_dfs = [self.expected.expected_prs1_df, self.expected.expected_prs2_df, expected_pheno_df,
                        self.expected.expected_pcs_df]

        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_sex(self, mock_ukb_pret):
        # NB: Sex should always be returned in the RAP so this situation should not occur there
        mock_ukb_release_prs_file_df_without_sex = mock_ukb_release_prs_file_df().drop(columns=['sex'])
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), mock_ukb_release_prs_file_df_without_sex.copy(),
                                       mock_ukb_release_prs_file_df_without_sex.copy(),
                                       mock_pheno_file(), mock_pheno_file()]
            with self.assertRaises(SystemExit):
                run_from_command_line_rap(self.args)

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_survival_covariates(self, mock_ukb_pret):
        mock_pheno_file_df = mock_pheno_file()
        mock_pheno_file_df = mock_pheno_file_df[[c for c in mock_pheno_file_df if c not in SURVIVAL_DATA_HEADERS]]
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), mock_ukb_release_prs_file_df(),
                                       mock_ukb_release_prs_file_df(), mock_pheno_file_df,
                                       mock_pheno_file_df]
            run_from_command_line_rap(self.args)

        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        expected_pheno_df = self.expected.expected_pheno_df_with_sex.copy()
        expected_pheno_df = expected_pheno_df[['AST', 'sex']]
        expected_dfs = [self.expected.expected_prs1_df, self.expected.expected_prs2_df, expected_pheno_df,
                        self.expected.expected_pcs_df]

        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_sex_in_phenotype_file_overwrites_sex_in_ukb_prs_output(self, mock_ukb_pret):
        ukb_release_prs_file = mock_ukb_release_prs_file_df()
        ukb_release_prs_file['sex'] = [1, 1, 1]
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), ukb_release_prs_file.copy(),
                                       ukb_release_prs_file.copy(), mock_pheno_file_with_sex(),
                                       mock_pheno_file_with_sex()]
            run_from_command_line_rap(self.args)
        simple_args = ('AST', 'user_prs', 'UKB_Enhanced_PRS_AST', '.', self.expected.ast_pheno_dict)
        self.assertEqual(simple_args, mock_ukb_pret.call_args[0][:5])
        for i, df in enumerate(mock_ukb_pret.call_args[0][5:]):
            pandas.testing.assert_frame_equal(df.sort_index(axis=1), self.expected.expected_dfs[i].sort_index(axis=1))

    @patch(f'{UKB_PRET_CLI_PATCH_PATH}.ukb_pret')
    def test_without_principal_components(self, mock_ukb_pret):
        # NB: Sex should always be returned in the RAP so this situation should not occur there
        mock_ukb_release_prs_file_df_without_sex = mock_ukb_release_prs_file_df().drop(columns=SUPPORTED_PC_HEADERS)
        with patch(f'{UKB_PRET_CLI_PATCH_PATH}.pandas.read_csv') as mock_pandas:
            mock_pandas.side_effect = [mock_user_prs_df(), mock_user_prs_df(),
                                       mock_ukb_enhanced_trait_codes(), mock_ukb_release_prs_file_df_without_sex.copy(),
                                       mock_ukb_release_prs_file_df_without_sex.copy(),
                                       mock_pheno_file(), mock_pheno_file()]
            with self.assertRaises(SystemExit):
                run_from_command_line_rap(self.args)
