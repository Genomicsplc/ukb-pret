from io import StringIO
import numpy
import pandas
import unittest
from unittest.mock import patch
import yaml

from ukb_pret._io import read_principal_components_file, read_pop_cluster_centers, read_prs_file, read_pheno_file,\
    read_ukb_prs_file
from ukb_pret.constants import PHENOTYPE_TRAIT_DEFS_PATH, POPULATION_CLUSTER_CENTRES_PATH


class TestInputs(unittest.TestCase):

    def test_traits_yaml_keys(self):
        with open(PHENOTYPE_TRAIT_DEFS_PATH, 'r') as f:
            traits_dict = yaml.safe_load(f)
        for trait_dict in traits_dict.values():
            self.assertTrue(all(x in trait_dict.keys() for x in ['follow_up_months', 'full_name', 'sex', 'scale']))

    def test_read_pop_cluster_centers(self):
        df = read_pop_cluster_centers(POPULATION_CLUSTER_CENTRES_PATH)
        self.assertCountEqual(set(df.columns), ('pc1', 'pc2', 'pc3', 'pc4'))
        self.assertEqual(str(df.index.name), 'population')

    @patch('builtins.open')
    def test_reading_valid_pc_file(self, mock_open):
        mock_open.return_value = StringIO("eid,pc1,pc2,pc3,pc4\nFK1,0.1,0.2,0.3,0.4\nFK2,0.6,0.7,0.8,0.9")
        df = read_principal_components_file('/a/fake/path')
        self.assertCountEqual(set(df.columns), ('pc1', 'pc2', 'pc3', 'pc4'))
        self.assertEqual(str(df.index.name), 'eid')

    @patch('builtins.open')
    def test_reading_pc_file_missing_columns(self, mock_open):
        mock_open.return_value = StringIO("eid,pc1,pc2,pc4\nFK1,0.1,0.2,0.4\nFK2,0.6,0.7,0.9")
        with self.assertRaises(SystemExit):
            _ = read_principal_components_file('/a/fake/path')

    @patch('builtins.open')
    def test_valid_prs_file(self, mock_open):
        mock_open.return_value = StringIO("eid,my_prs\nFK1,0.5\nFK2,-0.1")
        df = read_prs_file('/a/fake/path', 'my_prs')
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'], 'my_prs':  [0.5, -0.1]})

    @patch('builtins.open')
    def test_missing_header_prs_file(self, mock_open):
        mock_open.return_value = StringIO("eid\nFK1\nFK2")
        with self.assertRaises(SystemExit):
            _ = read_prs_file('/a/fake/path', 'my_prs')

    @patch('builtins.open')
    def test_extra_header_prs_file(self, mock_open):
        mock_open.return_value = StringIO("eid,my_prs,other_data\nFK1,0.5,12376\nFK2,-0.1,4368")
        with self.assertRaises(SystemExit):
            _ = read_prs_file('/a/fake/path', 'my_prs')

    @patch('builtins.open')
    def test_valid_binary_pheno_file(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno\nFK1,0\nFK2,1")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'other_stuff': ''})
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'], 'my_pheno':  [0, 1]})

    @patch('builtins.open')
    def test_valid_quantitative_pheno_file(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno\nFK1,0.34\nFK2,-9.234")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'other_stuff': ''})
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'], 'my_pheno':  [0.34, -9.234]})

    @patch('builtins.open')
    def test_valid_binary_pheno_file_with_sex(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,sex\nFK1,0,0\nFK2,1,0")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'other_stuff': ''})
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'], 'my_pheno':  [0, 1], 'sex': [0, 0]})

    @patch('builtins.open')
    def test_valid_binary_pheno_file_without_sex(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno\nFK1,0\nFK2,1")
        _ = read_pheno_file('/a/fake/path', 'my_pheno', {'other_stuff': ''})

    @patch('builtins.open')
    def test_binary_pheno_file_with_incompatible_header_subset(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,age_at_first_assessment,date_of_diagnosis,"
                                          "date_of_first_assessment,incident\n"
                                          "FK1,0,42.01,,2007-08-08,0\n"
                                          "FK2,1,51.82,2009-07-07,2010-06-06,0")
        with self.assertRaises(SystemExit):
            _ = read_pheno_file('/a/fake/path', 'my_pheno', {'': ''})

    @patch('builtins.open')
    def test_binary_pheno_file_with_extra_headers(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,age_at_first_assessment,date_of_diagnosis,"
                                          "date_of_first_assessment,date_of_death,incident,extra_header\n"
                                          "FK1,0,42.01,,2007-08-08,,0,123445\n"
                                          "FK2,1,51.82,2009-07-07,2010-06-06,,0,asdf")
        with self.assertRaises(SystemExit):
            _ = read_pheno_file('/a/fake/path', 'my_pheno', {'': ''})

    @patch('builtins.open')
    def test_valid_quantitative_pheno_file_with_sex(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,sex\nFK1,0.45,1\nFK2,1.123,1")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'other_stuff': ''})
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'], 'my_pheno':  [0.45, 1.123],
                                                            'sex': [1, 1]})

    @patch('builtins.open')
    def test_valid_ukb_prs_file(self, mock_open):
        mock_open.return_value = StringIO("eid,sex,p26211,pc1,pc2,pc3,pc4\nFK1,1,1,0.1,0.2,0.3,0.4\n"
                                          "FK2,1,0,-0.1,-0.2,-0.3,-0.4")
        df = read_ukb_prs_file('/a/fake/path', 26211)
        self.assertEqual(df.reset_index().to_dict('list'), {'eid': ['FK1', 'FK2'],
                                                            'p26211': [1, 0],
                                                            'pc1': [0.1, -0.1],
                                                            'pc2': [0.2, -0.2],
                                                            'pc3': [0.3, -0.3],
                                                            'pc4': [0.4, -0.4],
                                                            'sex': [1, 1]})

    @patch('builtins.open')
    def test_ukb_prs_file_with_extra_column(self, mock_open):
        mock_open.return_value = StringIO("eid,sex,p26211,pc1,pc2,pc3,pc4,other\nFK1,1,1,0.1,0.2,0.3,0.4,23\n"
                                          "FK2,1,0,-0.1,-0.2,-0.3,-0.4,23485")
        with self.assertRaises(SystemExit):
            _ = read_ukb_prs_file('/a/fake/path', 26211)

    @patch('builtins.open')
    def test_valid_binary_pheno_file_with_survival_data(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,age_at_first_assessment,date_of_diagnosis,"
                                          "date_of_first_assessment,date_of_death,incident\n"
                                          "FK1,0,42.01,,2007-08-08,,0\n"
                                          "FK2,1,51.82,2009-07-07,2010-06-06,,0")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'follow_up_months': 120})
        self.assertCountEqual(df.reset_index().to_dict('list'),
                              {'age_at_censoring': [52.01, 61.56246575342466],
                               'age_at_diagnosis': [numpy.nan, 50.904931506849316],
                               'age_at_first_assessment': [42.01, 51.82],
                               'age_event': [52.01, 50.904931506849316],
                               'censoring_date':
                                   [pandas.Timestamp('2017-08-05 00:00:00'),
                                    pandas.Timestamp('2020-03-01 00:00:00')],
                               'censoring_date_follow_up':
                                   [pandas.Timestamp('2017-08-05 00:00:00'),
                                    pandas.Timestamp('2020-06-03 00:00:00')],
                               'date_of_death': [pandas.NaT, pandas.NaT],
                               'date_of_diagnosis': [
                                   pandas.NaT, pandas.Timestamp('2009-07-07 00:00:00')],
                               'date_of_first_assessment': [
                                   pandas.Timestamp('2007-08-08 00:00:00'),
                                   pandas.Timestamp('2010-06-06 00:00:00')],
                               'date_of_last_encounter': [
                                   pandas.Timestamp('2017-08-05 00:00:00'),
                                   pandas.Timestamp('2020-03-01 00:00:00')],
                               'eid': ['FK1', 'FK2'],
                               'incident': [0.0, numpy.nan],
                               'my_pheno': [0, 1],
                               'prevalent': [0.0, 1.0],
                               'time_survival': [10.0, numpy.nan]})

    @patch('builtins.open')
    def test_valid_binary_pheno_file_with_survival_data_with_sex(self, mock_open):
        mock_open.return_value = StringIO("eid,my_pheno,sex,age_at_first_assessment,date_of_diagnosis,"
                                          "date_of_first_assessment,date_of_death,incident\n"
                                          "FK1,0,0,42.01,,2007-08-08,,0\n"
                                          "FK2,1,0,51.82,2009-07-07,2010-06-06,,0")
        df = read_pheno_file('/a/fake/path', 'my_pheno', {'follow_up_months': 120})
        self.assertCountEqual(df.reset_index().to_dict('list'),
                              {'age_at_censoring': [52.01, 61.56246575342466],
                               'age_at_diagnosis': [numpy.nan, 50.904931506849316],
                               'age_at_first_assessment': [42.01, 51.82],
                               'age_event': [52.01, 50.904931506849316],
                               'censoring_date':
                                   [pandas.Timestamp('2017-08-05 00:00:00'),
                                    pandas.Timestamp('2020-03-01 00:00:00')],
                               'censoring_date_follow_up':
                                   [pandas.Timestamp('2017-08-05 00:00:00'),
                                    pandas.Timestamp('2020-06-03 00:00:00')],
                               'date_of_death': [pandas.NaT, pandas.NaT],
                               'date_of_diagnosis': [
                                   pandas.NaT, pandas.Timestamp('2009-07-07 00:00:00')],
                               'date_of_first_assessment': [
                                   pandas.Timestamp('2007-08-08 00:00:00'),
                                   pandas.Timestamp('2010-06-06 00:00:00')],
                               'date_of_last_encounter': [
                                   pandas.Timestamp('2017-08-05 00:00:00'),
                                   pandas.Timestamp('2020-03-01 00:00:00')],
                               'eid': ['FK1', 'FK2'],
                               'incident': [0.0, numpy.nan],
                               'my_pheno': [0, 1],
                               'prevalent': [0.0, 1.0],
                               'sex': [0, 0],
                               'time_survival': [10.0, numpy.nan]})
