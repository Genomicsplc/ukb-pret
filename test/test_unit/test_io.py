from io import StringIO
import unittest
from unittest.mock import patch
import yaml

from ukb_pret._io import read_principal_components_file, read_pop_cluster_centers
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
        with self.assertRaises(KeyError):
            df = read_principal_components_file('/a/fake/path')
