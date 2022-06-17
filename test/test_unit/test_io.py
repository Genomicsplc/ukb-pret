from io import StringIO
import unittest
from unittest.mock import patch

from ukb_pret._io import read_principal_components_file


class TestInputs(unittest.TestCase):

    @patch('builtins.open')
    def test_reading_valid_pc_file(self, mock_open):
        mock_open.return_value = StringIO("eid,pc1,pc2,pc3,pc4\nFK1,0.1,0.2,0.3,0.4\nFK2,0.6,0.7,0.8,0.9")
        df = read_principal_components_file('/a/fake/path')
        self.assertCountEqual(set(df.columns), ('pc1', 'pc2', 'pc3', 'pc4'))

    @patch('builtins.open')
    def test_reading_pc_file_missing_columns(self, mock_open):
        mock_open.return_value = StringIO("eid,pc1,pc2,pc4\nFK1,0.1,0.2,0.4\nFK2,0.6,0.7,0.9")
        with self.assertRaises(KeyError):
            df = read_principal_components_file('/a/fake/path')
