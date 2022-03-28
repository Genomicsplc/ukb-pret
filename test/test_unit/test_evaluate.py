import numpy
import pandas
import unittest

from ukb_pret.error import UkbPretImportError
from ukb_pret.evaluate import _get_overlapping_results, _filter_null_values


def generate_prs_data():
    prs_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'prs_data': [0.4, -0.87, -2.1]}
    return pandas.DataFrame.from_dict(prs_dict).set_index('eid')


def generate_binary_pheno_data():
    bin_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'AST': [1, 0, 0]}
    return pandas.DataFrame.from_dict(bin_dict).set_index('eid')


def generate_quantitative_pheno_data():
    quant_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'LDL': [50.32, 90, -999999999.12345]}
    return pandas.DataFrame.from_dict(quant_dict).set_index('eid')


class TestEvaluate(unittest.TestCase):

    def test_overlapping_ids(self):
        prs_df = generate_prs_data()
        pheno_df = generate_binary_pheno_data()
        pheno_df = pheno_df.append(pandas.DataFrame.from_dict({'eid': ['FK4'], 'AST': [0]}).set_index('eid'))
        prs_df_intersect, pheno_df_intersect = _get_overlapping_results(prs_df, pheno_df)
        self.assertDictEqual(prs_df_intersect.to_dict(), generate_prs_data().to_dict())
        self.assertDictEqual(pheno_df_intersect.to_dict(), generate_binary_pheno_data().to_dict())

    def test_duplicate_column_names(self):
        prs_df = generate_prs_data()
        pheno_df = generate_binary_pheno_data()
        prs_df = prs_df.rename(columns={'prs_data': 'AST'})
        with self.assertRaises(UkbPretImportError):
            _, _ = _get_overlapping_results(prs_df, pheno_df)

    def test_removal_of_nulls(self):
        prs_df = generate_prs_data()
        prs_df = prs_df.append(pandas.DataFrame.from_dict({'eid': ['FK4', 'FK5'],
                                                           'prs_data': [None, numpy.nan]}).set_index('eid'))
        pheno_df = generate_binary_pheno_data()
        pheno_df = pheno_df.append(pandas.DataFrame.from_dict({'eid': ['FK4', 'FK5'],
                                                               'AST': [numpy.nan, None]}).set_index('eid'))

        out_prs, out_pheno, n_prs, n_pheno = _filter_null_values(prs_df, pheno_df, 'AST')
        self.assertEqual(n_prs, 2)
        self.assertEqual(n_pheno, 2)
        self.assertDictEqual(out_prs.to_dict(), generate_prs_data().to_dict())
        self.assertDictEqual(out_pheno.to_dict(), generate_binary_pheno_data().to_dict())
