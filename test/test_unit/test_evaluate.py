import numpy
import os
import pandas
from tempfile import TemporaryDirectory
import unittest

from ukb_pret._io import load_phenotype_dictionary
from ukb_pret.error import UkbPretImportError
from ukb_pret.evaluate import _get_overlapping_results, _filter_null_values, _infer_ancestry_from_pcs, evaluate_prs


def generate_prs_data():
    prs_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'prs_data': [0.4, -0.87, -2.1]}
    return pandas.DataFrame.from_dict(prs_dict).set_index('eid')


def generate_binary_pheno_data():
    bin_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'AST': [1, 0, 0]}
    return pandas.DataFrame.from_dict(bin_dict).set_index('eid')


def generate_quantitative_pheno_data():
    quant_dict = {'eid': ['FK1', 'FK2', 'FK3'], 'LDL': [50.32, 90, -999999999.12345]}
    return pandas.DataFrame.from_dict(quant_dict).set_index('eid')


def generate_pc_data():
    pc_dict = {'eid': ['FK1', 'FK2', 'FK3'],
               'pc1': [0.5, 0.5, 0.5],
               'pc2': [-0.3, 0.3, 0.7],
               'pc3': [0.1, 0.3, 0.5],
               'pc4': [0.7, -0.7, 0.0],
               }
    return pandas.DataFrame.from_dict(pc_dict).set_index('eid')


def generate_population_clusters():
    return pandas.DataFrame({'population': ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
                             'pc1': [0.5, 0.4, 0.3, 0.2, 0.1],
                             'pc2': [0.1, 0.2, 0.3, 0.4, 0.5],
                             'pc3': [-0.5, -0.4, -0.3, -0.2, -0.1],
                             'pc4': [-0.1, -0.2, -0.3, -0.4, -0.5]}).set_index('population')


def generate_sex_and_ukb_testing_data():
    pc_dict = {'eid': ['FK1', 'FK2', 'FK3'],
               'sex': [0, 0, 1],
               'in_ukb_wbu_testing': [0, 1, 1]
               }
    return pandas.DataFrame.from_dict(pc_dict).set_index('eid')


def generate_custom_df(*args):
    return pandas.concat([x() for x in args], axis=1)


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

    def test_infer_ancestry_from_pcs(self):
        df = generate_custom_df(generate_prs_data, generate_binary_pheno_data, generate_pc_data)
        pop_clusters = generate_population_clusters()
        new_df = _infer_ancestry_from_pcs(df, pop_clusters)
        expected_ancestries = {'FK1': 'AFR', 'FK2': 'SAS', 'FK3': 'SAS'}
        self.assertDictEqual(new_df['ancestry'].to_dict(), expected_ancestries)

    def test_evaluate_prs(self):
        df = generate_custom_df(generate_prs_data, generate_quantitative_pheno_data, generate_pc_data,
                                generate_sex_and_ukb_testing_data)
        pop_clusters = generate_population_clusters()
        new_df = _infer_ancestry_from_pcs(df, pop_clusters)

        # Increasing the size of the df to allow prs binning by PCs
        larger_df = pandas.concat([new_df]*40, ignore_index=True)
        larger_df['ancestry'] = ['AFR', 'EAS', 'EUR', 'SAS'] * 30
        larger_df.index = larger_df.index.astype(str)
        traits_yaml = load_phenotype_dictionary('LDL_SF')
        tmp_dir = TemporaryDirectory()
        os.mkdir(os.path.join(tmp_dir.name, 'plots'))
        eval_dict, cross_ancestry_eval_dict = evaluate_prs(larger_df, 'LDL_SF', 'prs_data', traits_yaml, tmp_dir.name)

        self.assertTrue(os.path.isdir(os.path.join(tmp_dir.name, 'plots')))
        self.assertTrue(os.path.isdir(os.path.join(tmp_dir.name, 'plots', 'prs_data_cross_ancestry')))
        expected_files = ['prs_data_prs_hist.png', 'prs_data_prs_box_plot.png'] + \
                         [f'prs_data_pc{i}_by_ancestry.png' for i in range(1, 5)]
        self.assertTrue(all(x in expected_files for x in os.listdir(os.path.join(tmp_dir.name, 'plots',
                                                                                 'prs_data_cross_ancestry'))))
        self.assertIsNotNone(eval_dict)
        self.assertIsNotNone(cross_ancestry_eval_dict)
