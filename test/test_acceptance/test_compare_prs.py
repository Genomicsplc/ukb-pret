import numpy
import os
import pandas
from tempfile import TemporaryDirectory
import unittest

from ukb_pret.evaluate import compare_prs
from ukb_pret._io import _add_survival_vars

PHENOTYPE_DICTIONARY = {'ABT': {'follow_up_months': 120, 'full_name': 'A Binary Trait',
                                'sex': 'both', 'scale': 'binary'},
                        'AQT': {'follow_up_months': 120, 'full_name': 'A Quantitative Trait',
                                'sex': 'both', 'scale': 'quantitative'},
                        'ASSBT': {'follow_up_months': 120, 'full_name': 'A Single Sex Binary Trait',
                                  'sex': 'female', 'scale': 'binary'},
                        'ASSQT': {'follow_up_months': 120, 'full_name': 'A Single Sex Quantitative Trait',
                                  'sex': 'female', 'scale': 'quantitative'}
                        }

BINARY_METADATA_DICT = {'prs1_label': 'prs1',
                        'prs2_label': 'prs2',
                        'trait_code': 'ABT',
                        'phenotype_dictionary': PHENOTYPE_DICTIONARY,
                        'output_dir': None}


def create_prs1():
    prs_dict = {'eid': [f'FK{i}' for i in range(1, 1501)], 'prs1': numpy.random.normal(0, 1.1, 1500).tolist()}
    return pandas.DataFrame.from_dict(prs_dict).set_index('eid')


def create_prs2():
    prs_dict = {'eid': [f'FK{i}' for i in range(1, 1501)], 'prs2': numpy.random.normal(0, 0.9, 1500).tolist()}
    return pandas.DataFrame.from_dict(prs_dict).set_index('eid')


def create_base_binary_phenotype(trait_code: str):
    bin_dict = {'eid': [f'FK{i}' for i in range(1, 1501)], trait_code: [1, 0, 0] * 500}
    return pandas.DataFrame.from_dict(bin_dict).set_index('eid')


def create_metadata_dict(trait_code: str):
    return {'prs1_label': 'prs1',
            'prs2_label': 'prs2',
            'trait_code': trait_code,
            'phenotype_dictionary': PHENOTYPE_DICTIONARY[trait_code]}


def create_base_quantitative_phenotype(trait_code: str):
    quant_dict = {'eid': [f'FK{i}' for i in range(1, 1501)], trait_code: [0.34, 2.4, 1.34] * 500}
    return pandas.DataFrame.from_dict(quant_dict).set_index('eid')


def create_pcs():
    pcs_dict = {'eid': [f'FK{i}' for i in range(1, 1501)],
                'pc1': numpy.random.normal(0, 0.9, 1500),
                'pc2': numpy.random.normal(0, 0.9, 1500),
                'pc3': numpy.random.normal(0, 0.9, 1500),
                'pc4': numpy.random.normal(0, 0.9, 1500)
                }
    return pandas.DataFrame.from_dict(pcs_dict).set_index('eid')


def add_both_sexes_to_phenotype(pheno_df: pandas.DataFrame):
    pheno_df['sex'] = [1, 0, 0, 0, 1, 1, 1, 0, 0, 0] * 150
    return pheno_df


def add_single_sex_to_phenotype(pheno_df: pandas.DataFrame, sex: int):
    pheno_df['sex'] = sex
    return pheno_df


def add_survival_data_to_phenotype(pheno_df: pandas.DataFrame, trait_code: str):
    survival_dict = {'eid': [f'FK{i}' for i in range(1, 1501)],
                     'age_at_first_assessment': [58.23, 34.54, 56.43, 48.9, 30.2, 43.2, 54.3, 23.4, 78.5, 34.4] * 150,
                     'date_of_first_assessment': [pandas.Timestamp('2009-08-07 00:00:00'),
                                                  pandas.Timestamp('2010-02-03 00:00:00'),
                                                  pandas.Timestamp('2010-02-04 00:00:00')] * 500,
                     'date_of_diagnosis': [pandas.NaT,
                                           pandas.Timestamp('2011-02-03 00:00:00'),
                                           pandas.NaT] * 500,
                     'date_of_death': [pandas.NaT] * 1500,
                     'incident': [0.] * 1500}

    pheno_df = pheno_df.join(pandas.DataFrame.from_dict(survival_dict).set_index('eid'))
    return _add_survival_vars(pheno_df, trait_code)


class TestComparePrsBinary(unittest.TestCase):
    """
    Simple tests that check the evaluate_prs method finishes successfully for all supported input combinations
    Also checks for the minimum expected output (that evaluation_metrics.csv produces some discrimination metrics,
    and that at least one plot is produced)
    """

    @classmethod
    def setUpClass(cls):
        cls._temp_output_dir = TemporaryDirectory()

    def _build_standard_inputs(self, trait_code: str, test_name: str):
        numpy.random.seed(42)
        metadata = create_metadata_dict(trait_code)
        output_dir = os.path.join(self._temp_output_dir.name, test_name)
        os.mkdir(output_dir)
        os.mkdir(os.path.join(output_dir, 'plots'))
        metadata['output_dir'] = output_dir
        return create_prs1(), create_prs2(), create_base_binary_phenotype(trait_code), metadata

    def _expected_plots_no_pcs_no_survival(self, metadata_dict):
        for dir in ['prs1_All', 'prs2_All']:
            self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots', dir, 'roc_plot.png')))

    def _expected_plots_with_pcs_no_survival(self, metadata_dict):
        for label in ['prs1', 'prs2']:
            for anc in ['AFR', 'EAS', 'EUR', 'SAS']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots',
                                                            f'{label}_{anc}', 'roc_plot.png')))

                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_box_plot.png')))
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_hist.png')))
            for pc in ['pc1', 'pc2', 'pc3', 'pc4']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_{pc}_by_ancestry.png')))

    def _expected_plots_no_pcs_with_survival(self, metadata_dict):
        for dir in ['prs1_All', 'prs2_All']:
            self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots', dir, 'roc_plot.png')))
            self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots', dir, 'cum_inc_plot.png')))

    def _expected_plots_with_pcs_with_survival(self, metadata_dict):
        for label in ['prs1', 'prs2']:
            # Only 2 ancestries make it past thresholding for fixed seed
            for anc in ['AFR', 'EAS']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots',
                                                            f'{label}_{anc}', 'roc_plot.png')))
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots',
                                                            f'{label}_{anc}', 'cum_inc_plot.png')))

                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_box_plot.png')))
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_hist.png')))
            for pc in ['pc1', 'pc2', 'pc3', 'pc4']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_{pc}_by_ancestry.png')))

    def test_binary_no_sex_no_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_no_sex_no_pcs_no_survival')
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_no_survival(metadata_dict)

    def test_binary_one_sex_no_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSBT', 'binary_one_sex_no_pcs_no_survival')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_no_survival(metadata_dict)

    def test_binary_both_sexes_no_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_both_sexes_no_pcs_no_survival')
        pheno = add_both_sexes_to_phenotype(pheno)
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_no_survival(metadata_dict)

    def test_binary_no_sex_with_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_no_sex_with_pcs_no_survival')
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

    def test_binary_one_sex_with_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSBT', 'binary_one_sex_with_pcs_no_survival')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs_no_survival(metadata_dict)

    def test_binary_both_sexes_with_pcs_no_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_both_sexes_with_pcs_no_survival')
        pheno = add_both_sexes_to_phenotype(pheno)
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs_no_survival(metadata_dict)

    def test_binary_no_sex_no_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_no_sex_no_pcs_with_survival')
        pheno = add_survival_data_to_phenotype(pheno, 'ABT')
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_with_survival(metadata_dict)

    def test_binary_one_sex_no_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSBT', 'binary_one_sex_no_pcs_with_survival')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        pheno = add_survival_data_to_phenotype(pheno, 'ASSBT')
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_with_survival(metadata_dict)

    def test_binary_both_sexes_no_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_both_sexes_no_pcs_with_survival')
        pheno = add_both_sexes_to_phenotype(pheno)
        pheno = add_survival_data_to_phenotype(pheno, 'ABT')
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs_with_survival(metadata_dict)

    def test_binary_no_sex_with_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT', 'binary_no_sex_with_pcs_with_survival')
        pcs = create_pcs()
        pheno = add_survival_data_to_phenotype(pheno, 'ABT')
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs_with_survival(metadata_dict)

    def test_binary_one_sex_with_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSBT', 'binary_one_sex_with_pcs_with_survival')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        pheno = add_survival_data_to_phenotype(pheno, 'ASSBT')
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs_with_survival(metadata_dict)

    def test_binary_both_sexes_with_pcs_with_survival(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ABT',
                                                                       'binary_both_sexes_with_pcs_with_survival')
        pheno = add_both_sexes_to_phenotype(pheno)
        pcs = create_pcs()
        pheno = add_survival_data_to_phenotype(pheno, 'ABT')
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['auc'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs_with_survival(metadata_dict)


class TestComparePrsQuantitative(unittest.TestCase):
    """
    Simple tests that check the evaluate_prs method finishes successfully for all supported input combinations
    Also checks for the minimum expected output (that evaluation_metrics.csv produces some discrimination metrics,
    and that at least one plot is produced)
    """

    @classmethod
    def setUpClass(cls):
        cls._temp_output_dir = TemporaryDirectory()

    def _expected_plots_no_pcs(self, metadata_dict):
        for dir in ['prs1_All', 'prs2_All']:
            self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots', dir,
                                                        'box_plots_deciles.png')))
            self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots', dir,
                                                        'box_plots_risk.png')))

    def _expected_plots_with_pcs(self, metadata_dict):
        for label in ['prs1', 'prs2']:
            for anc in ['AFR', 'EAS', 'EUR', 'SAS']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots',
                                                            f'{label}_{anc}', 'box_plots_deciles.png')))
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'plots',
                                                            f'{label}_{anc}', 'box_plots_risk.png')))

                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_box_plot.png')))
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_prs_hist.png')))
            for pc in ['pc1', 'pc2', 'pc3', 'pc4']:
                self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'],
                                                            'plots', f'{label}_cross_ancestry',
                                                            f'{label}_{pc}_by_ancestry.png')))

    def _build_standard_inputs(self, trait_code: str, test_name: str):
        metadata = create_metadata_dict(trait_code)
        output_dir = os.path.join(self._temp_output_dir.name, test_name)
        os.mkdir(output_dir)
        os.mkdir(os.path.join(output_dir, 'plots'))
        metadata['output_dir'] = output_dir
        return create_prs1(), create_prs2(), create_base_quantitative_phenotype(trait_code), metadata

    def test_quantitative_no_sex_no_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('AQT', 'quant_no_sex_no_pcs')
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs(metadata_dict)

    def test_quantitative_one_sex_no_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSQT', 'quant_one_sex_no_pcs')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs(metadata_dict)

    def test_quantitative_both_sexes_no_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('AQT', 'quant_both_sexes_no_pcs')
        pheno = add_both_sexes_to_phenotype(pheno)
        compare_prs(prs1, prs2, pheno, metadata_dict, pandas.DataFrame())

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_no_pcs(metadata_dict)

    def test_quantitative_no_sex_with_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('AQT', 'quant_no_sex_with_pcs')
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs(metadata_dict)

    def test_quantitative_one_sex_with_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('ASSQT', 'quant_one_sex_with_pcs')
        pheno = add_single_sex_to_phenotype(pheno, 0)
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs(metadata_dict)

    def test_quantitative_both_sexes_with_pcs(self):
        prs1, prs2, pheno, metadata_dict = self._build_standard_inputs('AQT', 'quant_both_sexes_with_pcs')
        pheno = add_both_sexes_to_phenotype(pheno)
        pcs = create_pcs()
        compare_prs(prs1, prs2, pheno, metadata_dict, pcs)

        output_metrics = pandas.read_csv(os.path.join(metadata_dict['output_dir'], 'evaluation_metrics.csv'))
        self.assertFalse(output_metrics['rsq'].isnull().values.all())
        self.assertTrue(os.path.exists(os.path.join(metadata_dict['output_dir'], 'prs_evaluation_report.pdf')))
        self.assertNotEqual(len(os.listdir(os.path.join(metadata_dict['output_dir'], 'plots'))), 0)

        self._expected_plots_with_pcs(metadata_dict)
