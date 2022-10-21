"""
Generates a report comparing input PRS to the Genomics plc standard
"""
from datetime import datetime
from fpdf import FPDF
import os
import pandas

from .constants import REPORT_METRIC_MAPPINGS, SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_BINARY, \
    SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_SAMPLES_ONLY, BINARY_METRICS_TO_INCLUDE_IN_REPORT, \
    QUANTITATIVE_METRICS_TO_INCLUDE_IN_REPORT, REPORT_ANCESTRY_ORDER, \
    SUPPORTED_PC_HEADERS, REPORT_SEX_ORDER, FONT_DIR, MIN_EVALUATION_SAMPLES_PER_GROUP
from .error import UkbPretReportError

IMAGE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'img')

os.environ['FPDF_FONTPATH'] = FONT_DIR


"""
Example results packet:

test = {'phenotype_code': 'AST',
        'phenotype_description': 'Asthma',
        'flavour': 'binary',
        'evaluation': {
            'input_data': {
                'male': {
                    'AFR': {'metrics': {'discrimination': {'value': 0.5, 'lci': 0.4, 'uci': 0.6},
                                               'odds_ratio_per_prs': {'value': ...},
                                               'odds_ratio_per_std': {'value': ...}},
                                   'plots': {'roc_plot': {'path': '/path/to/plot', 'description': ''},
                                             'other_plot': {...}, ...},
                                   'prs_sumstats': {'mean': 0.5, 'mean_lci': 0.35, 'mean_uci': 0.58, 'std': 0.02},
                                   'sample_data': {'n_samples': 1200, 'n_cases': 540, 'n_controls': 660,
                                                   'n_prevalent': 100, 'n_incident: 500} (NB - just n_samples for quant)
                               },
                    'All': {'metrics': {'discrimination': {'value': 0.5, 'lci': 0.4, 'uci': 0.6},
                                               'odds_ratio_per_prs': {'value': ...},
                                               'odds_ratio_per_std': {'value': ...}},
                                   'plots': {'roc_plot': '/path/to/plot', ...},
                                   'prs_sumstats': {'mean': 0.5, 'mean_lci': 0.35, 'mean_uci': 0.58, 'std': 0.02},
                                   'sample_data': {...}
                               },
                    },
                'female': {...},
                'both': {...}
                },
            'Gplc': {
                'male': {
                    'EAS': {'metrics': {'discrimination': {'value': 0.5, 'lci': 0.4, 'uci': 0.6},
                                        'odds_ratio_per_prs': {'value': ...},
                                        'odds_ratio_per_std': {'value': ...}},
                            'plots': {'roc_plot': {'path': '/path/to/plot', 'description': ''},
                            'prs_sumstats': {'mean': 0.5, 'mean_lci': 0.35, 'mean_uci': 0.58, 'std': 0.02},
                            'sample_data': {...}
                           }, ...
                    },
                'female': {... }, 'both': {...}
            }
        },
        'cross_ancestry_evaluation': {
            'input_data': {
                'male': {
                    'plots': { 'pc_prs_by_ancestry': { 'path': '/a/path', description: ''}, ... }
                },
                'female': {...}, 'both': {...}
            }, ...
        'metadata':{
            'input_data': {'n_removed_prs': {...}, 'n_removed_pheno': {...}},
            'Gplc': {'n_removed_prs': {...}, 'n_removed_pheno': {...}}
        }
    }
"""


def generate_report(results, output_path: str):

    pdf = UkbPretReport(results)
    pdf.add_page()
    pdf.report_title()
    pdf.preamble()
    pdf.next_page()
    pdf.n_samples_table()
    pdf.next_page()
    pdf.metrics_table()
    pdf.next_page()

    pdf.eval_plots()

    pdf.output(output_path, 'F')
    cleanup_font_dir()


def cleanup_font_dir():
    to_cleanup = [fn for fn in os.listdir(FONT_DIR) if '.pkl' in fn]
    for fn in to_cleanup:
        os.remove(os.path.join(FONT_DIR, fn))


class UkbPretReport(FPDF):

    text_size = {'s': 10, 'r': 12, 'title1': 42, 'title2': 28, 'title3': 20, 'title4': 16, 'vs': 8}

    def __init__(self, results):
        super().__init__(font_cache_dir=FONT_DIR)

        self.results = results
        self.data_sources = list(self.results['evaluation'].keys())
        self.all_sexes = list(set([sex for source in self.data_sources for sex in self.results['evaluation'][source]]))
        self.sex_major = 'unspecified' if self.all_sexes == ['unspecified'] else results['sex_major']
        ordered_sexes = list()
        for sex in REPORT_SEX_ORDER:
            if sex in self.all_sexes:
                ordered_sexes.append(sex)
        self.all_sexes = ordered_sexes
        ancestries = set()
        for data_source in self.data_sources:
            for sex in self.all_sexes:
                for ancestry in self.results['evaluation'][data_source][sex].keys():
                    ancestries.add(ancestry)
        all_ancestries = []
        for a in REPORT_ANCESTRY_ORDER:
            if a in ancestries:
                all_ancestries.append(a)
        self.all_ancestries = all_ancestries

        # Add metrics
        supported_metrics = [set(eval_data['metrics']) for sex_data in
                             self.results['evaluation'][self.data_sources[0]].values()
                             for eval_data in sex_data.values()]
        if not all(m == supported_metrics[0] for m in supported_metrics):
            raise UkbPretReportError("Entries should contain the same metric fields")
        supported_metrics = list(supported_metrics[0])

        self.metrics = []
        if self.results['flavour'] == 'binary':
            for metric in BINARY_METRICS_TO_INCLUDE_IN_REPORT:
                if metric in supported_metrics:
                    self.metrics.append(metric)
        else:
            for metric in QUANTITATIVE_METRICS_TO_INCLUDE_IN_REPORT:
                if metric in supported_metrics:
                    self.metrics.append(metric)

        self.table_buffer = 1.2

        # Add Poppins to supported fonts
        font_filenames = [fn for fn in os.listdir(FONT_DIR) if '.pkl' not in fn]
        for font_filename in font_filenames:
            filename_without_suffix = os.path.splitext(font_filename)[0]
            font_name = 'Poppins' if filename_without_suffix == 'Poppins-Regular' else filename_without_suffix
            self.add_font(font_name, fname=os.path.join(FONT_DIR, font_filename), uni=True)

        self.current_page = 0

    def header(self):
        self.current_page += 1
        self.image(os.path.join(IMAGE_DIR, 'genomics_logo_2021.png'), x=0, y=0, w=self.w)
        self._spacing('r')

    def footer(self):
        self.set_y(-20)
        self.set_font('poppins-italic', '', self.text_size['vs'])
        self.cell(0, 0, "This report was generated using Genomics plc's UK Biobank PRS Release Evaluation Tool.",
                  align='C')
        self._spacing('vs')
        license_text = "This software is available to use under the licence "
        self.cell(0, 0, license_text, align='C')
        self.set_x((self.w / 2) + (self.get_string_width(license_text) / 2) * 0.98)
        self.set_text_color(0, 0, 255)
        self.set_font(style="U")
        self.cell(0, 0, "here", link="https://github.com/Genomicsplc/ukb-pret/blob/main/LICENCE")
        self.set_font('poppins')
        self.set_text_color(0)
        self.cell(0, 0, f'{self.current_page}', align='R')

    def next_page(self):
        self.add_page()

    def preamble(self):
        self.set_font('poppins', size=self.text_size['r'])
        self._simple_text("This report evaluates the Phenotype:")
        self.set_font('poppins-bold', size=self.text_size['r'])
        self._simple_text(f"    {self.results['phenotype_description']} ({self.results['phenotype_code']})")
        self.set_font('poppins', size=self.text_size['r'])
        self._simple_text("For the Polygenic Risk Scores (PRSs) described by the labels:")
        self.set_font('poppins-bold', size=self.text_size['r'])
        for dataset in self.results['evaluation'].keys():
            self._simple_text(f"    {dataset}")
        self.write_sample_summary()

    def _simple_text(self, text):
        self.cell(0, 0, text)
        self._spacing('s')

    def write_sample_summary(self):
        self.set_font('poppins-bold', size=self.text_size['title3'])
        self._spacing('n')
        self.cell(0, 0, 'Summary of the PRS sample inputs:')
        self._spacing('n')
        self.set_font('poppins', size=self.text_size['r'])

        for name, data in self.results['metadata'].items():
            self.set_font('poppins-bold', size=self.text_size['r'])
            self._simple_text(name + ':')
            self.set_font('poppins', size=self.text_size['r'])
            self._simple_text('  ' + data['n_input_samples']['description'])
            self._simple_text('  ' + data['n_removed_prs_tot']['description'])
            self._simple_text('  ' + data['n_removed_pheno']['description'])

        if self.results['metadata'][self.data_sources[0]]['n_included_samples']['value'] != \
                self.results['metadata'][self.data_sources[1]]['n_included_samples']['value']:
            raise AssertionError(f'The number of evaluated samples for {self.data_sources[0]} and '
                                 f'{self.data_sources[1]} should be the same.')

        self._spacing('s')
        self.set_font('poppins-bold', size=self.text_size['title4'])
        self.cell(0, 0, f"Total number of evaluated samples:"
                        f" {self.results['metadata'][self.data_sources[0]]['n_included_samples']['value']}")

    def report_title(self):
        self.ukb_pret_title('PRS evaluation report', 'title1', align='L')
        self.ukb_pret_title('Facilitated by the UK Biobank PRS release', 'title3', spacing_size='s', align='L')
        self.set_font('poppins-italic', size=self.text_size['r'])
        self._simple_text('Generated by the UK Biobank PRS Release Evaluation Tool')
        self._spacing('s')

    def ukb_pret_title(self, content: str, size: str = 'title1', spacing_size: str = 'n', align: str = 'L'):
        """Specify title of the report. The size can be set to ['l', 'r', 's']. """
        self.set_font('poppins-bold', '', self.text_size[size])
        self.multi_cell(0, 7, content, align=align)
        self._spacing(spacing_size)

    def table(self, data: pandas.DataFrame, title: str, table_width: float = None):
        """Populates a table of data with headings"""
        if not all(len(x) == len(data.iloc[0]) for _, x in data.iterrows()):
            raise AssertionError('Please ensure all rows in input DataFrame contain the same number of columns')

        if not table_width:
            table_width = (self.w - 2*self.l_margin) / len(data.columns)

        list_data = [[x for x in data.columns]]
        for i, row in data.iterrows():
            list_data.append([item for item in row])

        formatting = []
        for i, row in enumerate(list_data):
            if i == 0:
                self.set_font('poppins-bold', '', size=self.text_size['s'])
            else:
                self.set_font('poppins', '', size=self.text_size['s'])
            cell_thickness = float(self.font_size) + self.table_buffer

            row_height_lines = 1
            lines_in_row = []

            for item in row:
                output = self.multi_cell(table_width, cell_thickness, str(item), align='C', border=1, ln=3,
                                         split_only=True)
                lines_in_row.append(len(output))
                if len(output) > row_height_lines:
                    row_height_lines = len(output)

            formatting.append({'row_height_lines': row_height_lines, 'lines_in_row': lines_in_row, 'row': row,
                               'cell_thickness': cell_thickness})

        # Determine whether table will extend over next page
        vertical_whitespace = 0
        for formatting_dict in formatting:
            vertical_whitespace += (formatting_dict['row_height_lines'] * formatting_dict['cell_thickness'])
        if self.get_y() + vertical_whitespace > self.h - 20:
            self.next_page()

        for i, formatting_dict in enumerate(formatting):
            if i == 0:
                self.set_font('poppins-bold', '', size=self.text_size['s'])
            else:
                self.set_font('poppins', '', size=self.text_size['s'])

            for tlines, datum in zip(formatting_dict['lines_in_row'], formatting_dict['row']):
                text = datum.rstrip('\n') + (1 + formatting_dict['row_height_lines'] - tlines) * '\n'
                self.multi_cell(table_width, formatting_dict['cell_thickness'], text, align='C', ln=3, border=1)
            self.ln(formatting_dict['row_height_lines'] * formatting_dict['cell_thickness'])

    def plot(self, path: str):
        self.set_font('poppins', '')
        self.image(path, w=(self.w - 2*self.l_margin) * 0.7)
        self._spacing()

    def paragraph(self, text: str, size: str = 'r'):
        self.set_font('poppins', '', self.text_size[size])
        self.multi_cell(0, 6, text, ln=1)

    def metrics_table(self):
        self.ukb_pret_title('Evaluation Metrics', size='title2', spacing_size='s')

        description = 'Table containing evaluation metrics that describe the input Polygenic Risk Scores, given ' \
                      'alongside their lower and upper 95% confidence intervals (LCI and UCI respectively). \n' \
                      'The metrics include:'
        self.paragraph(description, size='s')
        if self.results['flavour'] == 'binary':
            self.paragraph('- Conditional Odds Ratio per 1 standard deviation in PRS (Conditional OR), which '
                           'represents the relative difference in disease risk per one standard deviation higher PRS.'
                           ' This was calculated using a logistic regression model that includes sex and age at first '
                           'assessment as covariates (or just age at first assessment for single sex analyses). '
                           'For ORs in different percentiles, please refer to the output CSV file.',
                           size='s')
            self.paragraph('- Conditional Hazard Ratio per 1 standard deviation in PRS (Conditional HR), '
                           'which is an alternative estimate of the difference in disease risk per one standard '
                           'deviation higher PRS. This was calculated using a Cox proportional hazards survival model '
                           'that includes sex and age at first assessment as covariates (or just age at '
                           'first assessment for single sex analyses).',
                           size='s')
            self.paragraph('- Area Under Curve (AUC) calculated from a Receiver Operator Characteristic (ROC) curve, '
                           'which quantifies the model\'s ability to distinguish between a true and false positive.',
                           size='s')
        else:
            self.paragraph('- R-squared (Rsq), which indicates the proportion of the variation in the quantitative '
                           'trait which can be attributed to the PRS.', size='s')
            self.paragraph('- Beta per 1 standard deviation in PRS (beta per std), which is an estimate of the change '
                           'in the trait measure per one standard deviation higher PRS assuming a linear correlation.',
                           size='s')

        self._spacing('s')

        table_dict = {'Ancestry': list(), 'PRS data source': list()}
        for metric in self.metrics:
            table_dict[f'{REPORT_METRIC_MAPPINGS[metric]} \n(LCI, UCI)'] = list()

        for ancestry in self.all_ancestries:
            for source in self.data_sources:
                if ancestry in self.results['evaluation'][source][self.sex_major].keys():
                    table_dict['Ancestry'].append(ancestry)
                    table_dict['PRS data source'].append(source)
                    for metric in self.metrics:
                        if not self.results['evaluation'][source][self.sex_major][ancestry]['metrics'][metric]['value']:
                            table_dict[f'{REPORT_METRIC_MAPPINGS[metric]} \n(LCI, UCI)'].append(
                                'N/A'
                            )
                        else:
                            table_dict[f'{REPORT_METRIC_MAPPINGS[metric]} \n(LCI, UCI)'].append(
                                self._format_metric(source, self.sex_major, ancestry, metric))

            # TODO: Add checks to ensure df contains what is expected...
        self.table(pandas.DataFrame.from_dict(table_dict), '')

    def n_samples_table(self):
        self.ukb_pret_title('Sample Summary table', size='title2', spacing_size='r')
        if self.results['flavour'] == 'binary' and self.results['survival_data']:
            table_dict = self._format_table_dict(SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_BINARY,
                                                 'sample_data', None)
            self.ukb_pret_title("The number of samples, cases and controls split by PRS, sex and ancestry",
                                'title3')
            self.set_font('poppins', size=self.text_size['r'])
            self._simple_text(f'Please note that strata where the number of cases/controls < '
                              f'{MIN_EVALUATION_SAMPLES_PER_GROUP} will not be evaluated')
        else:
            table_dict = self._format_table_dict(SAMPLE_NUMBER_TABLE_TO_RESULTS_MAPPING_SAMPLES_ONLY,
                                                 'sample_data', None)
            self.ukb_pret_title("Number of samples split by PRS, sex and ancestry", 'title3')
            self.set_font('poppins', size=self.text_size['r'])
            self._simple_text(f'Please note that strata where the number of samples < '
                              f'{MIN_EVALUATION_SAMPLES_PER_GROUP} will not be evaluated')
        self.table(pandas.DataFrame.from_dict(table_dict), f'')

    def _format_table_dict(self, schema: dict, metric: str, precision: int = None):
        table_dict = {'Ancestry': list(), 'Sex': list()}
        for key in schema.keys():
            table_dict[key] = list()
        for ancestry in self.all_ancestries:
            for sex in self.all_sexes:
                dummy_data_source = list(self.results['evaluation'].keys())[0]
                if sex in self.results['evaluation'][dummy_data_source].keys():
                    if ancestry in self.results['evaluation'][dummy_data_source][sex].keys():
                        table_dict['Ancestry'].append(ancestry)
                        table_dict['Sex'].append(sex)
                        for table_key, result_key in schema.items():
                            table_dict[table_key].append(
                                self._format_sumstats(dummy_data_source, sex, ancestry, metric, result_key, precision))
        return table_dict

    def eval_plots(self):
        self.ukb_pret_title('Evaluation plots', size='title2', spacing_size='vs')
        self.paragraph(f"Evaluation plots comparing user-provided data {self.data_sources[0]} against "
                       f"UK Biobank PRS Release {self.data_sources[1]}")
        if len(self.results['evaluation'].keys()) != 2:
            raise AssertionError('Evaluation report can only be generated for comparing 2 PRS')

        if self.results['flavour'] == 'binary':
            self.binary_eval_plots()
        else:
            self.quantitative_eval_plots()

        if self.all_ancestries != ['All']:
            self.qc_eval_plots()

    def binary_eval_plots(self):
        roc_description = f"Receiver Operating Characteristic (ROC) curve for phenotype " \
                          f"{self.results['phenotype_code']} " \
                          "in genetically inferred ancestry {}"
        cum_inc_description = f"Cumulative Incidence Plot for phenotype {self.results['phenotype_code']} " \
                              "in genetically inferred ancestry {} with sex = {}"

        for ancestry in self.all_ancestries:
            plots_made1 = self.do_comparison_plots(ancestry, 'roc_curve', roc_description.format(ancestry))
            plots_made2 = self.do_comparison_plots(ancestry, 'cum_inc', cum_inc_description.format(
                ancestry, self.sex_major))
            if plots_made1 or plots_made2:
                self.next_page()
            if self.sex_major == 'both':
                plots_made1 = self.do_comparison_plots(ancestry, 'cum_inc',
                                                       cum_inc_description.format(ancestry, 'female'), 'female')
                plots_made2 = self.do_comparison_plots(ancestry, 'cum_inc',
                                                       cum_inc_description.format(ancestry, 'male'), 'male')
                if plots_made1 or plots_made2:
                    self.next_page()

    def quantitative_eval_plots(self):
        box_plot_risk_description = f"Box plots in different PRS groups for phenotype " \
                                    f"{self.results['phenotype_code']} " \
                                     "in genetically inferred ancestry {}"
        box_plot_deciles_description = f"Box plots in PRS deciles for phenotype {self.results['phenotype_code']} " \
                                       "in genetically inferred ancestry {}"

        for ancestry in self.all_ancestries:
            self.do_comparison_plots(ancestry, 'box_plots_risk_groups', box_plot_risk_description.format(ancestry))
            self.do_comparison_plots(ancestry, 'box_plots_deciles', box_plot_deciles_description.format(ancestry))
            self.next_page()

    def qc_eval_plots(self):
        self.ukb_pret_title('Quality Control Report', size='title2', spacing_size='s')
        self.paragraph(f"QC plots comparing user-provided data {self.data_sources[0]} against "
                       f"UK Biobank PRS Release {self.data_sources[1]}")
        self._spacing('s')
        prs_hist_description = f"PRS density distribution for phenotype {self.results['phenotype_code']} split by " \
                               f"genetically inferred ancestry"
        prs_box_plot_description = f"PRS box plot distributions for phenotype {self.results['phenotype_code']} " \
                                   f"split by genetically inferred ancestry"

        self.do_comparison_plots_cross_ancestry('prs_hist', prs_hist_description)
        self.do_comparison_plots_cross_ancestry('prs_box_plot', prs_box_plot_description)
        self.next_page()

        if self._contains_pc_plots():

            generic_title = f"PRS against genetic principal components for phenotype " \
                            f"{self.results['phenotype_code']}"
            subtitle = "Each coloured point represents an individual, coloured by their genetically inferred " \
                       "ancestry grouping. Each black dot is an equal-width binned mean shown with a 95% confidence " \
                       "interval. Binned means with 10 or fewer constituent data points are not shown. " \
                       "The blue line is a simple linear regression. The black dashed line shows y=0."

            self.ukb_pret_title(generic_title, size='title3', spacing_size='s')
            self.paragraph(subtitle, size='s')
            self.do_comparison_plots_cross_ancestry('pc1_binned_prs', "PRS against genetic principal component 1")
            self.do_comparison_plots_cross_ancestry('pc2_binned_prs', "PRS against genetic principal component 2")
            self.next_page()

            self.ukb_pret_title(generic_title, size='title3', spacing_size='s')
            self.paragraph(subtitle, size='s')
            self.do_comparison_plots_cross_ancestry('pc3_binned_prs', "PRS against genetic principal component 3")
            self.do_comparison_plots_cross_ancestry('pc4_binned_prs', "PRS against genetic principal component 4")

    def cross_ancestry_eval_plots(self, source):

        if not all(len(sex_data.keys()) == 1 for sex_data in self.results['evaluation'][source].values()):
            self.ukb_pret_title('Cross-Ancestry evaluation plots', size='title3', spacing_size='s')
        for sex_data in self.results['cross_ancestry_evaluation'][source].values():
            for plot in sex_data['plots'].values():
                self.ukb_pret_title(plot['description'], size='title3', spacing_size='n')
                self.plot(plot['path'])
                self.next_page()

    def do_comparison_plots(self, ancestry: str, plot_results_key: str, title: str, sex: str = None):
        source1, source2 = self.results['evaluation'].keys()
        if sex is None:
            sex = self.sex_major

        if plot_results_key in self.results['evaluation'][source1][sex][ancestry]['plots'].keys() and \
                plot_results_key in self.results['evaluation'][source2][sex][ancestry]['plots'].keys():
            self.side_by_side_plots(self.results['evaluation'][source1][sex][ancestry]['plots'][plot_results_key],
                                    self.results['evaluation'][source2][sex][ancestry]['plots'][plot_results_key],
                                    title)
            return True
        else:
            return False

    def do_comparison_plots_cross_ancestry(self, plot_results_key: str, title: str, sex: str = None):

        source1, source2 = self.results['evaluation'].keys()
        if sex is None:
            sex = self.sex_major

        self.side_by_side_plots(self.results['cross_ancestry_evaluation'][source1][sex]['plots'][plot_results_key],
                                self.results['cross_ancestry_evaluation'][source2][sex]['plots'][plot_results_key],
                                title)

    def side_by_side_plots(self, plot_dict_1: dict, plot_dict_2: dict, title: str = None):

        if title is None:
            title = plot_dict_1['description']

        self.ukb_pret_title(title, size='r', spacing_size='s')
        self.set_font('poppins', '')
        prev_y = self.y
        self.image(plot_dict_1['path'], w=self.epw * 0.44)  # third page height, 2/5 page width
        self.set_y(prev_y)
        self.image(plot_dict_2['path'], w=self.epw * 0.44,
                   x=self.epw / 2)  # full page height, half page width, right half of the page

    def _contains_pc_plots(self, sex: str = None) -> bool:
        source1, source2 = self.results['evaluation'].keys()
        if sex is None:
            sex = self.sex_major

        pc_eval_file_names = [f'{p}_binned_prs' for p in SUPPORTED_PC_HEADERS]
        if all(p in self.results['cross_ancestry_evaluation'][source1][sex]['plots'].keys()
               for p in pc_eval_file_names) \
            and all(p in self.results['cross_ancestry_evaluation'][source2][sex]['plots'].keys()
                    for p in pc_eval_file_names):

            return True
        else:
            return False

    def _format_metric(self, source: str, sex: str, ancestry: str, metric: str = 'discrimination'):
        return (f"{self.results['evaluation'][source][sex][ancestry]['metrics'][metric]['value']:.3} \n"
                f"({self.results['evaluation'][source][sex][ancestry]['metrics'][metric]['lci']:.3},"
                f" {self.results['evaluation'][source][sex][ancestry]['metrics'][metric]['uci']:.3})")

    def _format_sumstats(self, source: str, sex: str, ancestry: str, metric: str, result_key: str,
                         precision: int = None):
        if precision is not None:
            return f"{str(self.results['evaluation'][source][sex][ancestry][metric][result_key]):.{precision}}"
        else:
            return f"{str(self.results['evaluation'][source][sex][ancestry][metric][result_key])}"

    def _spacing(self, size='r'):
        """Library for supported line spacing options"""
        size_dict = {'l': 40, 'r': 20, 'n': 12, 'sn': 10, 's': 8, 'vs': 4}
        self.ln(size_dict[size])
