"""
Scripts for plotting evaluation metrics using colouringbook
"""

import matplotlib.pyplot as pyplot
import numpy
import os
import pandas
import plotnine as pn
import pkg_resources
import seaborn

from lifelines import KaplanMeierFitter
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from typing import Dict
import yaml

from .constants import CUM_INC_X_LABEL_MAPPINGS, SUPPORTED_PC_HEADERS, PLOT_COLOURS, SEX_MAPPING

__all__ = ['theme_gplc', 'scale_colour_gplc', 'scale_fill_gplc']

yaml_stream = pkg_resources.resource_stream(__name__, './img/palette.yaml')
palette_gplc_hex = yaml.safe_load(yaml_stream)['light']


def plot_roc_curve(df_dict: Dict[str, pandas.DataFrame], output_path: str, phenotype_column: str, prs_column: str,
                   ancestry: str = 'All'):
    """
    df_dict: Dict of dataframes keyed by model name, where each df is indexed on IIDs and contains PRS,
    phenotype and ancestry data
    """
    lr = LogisticRegression(solver='liblinear')

    fpr_tpr_df = pandas.DataFrame()
    text_annot = pandas.DataFrame()

    top = 0.28
    for i, (model, df) in enumerate(df_dict.items()):
        top = top - 0.05
        temp = pandas.DataFrame()
        temp_text = pandas.DataFrame()
        input_df = df[[prs_column]]
        funct = lr.fit(input_df, df[phenotype_column])
        pred = pandas.DataFrame(funct.predict_proba(input_df), index=input_df.index)
        fpr, tpr, thresholds = roc_curve(df[phenotype_column], pred.iloc[:, 1])
        temp['fpr'] = fpr
        temp['tpr'] = tpr
        temp['model'] = model
        temp['ancestry'] = f'Ancestry = {ancestry}'
        auc = str(round(roc_auc_score(df[phenotype_column], pred.iloc[:, 1]), 3))
        temp_text.loc[i, 'label'] = 'AUC = {auc}'.format(auc=auc)
        temp_text.loc[i, 'model'] = model
        temp_text.loc[i, 'x'] = 0.5
        temp_text.loc[i, 'y'] = top
        fpr_tpr_df = fpr_tpr_df.append(temp)
        text_annot = text_annot.append(temp_text)

    p = (pn.ggplot() +
         pn.geom_line(fpr_tpr_df, pn.aes('fpr', 'tpr', color='model'), size=1, alpha=0.7) +
         pn.geom_abline(colour='grey', linetype='dashed') +
         theme_gplc() + scale_colour_gplc() +
         pn.facet_wrap('~ancestry', ncol=2) +
         pn.theme(legend_key=pn.element_blank(),
                  legend_position="none",
                  strip_background=pn.element_rect(fill="white"),
                  strip_text=pn.element_text(colour='black')) +
         pn.geom_text(pn.aes(x='x', y='y', color='model', label='label'), data=text_annot, ha='left', size=8) +
         pn.labs(title=prs_column,
                 x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)"))

    p.save(output_path, height=5, width=5, units='in', dpi=300, verbose=False)


def theme_gplc(base_size=12, backcol="white", xangle=0, legendposition="bottom", legendtitle="no", font="Poppins"):
    textcol = palette_gplc_hex['oxfordblue']
    linecol = palette_gplc_hex['oxfordblue']
    if backcol != "white":
        backcol = palette_gplc_hex['grey']
    if legendposition == "bottom":
        legend_position = (.5, .05)
        subplots_adjust = {'bottom': 0.17}
        legend_direction = "horizontal"
    else:
        legend_position = None
        subplots_adjust = None
        legend_direction = "vertical"
    if xangle == 0:
        xhjust = 0.5
    else:
        xhjust = 1
    newtheme = (pn.theme_bw(base_size=base_size, base_family=font) +
                pn.theme(
                    # set standards for all text
                    plot_title=pn.element_text(family=font, size=base_size+4, face="plain", color=textcol),
                    legend_title=pn.element_blank(),
                    legend_text=pn.element_text(family=font, size=base_size, colour=textcol),
                    axis_title=pn.element_text(family=font, size=base_size, face="plain", color=textcol, hjust=0.5),
                    axis_text=pn.element_text(family=font, size=base_size, color=textcol),
                    axis_text_x=pn.element_text(angle=xangle, hjust=xhjust),
                    # legend, default to bottom
                    legend_position=legend_position,
                    subplots_adjust=subplots_adjust,
                    legend_direction=legend_direction,
                    legend_background=pn.element_blank(),
                    legend_key=pn.element_blank(),
                    # ticks and axis line for x-axis, neither for y-axis
                    axis_ticks_major_y=pn.element_blank(),
                    axis_ticks_minor_y=pn.element_blank(),
                    axis_ticks_major_x=pn.element_line(color=linecol),
                    axis_line_y=pn.element_blank(),
                    axis_line_x=pn.element_line(colour=textcol, size=1.0),
                    # no minor gridlines, only horizontal major gridlines
                    panel_grid_minor=pn.element_blank(),
                    panel_grid_major_y=pn.element_line(color="#bebebe", size=0.5),
                    panel_grid_major_x=pn.element_blank(),
                    # set panel and/or plot background to be white
                    # (change backcol at the top of the function to change this)
                    panel_background=pn.element_rect(fill=backcol),  # just the area which contains the plotted data
                    plot_background=pn.element_rect(fill=backcol),  # the entire area
                    panel_border=pn.element_blank(),
                    # Settings for facet_wrap etc
                    strip_background=pn.element_blank(),
                    strip_text=pn.element_text(colour=textcol, size=base_size,  hjust=0.5))
                )
    if legendtitle == "yes":
        newtheme = newtheme + pn.theme(legend_title=pn.element_text(
            family=font, size=base_size, hjust=0.5, colour=textcol))
    return newtheme


class scale_colour_gplc(pn.scales.scale.scale_discrete):
    _aesthetics = ['color']
    na_value = '#7F7F7F'

    def __init__(self, **kwargs):
        self.palette = lambda n: list(palette_gplc_hex.values())[0:n]
        pn.scales.scale.scale_discrete.__init__(self, **kwargs)


class scale_fill_gplc(scale_colour_gplc):
    _aesthetics = ['fill']


def _filter_by_sex(df: pandas.DataFrame, sex: str):
    if sex != 'both':
        df = df.loc[df['sex'] == SEX_MAPPING[sex]]
    return df


def generate_binary_plots(df: pandas.DataFrame, trait_code: str, data_column_name: str,
                          sex: str, ancestry: str, phenotype_dict: dict, output_dir: str,):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plots_dict = dict()
    if sex == phenotype_dict['gender']:
        plots_dict['roc_curve'] = prepare_and_plot_roc(df, trait_code, data_column_name,
                                                       os.path.join(output_dir, 'roc_plot.png'), sex, ancestry)
    plots_dict['cum_inc'] = prepare_and_plot_cumulative_incidence(df, trait_code, data_column_name,
                                                                  os.path.join(output_dir, 'cum_inc_plot.png'), sex,
                                                                  ancestry)

    return plots_dict


def prepare_and_plot_roc(df: pandas.DataFrame, trait_code: str, data_column_name: str,
                         output_path: str, sex: str = 'both', ancestry: str = 'All'):
    """Prepares data for plotting ROC plot, creates the plot and returns the plots evaluation packet"""
    plot_roc_curve({'test_model': df}, output_path, trait_code, data_column_name, ancestry)
    description = f'ROC plot for phenotype {trait_code} for {data_column_name} (sex: {sex}. ancestry: {ancestry}.)'
    return {'path': output_path, 'description': description}


def prepare_and_plot_cumulative_incidence(df: pandas.DataFrame, trait_code: str,
                                          data_column_name: str, output_path: str, sex: str = 'both',
                                          ancestry: str = 'All'):
    """Prepares data for plotting cumulative incidence curves, creates the plot and returns the plots
    evaluation packet"""

    df_surv = _prepare_risk_groups(df, data_column_name)

    temporal_column = 'age_event'
    kmfdf = _prepare_survival_model(df_surv, temporal_column, trait_code)
    _plot_cumulative_incidence(kmfdf, trait_code, output_path, CUM_INC_X_LABEL_MAPPINGS[temporal_column],
                               data_column_name, 73)
    description = f'Cumulative incidence plot for phenotype {trait_code} for {data_column_name} ' \
                  f'(sex: {sex}. ancestry: {ancestry}.)'
    return {'path': output_path, 'description': description}


def _prepare_risk_groups(df: pandas.DataFrame, data_column_name: str):
    df['risk_group'] = numpy.nan
    prs_quantiles = list(df[data_column_name].quantile([0.03, 0.4, 0.6, 0.97, 0.99]))

    df.loc[(df[data_column_name] <= prs_quantiles[0]), 'risk_group'] = 'Top 3% protected'
    df.loc[((df[data_column_name] >= prs_quantiles[1]) &
            (df[data_column_name] <= prs_quantiles[2])), 'risk_group'] = 'Average risk'
    df.loc[((df[data_column_name] >= prs_quantiles[3]) &
            (df[data_column_name] <= prs_quantiles[4])), 'risk_group'] = 'Top 3% risk'
    df.loc[(df[data_column_name] > prs_quantiles[4]), 'risk_group'] = 'Top 3% risk'
    df = df[df.risk_group.notnull()]
    return df


def _prepare_survival_model(df_surv: pandas.DataFrame, temporal_column: str, event_column: str,
                            split_by: str = 'risk_group'):

    pandas.options.mode.chained_assignment = None
    kmfdf = pandas.DataFrame(columns=['time', 'risk', 'risk_lci', 'risk_uci', 'strata', 'ancestry'])
    for name, grouped_df in df_surv.groupby(split_by):
        kmf = KaplanMeierFitter()
        kmf.fit(grouped_df[temporal_column],
                grouped_df[event_column], label=name)
        tmp = kmf.cumulative_density_.join(kmf.confidence_interval_cumulative_density_)
        tmp = tmp.reset_index()
        tmp.columns = ['time', 'risk', 'risk_lci', 'risk_uci']
        tmp['strata'] = name
        kmfdf = pandas.concat([kmfdf, tmp])

    return kmfdf


def _plot_cumulative_incidence(kmfdf: pandas.DataFrame, trait_code: str, output_path: str, x_label: str,
                               prs_column_name: str, t_max: int = 73):
    def _format_percent(input_vals: list):
        return [f'{x * 100:.2g}%' for x in input_vals]

    t_min = kmfdf.query('risk > 0.005').groupby('strata').time.min().values
    if(len(t_min) != 3):
        t_min = numpy.append(t_min, kmfdf.query('time!=0').time.min())
    t_min = round(int(min(t_min) / 5)) * 5

    p = (pn.ggplot(data=kmfdf.query('time>=@t_min & time<=@t_max'),
                   mapping=pn.aes(x='time', y='risk', ymax='risk_uci', ymin='risk_lci', colour='strata')) +
         pn.geom_ribbon(alpha=0.3, colour=None, mapping=pn.aes(fill='strata')) +
         pn.geom_step(size=1) +
         theme_gplc() +
         pn.scale_y_continuous(labels=_format_percent) +
         pn.coord_cartesian(ylim=[0, kmfdf.query('time<=@t_max').risk.max() * 1.01]) +
         pn.scale_color_manual(values=['#110035', '#0d5eff', '#fc4183']) +
         pn.scale_fill_manual(values=['#110035', '#0d5eff', '#fc4183']) +
         pn.xlab(x_label) +
         pn.ylab(f"Percentage diagnosed with\n{trait_code}") +
         pn.labs(title=prs_column_name))

    fwidth, fheight = (6.5, 6)
    p.save(filename=output_path,
           height=fheight, width=fwidth, units='in', dpi=200, verbose=False)


def generate_quantitative_plots(df: pandas.DataFrame, trait_code: str,
                                data_column_name: str, output_dir: str, sex: str = 'both', ancestry: str = 'All'):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plots_dict = dict()
    plots_dict['box_plots_risk_groups'] = prepare_and_plot_box_risk_groups(df, trait_code,
                                                                           data_column_name,
                                                                           os.path.join(output_dir,
                                                                                        'box_plots_risk.png'),
                                                                           sex, ancestry)

    plots_dict['box_plots_deciles'] = prepare_and_plot_box_plot_deciles(df, trait_code,
                                                                        data_column_name,
                                                                        os.path.join(output_dir,
                                                                                     'box_plots_deciles.png'),
                                                                        sex, ancestry)
    return plots_dict


def prepare_and_plot_box_risk_groups(df: pandas.DataFrame, trait_code: str,
                                     data_column_name: str, output_path: str, sex: str = 'both',
                                     ancestry: str = 'All'):
    """Prepares data for plotting box plots for risk groups, creates the plot and returns the plots evaluation packet"""
    _plot_box_plot_risk_groups(df, output_path, data_column_name, trait_code)
    description = f'Box plots for PRS risk groups: phenotype {trait_code} for {data_column_name}' \
                  f' (sex: {sex}. ancestry: {ancestry}.)'
    return {'path': output_path, 'description': description}


def _plot_box_plot_risk_groups(df: pandas.DataFrame, output_path: str, data_column_name: str, trait_code: str,
                               ancestry: str = 'All'):
    """Method that generates a box plot that illustrates the phenotype distributions (on raw scale) for the
    population in the top and bottom 3% PRS distribution, and the average risk population (40th - 60th percentile)
    """

    prs_percentiles = list(df[data_column_name].quantile([0.03, 0.4, 0.6, 0.97]))
    df.loc[(df[data_column_name] <= prs_percentiles[0]), 'risk_group'] = 'Top 3% protected'
    df.loc[((df[data_column_name] > prs_percentiles[1]) & (df[data_column_name] <= prs_percentiles[2])),
           'risk_group'] = 'Average risk'
    df.loc[(df[data_column_name] > prs_percentiles[3]), 'risk_group'] = 'Top 3% risk'

    seaborn.set(font_scale=1.2)
    seaborn.set_style("ticks")

    pyplot.rcParams["axes.labelsize"] = 17

    ctpyplot = seaborn.catplot(x="risk_group", y=trait_code, data=df, kind="box",
                               order=['Top 3% protected', 'Average risk', 'Top 3% risk'],
                               palette=seaborn.color_palette(['#53b24a', '#3a7bb8', '#e0001e']))
    ctpyplot.set_axis_labels("", trait_code + ' values')
    ctpyplot.set_titles(ancestry, size=18)
    [pyplot.setp(ax.get_xticklabels(), rotation=-45, ha='left') for ax in ctpyplot.axes.flat]
    pyplot.title(data_column_name)
    pyplot.savefig(output_path, bbox_inches='tight', dpi=500)
    pyplot.close()
    pyplot.style.use('default')


def prepare_and_plot_box_plot_deciles(df: pandas.DataFrame, trait_code: str,
                                      data_column_name: str, output_path: str, sex: str = 'both',
                                      ancestry: str = 'All'):
    """Prepares data for plotting box plots for risk groups, creates the plot and returns the plots evaluation packet"""
    _plot_box_plot_deciles(df, output_path, data_column_name, trait_code)
    description = f'Box plots for PRS deciles: phenotype {trait_code} for {data_column_name} ' \
                  f'(sex: {sex}. ancestry: {ancestry}.)'
    return {'path': output_path, 'description': description}


def _plot_box_plot_deciles(df: pandas.DataFrame, output_path: str, data_column_name: str, trait_code: str,
                           ancestry: str = 'All'):
    decs = list(range(1, 11))
    try:
        df['decile'] = pandas.qcut(df[data_column_name], 10, labels=decs, duplicates='drop')
    except ValueError:
        df['decile'] = pandas.qcut(df[data_column_name].rank(method='first'), 10,
                                   labels=decs, duplicates='drop')

    seaborn.set_style("ticks")
    pyplot.rcParams["axes.labelsize"] = 16
    ctplt = seaborn.catplot(x="decile", y=trait_code, data=df, kind="box", palette='YlOrRd')
    ctplt.set_axis_labels("PRS deciles", trait_code + ' values')
    ctplt.set_titles(ancestry)
    pyplot.title(data_column_name)
    pyplot.savefig(output_path, bbox_inches='tight', dpi=500)
    pyplot.close()
    pyplot.style.use('default')


def generate_cross_ancestry_plots(df: pandas.DataFrame, prs_column_name: str, output_dir: str, sex: str = 'both'):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plots_dict = dict()
    plots_dict['prs_hist'] = prepare_and_plot_prs_histogram_per_ancestry(
        df, prs_column_name, os.path.join(output_dir, f'{prs_column_name}_prs_hist.png'), sex)

    plots_dict['prs_box_plot'] = prepare_and_plot_prs_box_plot_per_ancestry(
        df, prs_column_name, os.path.join(output_dir, f'{prs_column_name}_prs_box_plot.png'), sex)

    if all(p in df.columns for p in SUPPORTED_PC_HEADERS):
        for p in SUPPORTED_PC_HEADERS:
            plots_dict[f'{p}_binned_prs'] = prepare_and_plot_genetic_pcs_binned_by_prs(
                df, prs_column_name, p, os.path.join(output_dir, f'{prs_column_name}_{p}_by_ancestry.png'), sex)

    return plots_dict


def plot_prs_binned(df: pandas.DataFrame, prs_column_name: str, bin_by, zoom=None):
    """
    Can be used to either plot the predictors (PCs) against the uncentered PRS
    Or plot the predictors against the residuals (the centered PRS).
    Or plot the overall predicted prs against the residuals.
    """
    pcbins = _bin_prs(df, bin_by, prs=prs_column_name)

    alpha = 10 / len(df)**0.5
    p = pn.ggplot(df, pn.aes(bin_by, prs_column_name, colour="ancestry")) + \
        pn.geom_hline(yintercept=0, linetype='dashed') + \
        pn.geom_point(alpha=alpha) + \
        pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1})) + \
        pn.scale_colour_brewer(type="qual", palette="Set1") + \
        pn.geom_smooth(method='lm', colour='blue') + \
        pn.geom_pointrange(data=pcbins, mapping=pn.aes(f'bin_mean', 'prs_mean', ymin='prs_lci', ymax='prs_uci'),
                           colour='black') + \
        pn.labs(title=prs_column_name)
    if zoom is not None:
        p = p + pn.coord_cartesian(ylim=(-zoom, zoom))
    return p


def _bin_prs(df, bin_by: str, prs: str = 'prs_centered'):
    # Calculate ancestry in a pseudo-continuous (binned) method
    df[f'{bin_by}_bins'] = pandas.cut(df[f'{bin_by}'], 50)

    pcbins = df.groupby(f'{bin_by}_bins').agg(
                bin_mean=(f"{bin_by}", "mean"),
                N=(f'{bin_by}', lambda x: len(x)),
                prs_mean=(prs, "mean"),
                prs_sem=(prs, lambda x: numpy.std(x, ddof=1) / numpy.sqrt(len(x))))

    pcbins['prs_lci'] = pcbins['prs_mean'] - pcbins['prs_sem'] * 1.96
    pcbins['prs_uci'] = pcbins['prs_mean'] + pcbins['prs_sem'] * 1.96
    pcbins = pcbins.query("N>=10")
    return pcbins


def prepare_and_plot_genetic_pcs_binned_by_prs(df: pandas.DataFrame, prs_column_name: str, pc_column_name: str,
                                               output_path: str, sex: str = 'both'):
    _plot_genetic_pcs_binned_by_prs(df, prs_column_name, pc_column_name, output_path)
    description = f'Plot displaying the binned PRS against genetically inferred principal component {pc_column_name},' \
                  f' with data coloured by genetically inferred ancestry (sex: {sex})'
    return {'path': output_path, 'description': description}


def prepare_and_plot_prs_histogram_per_ancestry(df: pandas.DataFrame, prs_column_name: str, output_path: str,
                                                sex: str = 'both'):
    _plot_prs_histogram_per_ancestry(df, prs_column_name, output_path)
    description = f'PRS density distributions per ancestry (sex: {sex})'
    return {'path': output_path, 'description': description}


def prepare_and_plot_prs_box_plot_per_ancestry(df: pandas.DataFrame, prs_column_name: str, output_path: str,
                                               sex: str = 'both'):
    _plot_prs_box_plot_per_ancestry(df, prs_column_name, output_path)
    description = f'PRS box plot distributions for each ancestry (sex: {sex})'
    return {'path': output_path, 'description': description}


def _plot_genetic_pcs_binned_by_prs(df: pandas.DataFrame, prs_column_name: str, pc_column_name: str, output_path: str):
    p = plot_prs_binned(df, prs_column_name, pc_column_name)
    fwidth, fheight = (9.6, 7.2)
    p.save(filename=output_path, height=fheight, width=fwidth, units='in', dpi=200, verbose=False)


def _plot_prs_histogram_per_ancestry(df: pandas.DataFrame, prs_column_name: str, output_path: str):

    pyplot.rcParams['axes.labelsize'] = 16
    pyplot.clf()
    pyplot.figure(figsize=(10, 7))
    ancestries = list(df['ancestry'].unique())

    for anc, col in zip(ancestries, PLOT_COLOURS[0:len(ancestries)]):
        df_a = df[df['ancestry'] == anc]

        seaborn.kdeplot(df_a[prs_column_name],
                        linewidth=3,
                        color=col,
                        label=anc)
        pyplot.axvline(df_a[prs_column_name].mean(), linestyle="--", color='grey')

    pyplot.xticks(fontsize=16)
    pyplot.yticks(fontsize=16)
    pyplot.legend(prop={'size': 16}, title='Ancestry')
    pyplot.xlabel('PRS')
    pyplot.ylabel('Density')
    pyplot.title(prs_column_name)
    pyplot.savefig(output_path, bbox_inches='tight', dpi=500)
    pyplot.close()


def _plot_prs_box_plot_per_ancestry(df: pandas.DataFrame, prs_column_name: str, output_path: str):
    pyplot.rcParams['axes.labelsize'] = 16
    pyplot.clf()
    pyplot.figure(figsize=(10, 7))
    ancestries = list(df['ancestry'].unique())

    seaborn.set_style('ticks')
    ctplt = seaborn.catplot(x='ancestry', y=prs_column_name, data=df, kind='box',
                            palette={k: v for k, v in zip(ancestries, PLOT_COLOURS[0:len(ancestries)])})
    ctplt.set_axis_labels('Ancestry', 'PRS distribution')
    pyplot.title(prs_column_name)
    pyplot.savefig(output_path, bbox_inches='tight', dpi=500)
    pyplot.close()
    pyplot.style.use('default')
