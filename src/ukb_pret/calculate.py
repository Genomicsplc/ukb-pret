"""
Scripts for calculating evaluation metrics
"""

import numpy
import pandas
from scipy import stats
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score, r2_score
import statsmodels.api as sm
import statsmodels.formula.api as smf

from lifelines import CoxPHFitter

from .constants import SURVIVAL_DATA_HEADERS, ODDS_RATIO_THRESHOLD_LIST, ODDS_RATIO_TITLE


def calculate_binary_metrics(df: pandas.DataFrame, prs_column_name: str, trait_code: str):
    auc, auc_lci, auc_uci = calculate_auc(df[[prs_column_name]], df[trait_code])
    or_1sd, or_1sd_lci, or_1sd_uci = calculate_odds_ratio_1sd(df, trait_code, prs_column_name)

    metrics_dict = {'auc': {'value': auc, 'lci': auc_lci, 'uci': auc_uci},
                    'odds_ratio_1sd': {'value': or_1sd, 'lci': or_1sd_lci, 'uci': or_1sd_uci}}

    df_with_percentiles = add_percentile_data(df, prs_column_name)
    metrics_dict = calculate_odds_ratio_in_percentile(df_with_percentiles, trait_code, prs_column_name, metrics_dict)
    # Add conditional odds ratio
    if all(x in df.columns for x in ['age_at_first_assessment', 'sex']):
        if len(df['sex'].unique()) == 1:
            cor_1sd, cor_1sd_lci, cor_1sd_uci = calculate_conditional_odds_ratio_1sd(df, trait_code, prs_column_name,
                                                                                     ['age_at_first_assessment'])
        else:
            cor_1sd, cor_1sd_lci, cor_1sd_uci = calculate_conditional_odds_ratio_1sd(df, trait_code, prs_column_name,
                                                                                     ['age_at_first_assessment', 'sex'])
        metrics_dict['cond_odds_ratio_1sd'] = {'value': cor_1sd, 'lci': cor_1sd_lci, 'uci': cor_1sd_uci}

    if all(h in df.columns for h in SURVIVAL_DATA_HEADERS):
        hr_1sd, hr_1sd_lci, hr_1sd_uci = calculate_hazard_ratio_1sd(df, 'age_event', trait_code, prs_column_name)
        metrics_dict['hazard_ratio_1sd'] = {'value': hr_1sd, 'lci': hr_1sd_lci, 'uci': hr_1sd_uci}

        if all(x in df.columns for x in ['age_at_first_assessment', 'sex']):
            if len(df['sex'].unique()) == 1:
                chr_1sd, chr_1sd_lci, chr_1sd_uci = calculate_conditional_hazard_ratio_1sd(df, trait_code,
                                                                                           prs_column_name,
                                                                                           ['age_at_first_assessment'])
            else:
                chr_1sd, chr_1sd_lci, chr_1sd_uci = calculate_conditional_hazard_ratio_1sd(df, trait_code,
                                                                                           prs_column_name,
                                                                                           ['age_at_first_assessment',
                                                                                            'sex'])
            metrics_dict['cond_hazard_ratio_1sd'] = {'value': chr_1sd, 'lci': chr_1sd_lci, 'uci': chr_1sd_uci}

    return metrics_dict


def calculate_quantitative_metrics(df: pandas.DataFrame, prs_column_name: str,
                                   trait_code: str):
    rsq, rsq_lci, rsq_uci = calculate_rsq(df, prs_column_name, trait_code)
    beta_1sd, beta_1sd_lci, beta_1sd_uci = calculate_betas_1sd(df, prs_column_name, trait_code)
    metrics_dict = {'rsq': {'value': rsq, 'lci': rsq_lci, 'uci': rsq_uci},
                    'beta_1sd': {'value': beta_1sd, 'lci': beta_1sd_lci, 'uci': beta_1sd_uci}}
    return metrics_dict


def calculate_hazard_ratio_1sd(df: pandas.DataFrame, duration_col: str, event_col: str, prs_column_name: str,
                               covariate_cols: list = None):
    try:
        filter_to = [duration_col, event_col, prs_column_name]
        if covariate_cols is not None:
            filter_to.extend(covariate_cols)
        df = df[filter_to]
        cph = CoxPHFitter()
        hr = cph.fit(df, duration_col=duration_col, event_col=event_col)
        hazard_ratio_1sd = hr.summary['exp(coef)'][prs_column_name]
        hazard_ratio_1sd_lci = hr.summary['exp(coef) lower 95%'][prs_column_name]
        hazard_ratio_1sd_uci = hr.summary['exp(coef) upper 95%'][prs_column_name]
    except BaseException:
        hazard_ratio_1sd, hazard_ratio_1sd_lci, hazard_ratio_1sd_uci = numpy.nan, numpy.nan, numpy.nan

    return hazard_ratio_1sd, hazard_ratio_1sd_lci, hazard_ratio_1sd_uci


def add_percentile_data(df: pandas.DataFrame, prs_column_name: str):
    for thr in ODDS_RATIO_THRESHOLD_LIST:
        percentile = (100 - thr) / 100
        varname = prs_column_name + '_INtop{}'
        varname_med = prs_column_name + '_INtop{}med'
        varname_at = prs_column_name + '_ATtop{}'

        quant = list(df[prs_column_name].quantile([percentile]))
        quant_med = list(df[prs_column_name].quantile([0.40, 0.60, percentile]))
        quant_at = list(df[prs_column_name].quantile([(percentile - 0.005), (percentile + 0.005)]))

        df[varname.format(thr)] = numpy.where((df[prs_column_name] < quant[0]), 0, 1)
        df[varname_med.format(thr)] = numpy.where(
            (df[prs_column_name] >= quant_med[0]) & (df[prs_column_name] <= quant_med[1]), 0, numpy.nan)
        df[varname_med.format(thr)] = numpy.where(
            (df[prs_column_name] >= quant_med[2]), 1, df[varname_med.format(thr)])

        df[varname_at.format(thr)] = numpy.where((df[prs_column_name] < quant_at[0]), 0, numpy.nan)
        df[varname_at.format(thr)] = numpy.where(
            (df[prs_column_name] >= quant_at[0]) & (df[prs_column_name] <= quant_at[1]), 1,
            df[varname_at.format(thr)])
    return df


def calculate_odds_ratio_in_percentile(df: pandas.DataFrame, trait_code: str, prs_column_name: str, metrics_dict: dict):

    colpattern = prs_column_name + '_INtop{}'
    for t in ODDS_RATIO_THRESHOLD_LIST:
        column_name = colpattern.format(str(t))
        or_1sd, or_1sd_lci, or_1sd_uci = calculate_odds_ratio_1sd(df, trait_code, column_name)
        metrics_dict[ODDS_RATIO_TITLE.format(t)] = {'value': or_1sd, 'lci': or_1sd_lci, 'uci': or_1sd_uci}

    return metrics_dict


def calculate_odds_ratio_1sd(df: pandas.DataFrame, trait_code: str, prs_column_name: str):
    """Calculate the odds ratio and confidence intervals"""
    return _calculate_effect_size(df, trait_code, sm.families.Binomial(), lambda x: numpy.exp(x), prs_column_name)


def calculate_conditional_odds_ratio_1sd(df: pandas.DataFrame, trait_code: str, prs_column_name: str, covars: list):

    params = [prs_column_name] + covars
    return _calculate_effect_size(df, trait_code, sm.families.Binomial(), lambda x: numpy.exp(x), prs_column_name,
                                  f'{trait_code}~{"+".join(params)}')


def calculate_conditional_hazard_ratio_1sd(df: pandas.DataFrame, trait_code: str, prs_column_name: str, covars: list):
    return calculate_hazard_ratio_1sd(df, 'age_event', trait_code, prs_column_name, covars)


def calculate_betas_1sd(df: pandas.DataFrame, prs_column_name: str, trait_code: str,
                        family=sm.families.Gaussian()):
    """Calculate the betas and confidence intervals"""
    return _calculate_effect_size(df, trait_code, family, lambda x: x, prs_column_name)


def _calculate_effect_size(df, trait_code, family, transform, prs_column_name: str, formula: str = None):
    try:
        if formula is None:
            formula = f'{trait_code}~{prs_column_name}'
        model = smf.glm(formula=formula, data=df, family=family).fit()
        metric = transform(model.params[prs_column_name])
        lci = transform(model.conf_int()[0][prs_column_name])
        with numpy.errstate(over='ignore'):
            uci = transform(model.conf_int()[1][prs_column_name])
        return metric, lci, uci
    except BaseException:
        return numpy.nan, numpy.nan, numpy.nan


def calculate_auc(prs_df: pandas.DataFrame, pheno_series: pandas.Series):
    """
    Calculate AUC statistics and confidence intervals.
    """

    try:
        LR = LogisticRegression(solver='liblinear')
        f_m = LR.fit(prs_df, pheno_series)
        p_m = pandas.DataFrame(f_m.predict_proba(prs_df), index=prs_df.index)
        auc = roc_auc_score(pheno_series, p_m.iloc[:, 1])
        auc_lci, auc_uci = auc_ci(pheno_series, p_m.iloc[:, 1])
    except BaseException:
        auc, auc_lci, auc_uci = numpy.nan, numpy.nan, numpy.nan

    return auc, auc_lci, auc_uci


def auc_ci(actual: pandas.Series, pred: pandas.Series, ci: float = 0.95):
    """
    Computes AUC and confidence intervals using a fast implementation of delong (DOI: 10.1109/LSP.2014.2337313)
    """

    def midranks(w: numpy.ndarray):
        m = len(w)
        w = numpy.append(w, w[-1] + 1)
        t = numpy.zeros(m)
        i = 1
        while i <= m:
            j = i
            while w[j - 1] == w[i - 1]:
                j += 1
            for k in range(i, j):
                t[k - 1] = (i + j - 1) / 2
            i = j
        return t

    pred = pred.values
    actual = actual.values
    pred_sort = numpy.argsort(pred)
    pred = pred[pred_sort]
    actual = actual[pred_sort]
    pos = pred[actual == 1]
    neg = pred[actual == 0]
    npos = len(pos)
    nneg = len(neg)
    # the 'expensive' bit
    posranks = midranks(pos)
    negranks = midranks(neg)
    bothranks = midranks(pred)
    auc = bothranks[actual == 1].sum() / (npos * nneg) - (npos + 1.) / (2. * nneg)
    v01 = (bothranks[actual == 1] - posranks) / nneg
    v10 = 1. - (bothranks[actual == 0] - negranks) / npos
    var = (numpy.var(v01) / (npos - 1)) + (numpy.var(v10) / (nneg - 1))
    CI = stats.norm.ppf(
        [(1 - ci) / 2, (ci + (1 - ci) / 2)],
        loc=auc,
        scale=numpy.sqrt(var))
    CI = numpy.clip(CI, 0, 1)
    return CI[0], CI[1]


def rsq_ci(rsq: numpy.float, n: int, ci: float = 0.95):
    """Computes approximate confidence interval for R2, based on fishers transformation"""
    mu = numpy.log((1 + numpy.sqrt(rsq)) / (1 - numpy.sqrt(rsq)))
    delta = stats.norm.ppf(ci + (1 - ci) / 2) * numpy.sqrt(4 / n)
    limits = [mu - delta, mu + delta]
    CI = [((numpy.exp(li) - 1) / (numpy.exp(li) + 1))**2 for li in limits]
    CI = numpy.clip(CI, 0, 1)
    return CI[0], CI[1]


def calculate_rsq(df: pandas.DataFrame, prs_column_name: str, pheno_field: str):
    """ Calculate the R^2 values and confidence intervals"""

    funct = f'{pheno_field}~{prs_column_name}'
    try:
        model = smf.glm(formula=funct, data=df, family=sm.families.Gaussian()).fit()
        rsq = 1 - (model.deviance / model.null_deviance)
        if '+' not in funct:
            # approximation might not work if more than one predictor and correlated
            rsq_lci, rsq_uci = rsq_ci(rsq, df[[prs_column_name]].shape[0])
        else:
            rsq_lci, rsq_uci = numpy.nan, numpy.nan
    except BaseException:
        rsq, rsq_lci, rsq_uci = numpy.nan, numpy.nan, numpy.nan

    return rsq, rsq_lci, rsq_uci


def calculate_prs_sumstats(df: pandas.DataFrame, prs_column_name: str):

    mean = df[prs_column_name].mean()
    sem = df[prs_column_name].std() / numpy.sqrt(len(df[prs_column_name]))
    mean_lci = mean - sem * 1.96
    mean_uci = mean + sem * 1.96

    return {'mean': mean, 'mean_lci': mean_lci, 'mean_uci': mean_uci, 'std': df[prs_column_name].std()}


def calculate_sample_information(df: pandas.DataFrame):

    num_samples = len(df)
    if len(df[(df['prevalent'] == 1) & (df['incident'] == 1)]) > 0:
        raise AssertionError(f"Cannot have samples that are both prevalent and incident - "
                             f"{len(df[(df['prevalent'] == 1) & (df['incident'] == 1)])} samples found")
    num_prevalent = len(df[df['prevalent'] == 1])
    num_incident = len(df[df['incident'] == 1])
    num_cases = num_prevalent + num_incident
    return {'n_samples': num_samples, 'n_cases': num_cases, 'n_prevalent': num_prevalent,
            'n_incident': num_incident, 'n_controls': num_samples - num_cases}
