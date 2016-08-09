"""utilities for wrangling data"""
import re
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from bcmproteomics import ispec
# scratch

def aggregate_all_data(exps, cols=None):
    """
    Function for aggregating all data across a list of E2G experiments"""
    funcat_cols = ['GeneCapacity', 'GeneSymbol', 'GeneDescription', 'FunCats', 'GeneID']
    if cols is None:
        cols = [x for x in exps[0].df.columns if x not in funcat_cols + ['PeptidePrint']]
        cols.append('Sample')
        cols = list(set(cols))
    for exp in exps:
        if 'Sample' not in exp.df.columns:
            exp.df['Sample'] = repr(exp)
    df = pd.concat( [exp.df[ cols ] for exp in exps] )
    pivot_df = df.pivot(columns='Sample')
    funcats = ispec.get_funcats( df.index.tolist())
    return pivot_df.join(funcats)

def log_normalize(data, cols=None, how='log2', rename=False):
    """log normalize columns in a DataFrame or array, replacing -inf with 0 and shifting everything up by the min value
    how: how to normalize - log2, log10, or log1p
    cols is a list of columns to normalize. If left as None will go across all values
    rename : if True, add a prefix to each column with the transformation
    """
    if cols is None and isinstance(data, pd.DataFrame):
        cols = data.columns
    if how == 'log2':
        func = np.log2
    elif how == 'log10':
        func = np.log10
    elif how == 'log1p':
        func = np.log1p
    if isinstance(data, pd.DataFrame):
        out = data[cols].fillna(0).apply(func)
    else:
        out = func(data)
    minvalue = np.min(np.array(out)[np.array(out) != -np.inf])
    maxvalue = np.max(np.array(out)[np.array(out) != -np.inf])
    min_scale = 1.1 if minvalue < 0 else 0.9
    max_scale = 1.1 if minvalue > 0 else 0.9
    if minvalue == 0:
        minvalue = func(10**-3)
    if maxvalue == 0:
        maxvalue = func(10**3)
    # out = out + abs(minvalue)
    if isinstance(data, pd.DataFrame):
        out.replace(-np.inf, minvalue*min_scale, inplace=True)
        out.replace(np.inf, maxvalue*max_scale, inplace=True)
    else:
        np.maximum(out, minvalue * min_scale, out)
        np.minimum(out, maxvalue * max_scale, out)
    if rename and isinstance(data, pd.DataFrame):
        out.rename( columns={ x : '{}_{}'.format(how, x) for x in cols}, inplace=True)
    return out

def aggregate_data(exps, area_col='iBAQ_dstrAdj', normalize=None, tofilter=None,
                   label_by=None, funcats=False):
    """
    Aggregate data based on 1 column (default iBAQ_dstrAdj)
    Input is a list of E2G experiments.

    normalize : How to normalize the data. Acceptable values are None, 'log2', 'log10', log1p, and 'row'.
    to_filter : How to filter the data for each experiment. Acceptable values are
        None, 'strict', 'set1' ('set1' is strict and set1 only)
    label_by  : if set to a string that is an attribute, will rename each column appropriately
    funcats   : if set to True includes funcats for each gene
    Output is a joined DataFrame of an area column of your choice. (Default iBAQ_dstrAdj).
    Area column is row normalized by default."""
    for exp in exps: # for identification
        exp.df['Sample'] = repr(exp)
    if tofilter is None:
        df = pd.concat( [exp.df[['Sample', area_col]] for exp in exps] )
    elif tofilter == 'strict':
        df = pd.concat( [exp.strict()[['Sample', area_col]] for exp in exps] )
    elif tofilter == 'set1':
        df = pd.concat( [exp.strict(set1=True)[['Sample', area_col]] for exp in exps] )
    else:
        df = pd.concat( [exp.df[exp.df.FunCats.str.contains(tofilter)][['Sample',
                                                                area_col]] for exp in exps] )
    pivot_df = df.pivot(columns='Sample', values=area_col)

    if normalize is None:
        norm_df = pivot_df
    elif normalize == 'row':
        norm_df =  pivot_df.div(data.sum(axis=1), axis=0).dropna(axis=0, how='all')
    elif any(normalize == x for x in ('log2', 'log10', 'log1p')):
        if normalize == 'log2':
            norm_df = log_normalize(pivot_df)
        elif normalize == 'log10':
            norm_df = log_normalize(pivot_df, how='log10')
        elif normalize == 'log1p':
            norm_df = log_normalize(pivot_df, how='log1p')

    if label_by is not None:
        if not all([label_by in exp.__dict__ for exp in exps]) == True:
            raise AttributeError ("'{}' is an invalid attribute. "
                                  "Valid attributes are : {}".format(label_by,
                                                                     '\n'.join([x for x in exp.__dict__ if not x.startswith('_')])))
        to_rename = dict()
        for col in norm_df.columns:
            exp = [exp for exp in exps if repr(exp) == col][0]  # works since col == repr(exp)
            to_rename[col] = '{} {}'.format(col, exp.__getattribute__(label_by))

        norm_df.rename(columns=to_rename, inplace=True)
        pat = re.compile(r'(\d+_\d+_\d+)')
        norm_df = norm_df.reindex_axis(sorted(norm_df.columns,
                                              key=lambda x: pat.split(x)[::-1]), axis=1)
        # norm_df = norm_df[ [sorted(norm_df.columns, key=lambda x: (x[-1], x[-2]))]]

    if funcats:
        funcats = ispec.get_funcats(norm_df.index.tolist())
        norm_df = norm_df.join(funcats)

    return norm_df

def t_test(set1, set2, equal_var=True, alpha=0.05, multiple_correction=True, correction_type=None,
           area_col='iBAQ_dstrAdj', tofilter=None):
    """Input is two lists of E2G experiments to be compared
    Returns a DataFrame with statistical information for each geneid"""
    if correction_type is None:
        correction_type = 'fdr_tsbky'

    set1_len = len(set1)
    set2_len = len(set2)
    df = aggregate_data(set1+set2, area_col=area_col, tofilter=tofilter).fillna(0)
    df = df.reindex_axis( [repr(x) for x in set1+set2], axis=1)
    sample1 = df.columns[0:set1_len]
    sample2 = df.columns[set1_len::]
    df['mean1'] = df[sample1].mean(axis=1)
    df['mean2'] = df[sample2].mean(axis=1)
    df['std1'] = df[sample1].std(axis=1)
    df['std2'] = df[sample2].std(axis=1)
    df['fold_change'] = np.divide(df['mean2'], df['mean1'])
    numeric_max = df.fold_change.replace(np.inf,0).max()  # max that isn't infinity
    df['fold_change'] = df['fold_change'].replace(np.inf, numeric_max+.05*numeric_max)
    df['fold_change'].fillna(0, inplace=True)

    df['t_stat'], \
    df['p_value'] = \
                    list(zip(*df.apply(lambda x : stats.ttest_ind(x[sample1],
                                                                  x[sample2],
                                                                  equal_var=equal_var),
                    axis=1)))
    if multiple_correction:
        fdr_out = multipletests(df.p_value.fillna(1), alpha=alpha, method=correction_type)
        df['fdr_pass'] = fdr_out[0]
    else:
        df['fdr_pass'] = True
    return df

def volcanoplot(data, ax=None):
    """A nice way to make a volcano plot from the result given by t_test"""
    log2_fc = np.log2(data['fold_change'])
    numeric_min = log2_fc.replace(-np.inf, 0).min()
    log2_fc = log2_fc.replace(-np.inf, numeric_min*1.05)
    colors = ['#19ff19' if (abs(x)>2 and y<=0.05 and z) else '#ff1919'
                    for x,y,z in zip(log2_fc,
                                     data['p_value'],
                                     data['fdr_pass'])]
    if ax is None:
        fig, ax = plt.subplots()
    ax.scatter(log2_fc, -np.log10(data['p_value']).fillna(0), c=colors,
               alpha=.8)
    ax.set_xlabel('log$_2$ fold change')
    ax.set_ylabel('-log$_{10}$ p-value')
    return fig, ax




def make_heatmap(data, gene_ids=None, ax=None, order=None, cmap='YlOrRd', square=False):
    """A nice way to make heatmaps from the result given by aggregate_data"""
    pat = re.compile('\d+_\d+_\d+')
    if gene_ids is not None:
        mydata = data.loc[gene_ids][[x for x in data.columns if pat.match(x)]]
    else:
        mydata = data[[x for x in data.columns if pat.match(x)]]

    if ax is None:
        fig, ax = plt.subplots()

    ax = sb.heatmap(mydata, ax=ax, cmap=cmap,
                    linewidths=0.0, rasterized=True, square=square,)
    for tick in ax.axes.yaxis.get_ticklabels():
        tick.set_rotation(0)

    fig.autofmt_xdate()

    return fig, ax
