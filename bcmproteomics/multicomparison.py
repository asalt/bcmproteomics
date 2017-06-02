"""Script to run multiple pairwise comparisons"""
import os
import sys
import itertools
from datetime import datetime
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd
from bcmproteomics import ispec

def _get_genelists():
    """Need to grab these lists of genes once and only once
    """
    global hu_genes, mou_genes
    hu_genes = ispec.get_geneids(9606)
    mou_genes = ispec.get_geneids(10090)

def _gid_rank_info(geneid, d, n=1):
    """takes in a given geneid and a dictionary with a list of ranks for each geneid and returns
    useful information"""
    values = d[geneid]['rank']
    exp_desc = ', '.join(set(d[geneid]['desc']))
    usd_prob = '{:.2f}'.format(d[geneid]['prob'][0])  # always length 1
    count_tot = len(values)/n
    #count_tot = '{}/{}'.format(len(values),n)  # times 2 since we are looking at pairwise exps.
    exp_print = '__'.join(d[geneid]['exps'])
    return np.mean(values), count_tot, exp_desc, usd_prob, exp_print

def _gene_summary(df, d, exp_counter=1):
    """Returns average rank, ratio of count/total, description,
    usd probability, and experiment print
    """
    df['avg_rank'], df['count/total'], df['desc'], df['usd_prob'], df['exp_print']= \
                    list(zip(*df.apply(lambda x : _gid_rank_info(x['GeneID'], d, n=exp_counter), axis=1)))

    df.sort_values(by=['count/total','avg_rank'],ascending=[False,True], inplace=True)
    funcat_table = ispec.get_funcats(df.GeneID.tolist())
    df = pd.merge(df,funcat_table, how='inner', on='GeneID', copy=False)  # add the funcats
    df['FunCats'] = df['FunCats'].fillna('')
    return df

def _taxon_normalizer(df, ratio):
    global hu_genes, mou_genes
    area = 'iBAQ_dstrAdj'
    #ratio = tnormalize[df['_e2g_EXPRecNo'].values[0]]

    if df['GeneID'] in hu_genes:
        tid = 9606
    elif df['GeneID'] in mou_genes:
        tid = 10090
    else:
        tid = 0
    r = ratio.get(tid,1)
    if r == 0:
        r = 10**-3  # avoid divide by zero error
    return df[area]/r

def format_result(df, name=None):
    """
    """
    if name is None:
        name ='<test name>'

    def aggregate(grp, all_comps):
        usd_stamp = list()
        usd_stamp_prob = list()
        fold_change = list()
        psms = list()
        idset = list()
        for comp in all_comps:
            if comp not in grp['comparison'].values:
                usd_stamp.append('-')
                usd_stamp_prob.append('-')
                fold_change.append(0)
                psms.append(0)
                idset.append(0)
                continue
            data = grp[grp['comparison'] == comp].drop_duplicates(subset=['GeneID']).squeeze()
            usd_stamp.append(data['USD'])
            usd_stamp_prob.append(data['USD_prob'])
            fold_change.append(data['dlog_diBAQ'])
            psms.append(data['dPSMs'])
            idset.append(data['dIDSet'])
        return pd.Series({'USD_Stamp' : ''.join(map(str, usd_stamp)),
                          'USD_Probabilities' : '|'.join(map(str, usd_stamp_prob)),
                          'USD_dPSMs' : '|'.join(map(str, psms)),
                          'USD_dIDSet' : '|'.join(map(str, idset)),
                          'USD_iBAQ_Log10_FoldChanges'   : '|'.join(map(str,fold_change))}
        )

    def average_subset(row):
        class_ = row['USD_Class']
        if class_ == 'N':
            return pd.Series({'USD_Class_iBAQ_log10_FoldChange':0,
                              'USD_Class_dPSMs' : 0,
                              'USD_Class_dIDSet' : 0,
                              'USD_Class_Probability':0})
        indices = [index for index, value in enumerate(row['USD_Stamp']) if value == class_]
        results = dict()
        for col in ('USD_iBAQ_Log10_FoldChanges', 'USD_Probabilities', 'USD_dPSMs', 'USD_dIDSet'):
            changes = [float(x) for ix, x in enumerate(row[col].split('|'))
                       if ix in indices]
            results['avg_'+col] = np.mean(changes)
        return pd.Series({'USD_Class_iBAQ_log10_FoldChange': results['avg_USD_iBAQ_Log10_FoldChanges'],
                          'USD_Class_dPSMs' : results['avg_USD_dPSMs'],
                          'USD_Class_dIDSet' : results['avg_USD_dIDSet'],
                          'USD_Class_Probability': results['avg_USD_Probabilities']
        })

    def rank_class(df):
        """Rank on per class basis
        """
        rank_df = df.assign(ibaq_abs = lambda x: np.abs(x['USD_Class_iBAQ_log10_FoldChange']))
        rank_on = ['USD_Class_Probability', 'ibaq_abs', 'USD_dIDSet', 'USD_dPSMs']

        return (rank_df.sort_values(by=rank_on, ascending=False)
                .groupby('USD_Class')
                .cumcount()+1)




    all_comparisons = df['comparison'].unique()

    out = df.groupby('GeneID').apply(aggregate, all_comparisons)
    out['Test_Comparisons'] = '|'.join(all_comparisons)
    out['Test_Count'] = len(all_comparisons)
    out['Test_Name'] = name

    out['USD_Class'] = 'N'
    for s in 'USD':
        col = 'USD_{}_Count'.format(s)
        out[col] = out.USD_Stamp.str.count(s)
        agree_col = 'USD_{}_Agreement'.format(s)
        out[agree_col] = out[col] / out['Test_Count']
        out.loc[out[agree_col]>=.75, 'USD_Class'] = s

    avgs = out.apply(average_subset, axis=1)
    out_avgs = out.join(avgs)
    out_avgs['USD_Class_Rank'] = 0
    for s in 'UD':
        out_avgs.loc[out_avgs['USD_Class'] == s, 'USD_Class_Rank'] = rank_class(out_avgs.query('USD_Class == @s'))

    return out_avgs

    # probability and average fold change postionally

def _get_meta(grp):
    keys = ('GeneDescription', 'GeneSymbol')
    results = dict()
    for key in keys:
        try:
            result = grp[key].unique()[0]
        except IndexError:
            result = ''
        results[key] = result
    return pd.Series(results)

def _main(comparisons, ibaqnorm=None, tnormalize=None, desc='', seed=None, name=None):
    """Perform pairwise comparisons across multiple experiments. An average ranking
    of each gene is calculated
    """
    if tnormalize is not None:
        _get_genelists()  # load into memory once and only once

    # up_genes, down_genes = defaultdict(lambda :
    #                                    defaultdict(list)), defaultdict(lambda :
    #                                                                    defaultdict(list))

    dtype_dict={'e2g_GeneID': object, 'GeneID': object}

    results = list()
    exp_counter = len(comparisons)
    for exps in comparisons:
        ctrl = exps[0]
        treat = exps[1]
        print('Processing {} vs {}'.format(ctrl, treat))

        if tnormalize:
            for exp in [ctrl, treat]:
                exp_ratio = tnormalize[exp.recno]
                exp.df['ibaq_norm'] = exp.df.apply(lambda x: _taxon_normalizer(x, exp_ratio),
                                                   axis=1)
                exp.df.rename(columns={'iBAQ_dstrAdj':'iBAQ_dstrAdj_old'}, inplace=True)
                exp.df.rename(columns={'ibaq_norm':'iBAQ_dstrAdj'}, inplace=True)
        f = open(os.devnull, 'w')
        stdout = sys.stdout
        sys.stdout = f
        exp_join = ispec.join_exps(ctrl, treat, normalize=ibaqnorm, seed=seed)  # automatically does the machine learning
        sys.stdout = stdout
        f.close()

        repr_ = '{!r}:{!r}'.format(treat, ctrl)
        COLS = ['USD', 'USD_prob', 'dlog_diBAQ', 'GeneSymbol', 'GeneDescription',
                'dPSMs', 'dIDSet']
        result = exp_join.df[COLS].copy()
        result['comparison'] = repr_
        if 'GeneID' not in result.columns:
            result['GeneID'] = result.index
        results.append(result)
    df = (pd.concat(results)
          .reset_index()
          .where(lambda x: ~x['GeneID'].isnull()).dropna(how='all')
          # .query('~GeneID.isnull()')  # not working for some reason on some computers
          .assign(GeneID = lambda x : x['GeneID'].astype(int)))

    metadata = df.groupby('GeneID').apply(_get_meta)

    print('Formatting results...')
    formatted_df = format_result(df, name=name).join(metadata)
    ctrls_list = '|'.join( sorted(set(repr(x[0]) for x in comparisons)) )
    treat_list = '|'.join( sorted(set(repr(x[1]) for x in comparisons)) )
    formatted_df['Test_Controls'] = ctrls_list
    formatted_df['Test_Experiments']     = treat_list
    formatted_df['Test_Date']     = datetime.ctime(datetime.now())
    formatted_df.index = formatted_df.index.astype(int)

    return formatted_df

def _expconstructor(ctrls=None, samples=None, by_pairs=False, data_dir=None):
    """ Magically makes all the data numbers you wish to find
    """
    e2gs = ispec.async_get_datas(ctrls+samples, verbosity=1)
    ctrl_e2gs = filter(None, [x if len(x.df) > 0 else print('No data for', x, 'skipping..')
                 for x in e2gs[0:len(ctrls)]])
    sample_e2gs = filter(None, [x if len(x.df) > 0 else print('No data for', x, 'skipping..')
                   for x in e2gs[len(ctrls):]])

    # ctrls = [(x,1,1) if not isinstance(x,tuple) else x for x in ctrls]
    # samples = [(x,1,1) if not isinstance(x,tuple) else x for x in samples]
    # ctrl_e2gs = []
    # print('Loading data from iSPEC, please wait...')
    # for ctrl in ctrls:
    #     exp = ispec.E2G(*ctrl, data_dir=data_dir)
    #     if len(exp.df) != 0: # only pick exps that have data!
    #         ctrl_e2gs.append(exp)
    #     elif len(exp.df) == 0:
    #         print('No data for ', exp, '\nSkipping')
    # sample_e2gs = []
    # for sample in samples:
    #     exp = ispec.E2G(*sample, data_dir=data_dir)
    #     if len(exp.df) != 0: # only pick exps that have data!
    #         sample_e2gs.append(exp)
    #     elif len(exp.df) == 0:
    #         print('No data for ', exp, '\nSkipping')
    if by_pairs:
        return [(ctrl, sample) for ctrl, sample in zip(ctrl_e2gs, sample_e2gs)]

    pairs = [(ctrl, sample) for ctrl, sample in itertools.product(ctrl_e2gs, sample_e2gs)]
    return pairs

def multicomparison(ctrls=None, samples=None, description=None, ibaq_normalize=None,
                    name=None,
                    taxon_normalize=None, seed=None, by_pairs=False, data_dir=None):
    """

    Function to run combinations of pairwise comparisons of experimental results located
    in an iSPEC database. Returns two pandas dataframes, the first containing all gene products
    found be considered up in the second group at least once and the second containing all
    gene products found to be down at least once.

    :param ctrls: list of experiment records that belong to control group
                  format is [ (rec1, run1, search1), (rec2, run2, search2) ]
    :param samples: list of experiment records that belong to sample group
                    same format as ctrls
    :param description: (optional) description of the comparison. Defaults to today's date.
    :param ibaq_normalize: 'mean', 'median' or None (Optional)
                           method of normalization of iBAQ_dstrAdj
                           Has the side effect of adding ibaq_norm column to each input DataFrame.
    :param taxon_normalize: (optional) a dictionary which contains dictionaries of different taxon ratios
                            for a given experiment
    :param seed: (optional) set seed for random forest classifier (default None)
    :param by_pairs: (optional) only do comparisons by pairs based on the ordering of the experiments
    :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
    :returns: up_df, down_df
    :rtype: pd.DataFrame, pd.DataFrame


    :Example:

    >>> ctrls = [(12345,1,1), (12345,2,1)]
    >>> treatments = [12346, 12347]
    >>> up, down = multicomparison(ctrls, treatments, description='My Comparison 1')
    .. note::
        The tuple format is only necessary when specifying run and search numbers.
        The record number can be the only entry and the run and search numbers will default to 1.

    :Example:

    >>> taxon_normalize : {12345: {9606: 0.4, 10090: 0.6},
                          {12346: ... }
    >>> up, down = multicomparison(ctrls, treatments, description='My Comparison 2',
                                   taxon_normalize=taxon_normalize)

    .. note::
        Here experiment number 12345 has a human ratio of 0.4 and a mouse ratio of 0.6.
        If taxon_normalize is used, an entry for each experiment record must be included.
    """

    if not ctrls and samples:
        print('Both input lists must have at least 1 experiment!')
        return


    pairs = _expconstructor(ctrls, samples, by_pairs=by_pairs, data_dir=data_dir)
    if not pairs:
        print('No pairs of experiments to compare!')
        return

    if description is None:
        description = datetime.ctime(datetime.now())

    # up_df, down_df = _main(pairs, tnormalize=taxon_normalize, desc=description, seed=seed)
    # return (up_df, down_df)

    result = _main(pairs, tnormalize=taxon_normalize, desc=description, seed=seed, name=name)

    COLS = [('Test_Name', )]

    COL_MAPPING = OrderedDict()


    COL_MAPPING = {'m'}
    return result
