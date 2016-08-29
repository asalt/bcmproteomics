"""Script to run multiple pairwise comparisons"""
import itertools
from datetime import datetime
from collections import defaultdict
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

def _main(comparisons, ibaqnorm=None, tnormalize=None, desc='', seed=None):
    """Perform pairwise comparisons across multiple experiments. An average ranking
    of each gene is calculated
    """
    if tnormalize is not None:
        _get_genelists()  # load into memory once and only once

    up_genes, down_genes = defaultdict(lambda :
                                       defaultdict(list)), defaultdict(lambda :
                                                                       defaultdict(list))

    dtype_dict={'e2g_GeneID': object, 'GeneID': object}
    for exps in comparisons:
        exp_counter = len(comparisons)
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

        exp_join = ispec.join_exps(ctrl, treat,
                                   normalize=ibaqnorm, seed=seed)  # automatically does the machine learning
        if 'USD' not in exp_join.df.columns:
            continue

        df_U = exp_join.df[exp_join.df.USD=='U'].sort_values(by=['dlog_diBAQ','dPSMs'],ascending=[False,
                                                            False]).reset_index(drop=True)
        df_D = exp_join.df[exp_join.df.USD=='D'].sort_values(by=['dlog_diBAQ','dPSMs'],ascending=[True,
                                                            True]).reset_index(drop=True)
        df_U['ranks'] = df_U.index+1
        df_D['ranks'] = df_D.index+1
        for (ix1,U) , (ix2,D) in zip(df_U.iterrows(), df_D.iterrows()):
            up_genes[U.GeneID]['rank'].append(ix1+1)
            up_genes[U.GeneID]['desc'].append(desc)
            up_genes[U.GeneID]['prob'].append(U.USD_prob)
            up_genes[U.GeneID]['exps'].append('_'.join([str(ctrl.recno),
                                                        str(ctrl.runno),
                                                        str(treat.recno),
                                                        str(treat.runno),
            ]
            ))

            down_genes[D.GeneID]['rank'].append(ix2+1)
            down_genes[D.GeneID]['desc'].append(desc)
            down_genes[D.GeneID]['prob'].append(D.USD_prob)
            down_genes[D.GeneID]['exps'].append('_'.join([str(ctrl.recno),
                                                          str(ctrl.runno),
                                                          str(treat.recno),
                                                          str(treat.runno),
            ]
            ))

        df_U.index = df_U['GeneID'].astype('float64')
        df_D.index = df_D['GeneID'].astype('float64')
        df_UD = pd.concat([df_U, df_D])
        exp_join._df = exp_join._df.join(df_UD.ranks)
        #outname = '_'.join([x.strip(',') for x in exp_join.__repr__().split(' ')])
        #exp_join.df.to_csv(outpath+outname+'.tab', sep='\t', index=False,
        #                   columns=['GeneID','ranks','USD','USD_prob'])

        up_genes_df = pd.DataFrame({'GeneID': list(up_genes.keys())})
        down_genes_df = pd.DataFrame({'GeneID': list(down_genes.keys())})
        cols = ['GeneID','global_rank','avg_rank','count/total',
                'exp_print','desc','usd_prob','GeneSymbol','GeneDescription',
                'FunCats']
    if not up_genes_df.empty:
        up_df = _gene_summary(pd.DataFrame({'GeneID': list(up_genes.keys())}),
                                 up_genes, exp_counter=exp_counter)
    else:
        up_df = pd.DataFrame(columns=cols)
    if not down_genes_df.empty:
        down_df = _gene_summary(pd.DataFrame({'GeneID': list(down_genes.keys())}),
                               down_genes, exp_counter=exp_counter)
    else:
        down_df = pd.DataFrame(columns=cols)
    return (up_df, down_df)
    #u_name = expset+'_up_combined.tab'
    #d_name = expset+'_down_combined.tab'
    #up_df.to_csv(outpath+u_name, index=False, sep='\t', columns=cols, dtype=dtype_dict)
    #down_df.to_csv(outpath+d_name, index=False, sep='\t', columns=cols, dtype=dtype_dict)

def _expconstructor(ctrls=None, samples=None, by_pairs=False):
    """ Magically makes all the data numbers you wish to find
    """
    ctrls = [(x,1,1) if not isinstance(x,tuple) else x for x in ctrls]
    samples = [(x,1,1) if not isinstance(x,tuple) else x for x in samples]
    ctrl_e2gs = []
    conn = ispec.filedb_connect()
    if isinstance(conn, str):
        print('Unable to connect to iSPEC')
        print(conn)
        return None
    print('Loading data from iSPEC, please wait...')
    for ctrl in ctrls:
        exp = ispec.E2G(*ctrl)
        if len(exp.df) != 0: # only pick exps that have data!
            ctrl_e2gs.append(exp)
        elif len(exp.df) == 0:
            print('No data for ', exp, '\nSkipping')
    sample_e2gs = []
    for sample in samples:
        exp = ispec.E2G(*sample)
        if len(exp.df) != 0: # only pick exps that have data!
            sample_e2gs.append(exp)
        elif len(exp.df) == 0:
            print('No data for ', exp, '\nSkipping')
    if by_pairs:
        return [(ctrl, sample) for ctrl, sample in zip(ctrl_e2gs, sample_e2gs)]

    pairs = [(ctrl, sample) for ctrl, sample in itertools.product(ctrl_e2gs, sample_e2gs)]
    return pairs

def multicomparison(ctrls=None, samples=None, description=None, ibaq_normalize=None,
                    taxon_normalize=None, seed=None, by_pairs=False):
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


    pairs = _expconstructor(ctrls, samples, by_pairs=by_pairs)
    if not pairs:
        print('No pairs of experiments to compare!')
        return

    if description is None:
        description = datetime.ctime(datetime.now())

    up_df, down_df = _main(pairs, tnormalize=taxon_normalize, desc=description, seed=seed)
    return (up_df, down_df)
