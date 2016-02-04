"""Docstring
"""
from os import path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
def score_experiments(exp1_2, seed=None):
    """Classify gene products from joined E2G files
    based on supervised machine learning algorithm random forest"""

    clf = classifier_training(merge_cats=True, subsample=True, seed=seed)
    features, df_all = feature_grabber(exp1_2.df)
    df_all['USD'], df_all['USD_prob'] = clf.predict(features),\
                                        np.max(clf.predict_proba(features), axis=1,
                                               keepdims=False)

    colormap = {'DD':' #ff1919', 'D': '#b20000', 'S':' #333333',
                'SS': '#000000','U': '#00b200', 'UU': '#19ff19'}
    df_all['color'] = df_all['USD'].map(lambda x : colormap[x])
    exp1_2._df = df_all
    return exp1_2


def feature_grabber(df):
    """gets and adds features to a merged dataframe.
    The merged dataframe columns have default _x and _y extentions
    """
    #df = e2g_col_renamer(df)
    def normalize_pepts(df, col):
        old_value = df[col]
        max_value = max(df.get('GeneCapacity_x', 1), df.get('GeneCapacity_y', 1),
                        df.get('GeneCapacity',1), 1) # for backwards compatability
        return old_value/max_value
    # -----------------------Special replace for a certain column ---------------------------------#

    #df.replace(to_replace=replace_dict, inplace=True)
    #df[df['IDGroup_x']==0]['IDGroup_x']=9
    #df['IDGroup_x'].replace(0,9, inplace=True)
    # 1 is the best IDGroup. If no IDGroup (and thus no gene, then we need to make it the worst IDGroup, 9
    #df['IDGroup_y'].replace(0,9, inplace=True)
    #df['IDSet_x'].replace(0,4, inplace=True)  # new IDSet, 4, means doesn't exist
    #df['IDSet_y'].replace(0,4, inplace=True)
    df.fillna(0, inplace=True)
    #df['iBAQ_dstrAdj_x'] =  df['iBAQ_dstrAdj_x'].fillna(0)
    #df['iBAQ_dstrAdj_y'] = df['iBAQ_dstrAdj_y'].fillna(0)
    # ---------------------------------------------------------------------------------------------#
    feature_cols = ['PSMs', 'IDSet', 'IDGroup', 'PeptideCount_S',
                    'PeptideCount_u2g', 'iBAQ_dstrAdj']  # modify here to change features
    replace_dict = {'IDGroup': 9, 'IDSet': 4,}

    for col in replace_dict:
        replace_value = replace_dict[col]
        for v in ['_x','_y']:
            df.loc[df[col+v] == 0, col+v] = replace_value

    for col in feature_cols:
        #print(col)
        if col in ['iBAQ_dstrAdj']:
            df['d'+col] = (df[col+'_y']/df[col+'_x'])
        elif col not in ['IDGroup', 'IDSet']:
            dcol = 'd'+col
            df[dcol] = (df[col+'_y']-df[col+'_x'])
            df[dcol] = df.apply(normalize_pepts, args=(dcol,), axis=1)
        elif col in ['IDGroup', 'IDSet']:
            #pass
            df['d'+col] = np.exp(1/df[col+'_y']) - np.exp(1/df[col+'_x'])

    #  log of the change in the iBAQ
    df['dlog_diBAQ'] = np.log2(df.diBAQ_dstrAdj)
    #print(df['dlog_diBAQ'])
    localmin = df.dlog_diBAQ.replace(-np.inf, 0).min()
    localmax = df.dlog_diBAQ.replace(np.inf, 0).max()
    if localmin == 0:
        localmin = np.log2(-10**3)
    if localmax == 0:
        localmax = np.log2(10**3)

    df['dlog_diBAQ'].replace(-np.inf, localmin*1.1, inplace=True)
    df['dlog_diBAQ'].replace(np.inf, localmax*1.1, inplace=True)
    feature_cols[-1] = 'log_diBAQ'
    #  log of the change in the iBAQ
    #feature_cols.remove('e2g_IDGroup')
    df.replace(to_replace=[np.inf], value=1e3, inplace=True)
    df.fillna(0, inplace=True)  # 0/0 == NaN, for example
    features = np.array(df[['d'+col for col in feature_cols]])
    return features, df


def e2g_col_renamer(df, returndf=True):
    ''' takes in a dataframe from an e2g file and returns the
    same dataframe with the e2g_ extension removed from many
    of the columns.

    Requires an unmodified (not merged, joined, etc..) dataframe'''
    df.rename(columns={k: k.split('e2g_')[1] for k in
                       [e2gcol for e2gcol in df.columns if e2gcol.startswith('e2g_')]},
              inplace=True)
    df.rename(columns={'n_iBAQ_dstrAdj': 'iBAQ_dstrAdj',
                       'n_iBAQ_dstrAdj_x': 'iBAQ_dstrAdj_x',
                       'n_iBAQ_dstrAdj_y': 'iBAQ_dstrAdj_y',
    },
              inplace=True)
    if returndf:
        return df

def get_training_data(merge_cats=True, subsample=True):
    ''' get training data'''
    training_dir = path.join(path.dirname(__file__), 'training_data')


    if not merge_cats:
        subsample = False  # subsample cannot be true if merge_cats is not true
    pairwise_training = [('12452_1_1_1_none_nolabel_e2g.txt', '12453_1_1_1_none_nolabel_e2g.txt'),
                         ('30223_1_1_1_none_nolabel_e2g.txt', '30224_1_1_1_none_nolabel_e2g.txt'),
                         ('30100_1_1_1_none_nolabel_e2g.txt', '30101_1_1_1_none_nolabel_e2g.txt'),
    ]
    pairwise_classes = ['gene2usd_PMA_12452_12453_run1.txt',
                        '30223_30224_USD.txt',
                        '30100_30101_USD.txt',
                        ]
    # All of this is data wrangling
    dfs = []
    for comparison in zip(pairwise_training, pairwise_classes):

        (f1, f2), train_data = comparison
        df1 = pd.read_table(path.join(training_dir, f1))
        df2 = pd.read_table(path.join(training_dir, f2))
        df1 = e2g_col_renamer(df1)
        df2 = e2g_col_renamer(df2)

        df_labels = pd.read_table(path.join(training_dir, train_data))
        df_labels['USD'] = df_labels.USD.astype('category', categories=['DD','D','S','SS','U','UU'],
                                                ordered=True)
        if merge_cats:
            for dup in ['DD', 'SS', 'UU']:
                df_labels['USD'] = df_labels['USD'].replace(dup, dup[0])
            # get rid of redundant categories, as redundant categories are redundant
        df1_2 = pd.merge(df1, df2, on='GeneID', how='outer')
        df1_2.rename(columns={'e2g_GeneID':'GeneID'}, inplace=True)
        temp_df = pd.merge(df1_2, df_labels, on='GeneID')  # only has columns that have assignment
        dfs.append(temp_df)

    df = pd.concat(dfs)
    #    feature_cols = ['e2g_PSMs', 'e2g_IDSet','e2g_IDGroup', 'e2g_PeptideCount_S', 'e2g_PeptideCount_u2g', \
    #            'e2g_n_iBAQ_dstrAdj']  # modify here to change your features
    if subsample:
        sample_val = max(len(df[df.USD == 'D']), len(df[df.USD == 'U']))
        temp_df = df[df.USD == 'S'].sample(n=sample_val, random_state=42)
        temp_df = pd.concat([temp_df, df[df.USD.isin(['D', 'U'])]])
        df = temp_df

    features, df = feature_grabber(df)

    labels = np.array(df['USD'])
    if not merge_cats:
        df['USD'] = df.USD.astype('category', categories=['DD', 'D', 'S', 'SS', 'U', 'UU'],
                                  ordered=True)
    elif merge_cats:
        df['USD'] = df.USD.astype('category', categories=['D', 'S', 'U'],
                                  ordered=True)

    # Data wrangling complete
    return features, labels, df

def get_classifier(seed=None):
    clf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='entropy',
            max_depth=6, max_features=None, max_leaf_nodes=None,
            min_samples_leaf=1, min_samples_split=6,
            min_weight_fraction_leaf=0.0, n_estimators=10, n_jobs=1,
            oob_score=False, random_state=seed, verbose=0,
                                 warm_start=False)
    return clf

def classifier_training(merge_cats=True, subsample=True, seed=None):

    clf = get_classifier(seed=seed)
    features, labels, df = get_training_data(merge_cats=merge_cats, subsample=subsample)
    clf.fit(features, labels)
    return clf


def plot_confusion_matrix(cm, labels, title='Confusion matrix', cmap=plt.cm.Blues, savefig=False, figname=None,
                          **kwargs):
    """ Takes in your n by n dimensional confusion matrix, the appropriate labels, 
    and optional title and color map.
    Default uses matplotlib imported as plt"""
    
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    print(labels)
    plt.title(title)
    cbar = plt.colorbar()
    cbar.solids.set_edgecolor("face")  # alternatively use rasterized=True somewhere
    tick_marks = np.arange(len(labels))
    plt.xticks(tick_marks, labels, rotation=45)
    plt.yticks(tick_marks, labels)
    #plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tick_params(axis='both',
                    bottom='off',
                    top='off',
                    left='off',
                    right='off',
                    )
    if savefig:
            plt.savefig(figname, **kwargs)
            print('Saving figure : ', figname)
    else:
        plt.show()

def _feature_importance(forest, features=None, plot=False, savefig=False, figname=None, **kwargs):
    """Takes in a trained forest classifier.
    Returns feature importances and makes a nice plot if requested.
    -----

    Taken from :
    http://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html
    """
    importances = forest.feature_importances_  #returns NotFittedError if not fitted yet
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)
    if features is not None:
        labels = [features[i] for i in indices]
    else:
        labels = indices
    # Print the feature ranking
    print("Feature ranking:")
    for f in range(X.shape[1]):
        print("{}. feature {} ({:.4f})".format(f + 1, labels[f], importances[indices[f]]))
    if plot:
        cmap = cm.Set1
        c = [cmap(n/12) for n in range(1,12)]
        #c.reverse()
        fig, ax = plt.subplots()
        ax.set_title = ("Feature import ances")
        ax.barh(range(X.shape[1]), importances[indices],
                color=c[2], xerr=std[indices], align='center',
                alpha=0.5)
        ax.set_yticks(range(X.shape[1]))
        ax.set_yticklabels(labels)
        ax.set_ylim([-1, X.shape[1]])
        ax.set_xlabel('Feature importance')

        xlim = ax.get_xlim()
        ax.set_xlim([0,xlim[1]]) # force xmin to be 0.0
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='x', direction='out', width=1) # set width equal to spine width

        ax.yaxis.set_ticks_position('none')
        ax.yaxis.labelpad = 20
        ax.tick_params(axis='y', pad=15)  # so not too close to spine
        
        for spine in ['top', 'right',]:
            ax.spines[spine].set_visible(False)
        if savefig:
            plt.savefig(figname, **kwargs)

def _average_class_accuracy(clf,features,labels, n=None, plot=False, seed=None, cats=None, savefig=False, figname=None, **kwargs):
    '''calculates and prints an average confusion matrix from the given trained classifier clf'''
    if n is None:
        n = 200
    if seed is not None:
        print('Warning, seed is fixed, results from each iteration will be identical.')

    cm_matrices, tot_acc = [], []
    for _ in range(n):
        features_train, features_test, labels_train, labels_test = \
        train_test_split(features, labels, test_size=0.33, random_state=seed)
        clf.fit(features_train, labels_train)
        y_true, y_pred = labels_test, clf.predict(features_test)
        # make a normalized confusion matrix
        cm = confusion_matrix(y_true, y_pred)
        cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        #print(x,cm_normalized)
        cm_matrices.append(cm_normalized)
        tot_acc.append(clf.score(features_test, labels_test))
    cm_normalized = np.mean(np.array(cm_matrices), axis=0)
    mean_acc = np.mean(tot_acc)
    np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})
    print('Average accuracy is : {}'.format(mean_acc))
    print()
    print('Normalized confusion matrix')
    print(cm_normalized)
    if plot:
        if cats is None:
            cats = [str(i+1) for i in range(len(cm_normalized))]
        plot_confusion_matrix(cm_normalized, labels=cats,
                              title='Normalized confusion matrix',
                              savefig=savefig, figname=figname,
                              **kwargs)

    return cm_normalized

def test_training_data(plot=False, n=None, seed=None):
    """Test the current training data set to see its performance.
    """
    clf = get_classifier(seed=seed)
    features, labels, df = get_training_data(merge_cats=True, subsample=True)
    categories = df.USD.cat.categories.tolist()
    cm = _average_class_accuracy(clf, features, labels, n=n, seed=seed,
                                 plot=plot, cats=categories)
    feature_cols = ['PSMs', 'IDSet', 'IDGroup', 'Peptide\nCount_S',
                    'Peptide\nCount_u2g', 'iBAQ\ndstrAdj']  # modify here to change features
    _feature_importance(clf, feature_cols)
    return cm
