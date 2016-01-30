"""Functions and classes used for querying and manipulating data in
BCM Proteomics iSPEC.

"""
from __future__ import print_function
import numpy as np
import pandas as pd
from collections import OrderedDict
import pyodbc
from getpass import getpass
from bcmproteomics.e2gitems import e2gcolumns
from bcmproteomics.classify import score_experiments

user = None
pw   = None
url  = None
database = None
params = {'user': user, 'pw': pw, 'url': url, 'database': database}


class E2G:
    """ An object for working with iSPEC Proteomics data at BCM

    Attributes
    ----------
    df : pandas DataFrame of the experimental data

    recno : experimental record number

    runno : experimental run number

    sample : experimental sample cell/tissue

    treatment : any applied treatment

    exptype : type of experiment (Profile, Affinity, ...)

    genotype : genotype for the experiment
    ----------

    Initialize a new experiment object and grab your data:
    >>> import BCM_proteomics as BCM
    >>> BCM.user = 'username'
    >>> BCM.pw = 'password'
    >>> exp = BCM.E2G(recno=12345, runno=1) # note runno defaults to 1

    ----------

    Multiple E2G instances can be joined using join_exps:
    >>> from bcmproteomcs import ispec
    >>> exp1 = ispec.E2G(12345,1)
    >>> exp2 = ispec.E2G(12345,2)
    >>> exp1_2 = ispec.join_exps(exp1, exp2)

    note that many methods on E2G are (currently) unavailable with a joined E2G instance
    ----------
    """
    def __init__(self, recno=None, runno=None, searchno=None):
        self.recno = recno
        self.runno = runno
        self.searchno = searchno
        #self.techrepno =
        self.sample = ''
        self.treatment = ''
        self.exptype = ''
        self.genotype = ''
        self.added_by = ''
        if recno is not None:
            self.get_exprun(recno, runno, searchno) # self._df gets set through here
        else:
            self._df = pd.DataFrame() # else set as empty dataframe
        self._joined = False
        self.ibaq_normalize = None

    def __repr__(self):
        return 'Record number {}, run number {}'.format(self.recno, self.runno)

    @property
    def df(self):
        """Acccess (potentially in the future) dictionary of pandas DataFrames for record sections."""
        return self._df

    def get_exprun(self, recno=None, runno=1, searchno=1):
        """queries iSPEC database and grabs the gene product table which is stored in self.df
        """
        if not recno:
            raise ValueError('recno must be specified')
        if runno is None:
            runno = 1  # make sure default runno is 1
        if searchno is None:
            searchno = 1

        conn = filedb_connect()
        if isinstance(conn, str):
            return # failed to make connection, user is informed in filedb_connect()

        sql = ("SELECT {} from iSPEC_BCM.iSPEC_exp2gene where "
               "e2g_EXPRecNo={} "
               "AND e2g_EXPRunNo={} "
               "AND e2g_EXPSearchNo={}").format(', '.join(e2gcolumns),
                                                recno,
                                                runno,
                                                searchno)

        self._df = self._construct_df(sql, conn)
        sql_description = ("SELECT exp_EXPClass, exp_Extract_CellTissue, exp_Exp_Description,"
                           "exp_Extract_Treatment, exp_IDENTIFIER, exp_Extract_No, "
                           "exp_Extract_Genotype, exp_AddedBy, exp_Digest_Type, exp_Digest_Enzyme "
                           "from iSPEC_BCM.iSPEC_Experiments where exp_EXPRecNo={}").format(recno)
        info = pd.read_sql(sql_description, conn).to_dict('list') # a 1 row dataframe
        self.sample = ''.join(item for item in info.get('exp_Extract_CellTissue', '') if item)
        self.treatment = ''.join(item for item in info.get('exp_Extract_Treatment', '') if item)
        self.exptype = ''.join(item for item in info.get('exp_EXPClass', '') if item)
        self.description = [item for items in info['exp_Exp_Description'] for item in items.splitlines()]
        self.genotype = ''.join(item for item in info.get('exp_Extract_Genotype', '') if item)
        self.added_by = ''.join(item for item in info.get('exp_AddedBy', '') if item)
        self.digest_type = ''.join(item for item in info.get('exp_Digest_Type', '') if item)
        self.digest_enzyme = ''.join(item for item in info.get('exp_Digest_Enzyme', '') if item)
        #print( ''.join(int(item) for item in info.get('exp_Extract_No', 0) if str(item).isdigit()))
        self.extract_no = ''.join(str(item) for item in info.get('exp_Extract_No', 0) if item)
        self.recno = recno
        self.runno = runno
        self.searchno = searchno

        if self.extract_no:
            sql_extractdata = ("SELECT ext_Fractions, ext_Protocol FROM "
                               "iSPEC_BCM.iSPEC_Extracts "
                               "WHERE ext_ExtRecNo={}").format(self.extract_no)
            extract_info = pd.read_sql(sql_extractdata, conn).to_dict('list') # a 1 row dataframe
            self.extract_fractions = ''.join(extract_info.get('ext_Fractions', ''))
            self.extract_protocol = ''.join(extract_info.get('ext_Protocol', '')).replace('\r', ', ')
        conn.close()
        return self

    @staticmethod
    def _construct_df(sql, conn):
        """Construct a pandas dataframe from data in the iSPEC.
        """
        df = pd.read_sql(sql, conn, index_col='e2g_GeneID')
        df.rename(columns={k: k.split('e2g_')[1] for k in 
               [e2gcol for e2gcol in df.columns if e2gcol.startswith('e2g_')]},
             inplace=True)
        df.rename(columns={'n_iBAQ_dstrAdj': 'iBAQ_dstrAdj',},
                  inplace=True)

        df.index.rename('GeneID',inplace=True)
        #df.rename(columns={'e2g_GeneID':'GeneID'}, inplace=True)
        geneidlist = [str(int(x)) for x in df.index.tolist() if not np.isnan(x)]
        genesql  = "Select gene_GeneID, gene_u2gPeptiBAQAveCount, "\
                   "gene_GeneSymbol, gene_GeneDescription, "\
                   "gene_FunCats " \
                   "from iSPEC_BCM.iSPEC_Genes "\
                   "where gene_GeneID in ({})".format(', '.join(geneidlist))
        if geneidlist:
            genedf = pd.read_sql(genesql, conn, index_col='gene_GeneID')  # all headers start with gene_
            generename = {c: c.split('gene_')[1] for c in genedf.columns}
            genedf.rename(columns=generename, inplace=True)
            genedf.index.rename('GeneID', inplace=True)
            genedf.rename(columns={'u2gPeptIBAQAveCount':'GeneCapacity'}, inplace=True)
            #df['e2g_GeneCapacity'] = df['e2g_nGPArea_Sum_dstrAdj']/df['e2g_nGPArea_Sum_dstrAdj']
            #funcat_table = pd.read_table('E:/reference_databases/iSPEC_genes_Aug_2015.tab')
            #df = pd.merge(df, funcat_table, how='left', on='GeneID')
            #df = pd.merge(df, genedf, on='GeneID')
            df = df.join(genedf)
            df['FunCats'].fillna('', inplace=True)
        df['GeneID'] = df.index  # put index in its own column as well for easy access    
        return df

    def strict(self, df=None):
        """Returns a strict selection of gene products from the dataframe"""

        if self._joined:
            return df[((df['IDSet_x'] < 3) & (df['IDGroup_x'] <= 3) |
                      (df['IDSet_y'] < 3) & (df['IDGroup_y'] <= 3))]
                       
        if df is None:
            df = self.df
            
        return df[(df['IDSet'] < 3) &
                  (df['IDGroup'] <= 3)]

    @property
    def tfs_coRs(self):
        """Return gene products annotated as a DBTF or CoRegulator"""
        return self.df[self.df['FunCats'].str.contains('DBTF|TC|DBTC|CBTC')]

    @property
    def tfs(self):
        """Return gene products annotated as a DBTF"""
        return self.df[self.df['FunCats'].str.contains('DBTF')]
                       
    @property
    def kinases(self):
        """Return gene products annotated as a kinase"""
        return self.df[self.df['FunCats'].str.contains('KI')]

    def USD_selector(self, df=None):
        """Returns U, S, D DataFrames for a joined E2G instance"""
        if self._joined is False:
            raise ValueError('Not a joined E2G instance.')
        if df is None:
            df = self.df
        ret = list()
        for cat in list('USD'):
            ret.append(df[df.USD == cat])
        return ret

    def summary(self):
        """Print out a quick summary of the types of gene products observed in the experiment.
        """

        if len(self._df) == 0:
            raise ValueError('DataFrame is empty!')

        searchstr = OrderedDict({'TFs': 'DBTF',
                                 'Tfs and CoRs': 'DBTF|TC|DBTC|CBTC',
                                 'Kinases': 'KI',
                                 }
                                )
        print('\nSummary for record {}, run {}.\n'.format(self.recno, self.runno, ))

        if self._joined:    
            dfU, dfS, dfD = self.USD_selector()
            print('Total gene products found : '
                  '{} same, {} up, {} down'.format(len(dfS), len(dfU), len(dfD)))
            print('Total strict gene products found : '
                  '{} same, {} up, {} down'.format(len(self.strict(dfS)),
                                                   len(self.strict(dfU)),
                                                   len(self.strict(dfD))))
            print('-'*15,'\n')
            for s in searchstr:
                df_subsection = self.df[(self.df['FunCats'].str.contains(searchstr[s]))]
                dfU, dfS, dfD = self.USD_selector(df_subsection)
                print('Total {} gene products found : '
                      '{} same, {} up, {} down'.format(s, len(dfS), len(dfU), len(dfD)))
                print('Total {} strict gene products found : '
                      '{} same, {} up, {} down'.format(s,
                                                       len(self.strict(dfS)),
                                                       len(self.strict(dfU)),
                                                       len(self.strict(dfD))))
                print('-'*15,'\n')
            return

        print('Total gene products found : {}'\
              '\n\tiBAQ total : {:4f}.'.format(len(self._df),  self._df.iBAQ_dstrAdj.sum()))

        df_strict = self.strict()
        print('Total strict gene products found : {}'\
              '\n\tiBAQ total : {:4f}.'.format(len(df_strict),
                                               df_strict.iBAQ_dstrAdj.sum()))
        print('-'*15,'\n')
        for s in searchstr:
            df_subsection = self.df[(self.df['FunCats'].str.contains(searchstr[s]))]
            print('Total {} gene products found : {}'\
                  '\n\tiBAQ total : {:4f}'.format(s, len(df_subsection),
                                                  df_subsection.iBAQ_dstrAdj.sum()))

            df_sub_strict = self.strict(df_subsection)
            print('Total {} strict gene products found : {}'\
                  '\n\tiBAQ total : {:4f}.'.format(s, len(df_sub_strict),
                                                   df_sub_strict.iBAQ_dstrAdj.sum()))
            print('-'*15,'\n')



def _getlogin():
    """Checks if the username and password have been established for the current python session.
    If username or password is undefined, will prompt the user.

    Defaults to bcmproteomics.
    """
    servers = {'bcmproteomics': '10.16.2.74',
               'jun lab': '10.13.14.171',
               }
    databases = {'10.16.2.74': ['iSPEC_BCM', 'iSPEC_BCM_psms'],
                 '10.13.14.171': ['iSPEC'],
                 }
    if params.get('user') is None:
        print('Username is not set')
        user = input('Enter your iSPEC username : ')
        params['user'] = user
    if params.get('pw') is None:
        print('Password is not set')
        pw = getpass('Enter your password : ')
        params['pw'] = pw
    if params.get('url') is None:
        print('iSPEC is not set, the options are:')
        print(*[server for server in servers], sep='\n')
        server = input('Select an iSPEC : ').strip()
        params['url'] = servers.get(server, '10.16.2.74')
    elif 'url' in params:
        params['url'] = servers.get(url, '10.16.2.74')

    if params.get('database') is None:
        server_url = params['url']
        if len(databases.get(server_url, [])) == 1:
            params['database'] = databases[server][0]
        else:
            print('iSPEC database is not set, the options are:')
            print(*[database for database in databases[server_url]], sep='\t')
            db = input('Select an iSPEC database: ').strip()
            params['database'] = databases.get(db, 'iSPEC_BCM')

    return params

def filedb_connect():
    """Establish a connection to the BCM Proteomics FileMaker Server.

    Returns a database connection, which can then be used to execute sql statements.

    If the password has not been established yet, user will be prompted.
    Password is only stored in memory for the duration of the python session.
    """
    
    params = _getlogin()

    driver    = 'DRIVER={FileMaker ODBC};'

    login_info = 'UID={user};PWD={pw}'.format(**params)

    server_info = 'SERVER={url};'.format(**params)

    database_info = 'DATABASE={database};'.format(**params)

    conn_string = "{0}{1}{2}{3}".format(driver, server_info, database_info, login_info)
    try:
        conn = pyodbc.connect(conn_string)
    except pyodbc.Error as e:  # any ODBC error is returned to python as "Error"
        error = e.__str__()
        if 'password incorrect' in error.lower():
            print('Invalid password.')
        elif 'account' in error.lower():
            print('Incorrect account.')
        elif 'failed to connect' in error.lower():
            print('Error making connection, check the address.')
        else:
            print('Error with username or password')  # maybe raise an error instead?
        params.clear()
        return error

    return conn

def get_funcats(geneidlist):
    """Takes a list of gene ids and returns the funcats table from iSPEC in the form of a 
    pandas DataFrame.

    The index of the dataframe is the geneids (as the pandas default type of float64.
    For convienence, the geneids are also found within their own column as an object type.

    The funcats column of the DataFrame has had each null value filled with an empty string, 
    which allows for easy string searching with pandas.
    """
    if not isinstance(geneidlist, list) or len(geneidlist) == 0:
        raise TypeError('Input must be a list with at least 1 element')
    # need to make sure every element in the geneidlist is a string
    conn = filedb_connect()
    if isinstance(conn, str):
        return conn # failed to make connection, user is informed in filedb_connect()

    genesql  = "Select gene_GeneID, gene_u2gPeptiBAQAveCount, "\
               "gene_GeneSymbol, gene_GeneDescription, "\
               "gene_FunCats " \
               "from iSPEC_BCM.iSPEC_Genes "\
               "where gene_GeneID in ({})".format(', '.join([str(x) for x in geneidlist]))
                                                             

    genedf = pd.read_sql(genesql, conn, index_col='gene_GeneID')  # all headers start with gene_
    generename = {c: c.split('gene_')[1] for c in genedf.columns}
    genedf.rename(columns=generename, inplace=True)
    genedf.index.rename('GeneID', inplace=True)
    genedf.rename(columns={'u2gPeptIBAQAveCount':'GeneCapacity'}, inplace=True)
    genedf['GeneID'] = [str(int(x)) for x in genedf.index.tolist()]
    genedf['FunCats'].fillna('', inplace=True)
    conn.close()
    return genedf

def get_geneids(taxonid):
    """Get all gene ids for a given taxon"""
    if not isinstance(taxonid, int):
        try:
            taxonid = int(taxonid.strip())
        except:
            raise TypeError('Input must be a taxon id')

    conn = filedb_connect()
    if isinstance(conn, str):
        return conn # failed to make connection, user is informed in filedb_connect()
    sql = "SELECT gene_GeneID from iSPEC_BCM.iSPEC_Genes "\
           "WHERE gene_TaxonID={}".format(taxonid)
    cursor = conn.execute(sql)
    gids = cursor.fetchall()
    genes = [int(g) for gid in gids for g in gid]
    return genes

def join_exps(exp1, exp2, normalize=None, seed=None):
    """Nice outer join two experiments based on their geneids. Useful for comparing between experiments.
    Keeps gene information as well.
    
    Parameters
    ----------
    exp1, exp2 : non jonied E2G objects

    normalize  : 'mean', 'median', or None (Optional)
                  method of normalization of iBAQ_dstrAdj
                  Has the side effect of adding ibaq_norm column to each input DataFrame.

    seed       : set state for random forest classifier (Optional)

    Optional seed argument sets seed for random forest classifier.
    """
    if any(type(exp) is not E2G for exp in [exp1, exp2]):
        raise  TypeError('Incorrect input type')

    
    if normalize is not None:
        for exp in [exp1, exp2]:
            if normalize.strip() == 'mean':
                exp.df['ibaq_norm'] = exp.df.iBAQ_dstrAdj/exp.df.iBAQ_dstrAdj.mean()
            elif normalize.strip() == 'median':
                exp.df['ibaq_norm'] = exp.df.iBAQ_dstrAdj/exp.df.iBAQ_dstrAdj.median()
            exp.df.rename(columns={'iBAQ_dstrAdj':'iBAQ_dstrAdj_raw'}, copy=True)
            exp.df.rename(columns={'ibaq_norm':'iBAQ_dstrAdj'}, copy=True)
            
    funcat_cols = ['GeneCapacity', 'GeneSymbol', 'GeneDescription', 'FunCats', 'GeneID']
    funcat1 = exp1.df[funcat_cols]
    funcat2 = exp2.df[funcat_cols]
    funcats = pd.concat([funcat1, funcat2]).drop_duplicates()

    joinexp = E2G()
    for attr in [key for key in joinexp.__dict__.keys() if
                 not key.startswith('_')]:
        
        value = str(exp1.__getattribute__(attr)) + \
                 ' vs ' + \
                 str(exp2.__getattribute__(attr))
        joinexp.__setattr__(attr, value)
        
    nofuncats = [col for col in exp1.df.columns if col not in funcat_cols]
    joinexp._df = exp1.df[nofuncats].join(exp2.df[nofuncats],
                                          lsuffix = '_x', rsuffix='_y',
                                          how='outer', )
    joinexp._df = joinexp._df.join(funcats, how='left')
    joinexp._df['GeneID'] = [str(int(x)) if not np.isnan(x) else x
                             for x in joinexp.df.index] # for convienence
    joinexp._joined = True
    joinexp.ibaq_normalize = normalize
    try:
        score_experiments(joinexp, seed=seed)
    except Exception as e:
        print('Error scoring experiments')
        print(e)
    return joinexp
