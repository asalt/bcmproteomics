
import os
import json
from functools import wraps
import logging
from logging.handlers import TimedRotatingFileHandler, RotatingFileHandler

from flask import (Flask, request, Response, render_template,
                   redirect, url_for)
from flask_login import (login_user, logout_user, current_user,
                         LoginManager, UserMixin, login_required)
from flask_cache import Cache
from click import get_app_dir

from bcmproteomics import ispec



DATADIR = None
DATADIR = get_app_dir('local-ispec', roaming=False, force_posix=True)
LOGDIR = os.path.join(DATADIR, 'logs')
if not os.path.exists(LOGDIR):
    os.mkdir(LOGDIR)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler = TimedRotatingFileHandler(os.path.join(LOGDIR, 'bcmproteomics.log'), when='midnight', backupCount=60)
handler.setLevel(logging.INFO)
handler.setFormatter(formatter)


print('Data stored locally at', DATADIR)
if not os.path.exists(DATADIR):
    os.mkdir(DATADIR)
server = {'bcmproteomics': '10.16.2.74',
          'jun lab': '10.13.14.171',
}

app = Flask('bcmproteomics')
app.config['CACHE_TYPE'] = 'simple'
app.cache = Cache(app)
app.logger.setLevel(logging.INFO)
app.logger.addHandler(handler)
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'
login_manager.login_message_category = 'warning'

# def get_ispec_params():
#     # Set globably once for database access
#     params = ispec.Conf()
#     params.config['iSPEC']['user'] = 'flask_login'
#     params.config['iSPEC']['pw'] = 'flask_login'
#     params.config['iSPEC']['database'] = 'iSPEC_BCM'
#     params.config['iSPEC']['url'] = '10.16.2.74'
#     # ispec.params['website'] = '10.16.3.148:5000'
#     return params

def check_auth(username, password):
    """This function is called to check if a username /
    password combination is valid.
    """
    app.logger.info('{} is trying to login.'.format(username))
    # params.config['iSPEC']['user'] = username
    # params.config['iSPEC']['pw'] = password
    # params.config['iSPEC']['database'] = 'iSPEC_BCM'
    test_params = dict(user=username, pw=password, database=ispec.login_params.database
                       url=ispec.login_params.url
    )
    # ispec.params['website'] = '10.16.3.148:5000'

    conn = ispec.filedb_connect(params=test_params)
    if isinstance(conn, str):
        # app.logger.warning('{} is unable to register to {}.'.format(username, ispec_db))
        # ispec.params['database'] = 'iSPEC_BCM'
        # ispec.params['url'] = '10.16.2.74'
        app.logger.info('{} failed with {}.'.format(username, conn))
        return False
    app.logger.info('{} successfully logged in.'.format(username))
    return True


def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response(
    'Could not verify your access level for that URL.\n'
    'You have to login with proper credentials', 401,
    {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)
    return decorated

@app.route('/', methods=['GET', 'POST'])
def home():
    return redirect(url_for('login'))


@app.route('/login', methods=['GET', 'POST'])
@requires_auth
def login():
    return Response(
        'Successful login.',
        200,
        {'WWW-Authenticate': 'Basic realm="Login Required"'})

@app.route('/api/data/<rec>_<run>_<search>/<typeof>/<presplit>')
@requires_auth
def data(rec=None, run=1, search=1, typeof='e2g', presplit=0):
    """
    typeof : type of data to return. Options include e2g or psms"""

    if typeof == 'e2g':
        exp = get_e2g_exp(rec, run, search)
    elif typeof == 'psms':
        exp = get_psms(rec, run, search, presplit=presplit)
    # json_data = exp.df.to_json()
    # return render_template('expapi.html', json=json_data)
    def normalize(value, precision=7):
        precision_fmt = '{}f'.format(precision)
        if isinstance(value, str):
            return value
        if isinstance(value, int):
            return str(value)
        try:
            return format(value, precision_fmt)
        except TypeError:
            print('Type error with ', value)
            return str(value)
    def gen(chunksize=15000):
        # for col in exp.df.columns:
        #     yield exp.df[col].to_json()
        yield '\t'.join(exp.df.columns) + '\n'
        ix = 0
        total_rows = len(exp.df)
       	while ix < total_rows:
            stop = ix + chunksize
            chunk = exp.df.iloc[ix:stop]
            yield chunk.to_csv(header=False, sep='\t')
            ix = stop

    # app.logger.info('{!r} has {} lines'.format(exp, len(exp.df)))
    # print('{!r} has {} lines'.format(exp, len(exp.df)))
    return Response(gen(), mimetype='text/csv')

@app.cache.memoize(timeout=1000)
def get_e2g_exp(rec, run=1, search=1):
    # ispec.params = get_ispec_params()
    exp = ispec.E2G(rec, run, search, data_dir=DATADIR)
    if 'GeneID' not in exp.df.columns:
        exp.df['GeneID'] = exp.df.index
    if not exp.df.index.is_unique:
        exp.df.reset_index(drop=True, inplace=True)
    return exp

@app.cache.memoize(timeout=1000)
def get_exp(rec, run=1, search=1):
    exp = ispec.Experiment(rec, run, search)
    return exp

@app.cache.memoize(timeout=1000)
def get_psms(rec, run=1, search=1, presplit=False):
    if isinstance(presplit, str) and presplit.lower() == 'true' or presplit == 1:
        presplit_ = True
    else:
        presplit_ = False
    exp = ispec.PSMs(rec, run, search, data_dir=DATADIR, presplit=presplit_)
    return exp

@app.route('/api/funcats/<gids>')
@requires_auth
def funcats(gids):
    print(gids)
    gidlist = gids.split(',')
    funcats = ispec.get_funcats(gidlist)
    return app.response_class(response=funcats.to_json(),
                              status=200,
                              mimetype='application/json')

@app.route('/api/meta/<rec>_<run>_<search>')
@requires_auth
def meta(rec=None, run=1, search=1):
    exp = get_exp(rec, run, search)
    d = dict()
    for attr in [x for x in exp.__dict__.keys() if not x.startswith('_')]:
        d[attr] = exp.__getattribute__(attr)
    json_data = json.dumps(d)
    return app.response_class(response=json_data,
                              status=200,
                              mimetype='application/json')

@app.route('/api/geneids/<int:taxonid>')
@requires_auth
def geneids(taxonid):
    geneids = ispec.get_geneids(taxonid)
    json_data = json.dumps(geneids)
    return app.response_class(response=json_data,
                              status=200,
    mimetype='application/json')

if __name__ == '__main__':
    # app.run(host='0.0.0.0', debug=False, threaded=True)
    #app.run(host='0.0.0.0', port=6000, debug=False, threaded=True)
    from tornado.wsgi import WSGIContainer
    from tornado.httpserver import HTTPServer
    from tornado.ioloop import IOLoop

    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(6000)
    IOLoop.instance().start()
