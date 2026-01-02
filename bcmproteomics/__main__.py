
import os
import json
from datetime import datetime, timezone
from functools import wraps
import logging
from logging.handlers import TimedRotatingFileHandler, RotatingFileHandler

from flask import (Flask, request, Response, render_template,
                   redirect, url_for)
from flask_login import (login_user, logout_user, current_user,
                         LoginManager, UserMixin, login_required)
#from flask_cache import Cache
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
server = {'bcmproteomics': '10.16.1.44',
          'jun lab': '10.13.14.171',
}

app = Flask('bcmproteomics')
app.config['CACHE_TYPE'] = 'simple'
#app.cache = Cache(app)
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
    test_params = dict(user=username, pw=password, database=ispec.login_params.database,
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


def _json_response(payload, status=200):
    return app.response_class(
        response=json.dumps(payload, default=str),
        status=status,
        mimetype='application/json',
    )


def _truthy(value):
    if value is None:
        return False
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def _normalize_triplet(rec, run, search):
    try:
        recno = int(float(rec))
        runno = int(float(run))
        searchno = int(float(search))
    except (TypeError, ValueError):
        return None
    return recno, runno, searchno


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _file_info(path):
    try:
        st = os.stat(path)
    except OSError:
        return None
    return dict(
        name=os.path.basename(path),
        size_bytes=st.st_size,
        mtime=datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
    )


def _read_metadata_file(path):
    try:
        with open(path, "r") as handle:
            payload = json.load(handle)
    except OSError:
        return None

    try:
        return json.loads(payload)
    except (TypeError, json.JSONDecodeError):
        return payload


def _df_to_tsv_stream(df, chunksize=15000):
    yield '\t'.join([str(c) for c in df.columns]) + '\n'
    ix = 0
    total_rows = len(df)
    while ix < total_rows:
        stop = ix + chunksize
        chunk = df.iloc[ix:stop]
        yield chunk.to_csv(header=False, sep='\t', index=False)
        ix = stop


def _json_safe(value):
    if value is None:
        return None
    try:
        # NaN is the only value not equal to itself
        if value != value:
            return None
    except Exception:
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def _request_db_params():
    auth = request.authorization
    if not auth:
        return None
    return dict(
        user=auth.username,
        pw=auth.password,
        database=ispec.login_params.database,
        url=ispec.login_params.url,
    )


def _odbc_connect_from_request():
    params = _request_db_params()
    if not params:
        return None, "Missing request authorization"
    conn = ispec.filedb_connect(params=params)
    if isinstance(conn, str):
        return None, conn
    return conn, None


def _quote_ident(identifier):
    if identifier is None:
        raise ValueError("identifier must be non-null")
    ident = str(identifier)
    return '"' + ident.replace('"', '""') + '"'


def _resolve_table_name(cursor, requested_name):
    if not requested_name:
        return None
    try:
        tables = [row[2] for row in cursor.tables()]
    except Exception:
        return None

    mapping = {str(name).strip().lower(): name for name in tables if name is not None}
    return mapping.get(str(requested_name).strip().lower())

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

#@app.cache.memoize(timeout=1000)
def get_e2g_exp(rec, run=1, search=1):
    # ispec.params = get_ispec_params()
    exp = ispec.E2G(rec, run, search, data_dir=DATADIR)
    if 'GeneID' not in exp.df.columns:
        exp.df['GeneID'] = exp.df.index
    if not exp.df.index.is_unique:
        exp.df.reset_index(drop=True, inplace=True)
    return exp

#@app.cache.memoize(timeout=1000)
def get_exp(rec, run=1, search=1):
    exp = ispec.Experiment(rec, run, search)
    return exp

#@app.cache.memoize(timeout=1000)
def get_psms(rec, run=1, search=1, presplit=False):
    if isinstance(presplit, str) and presplit.lower() == 'true' or presplit == 1:
        presplit_ = True
    else:
        presplit_ = False
    exp = ispec.PSMs(rec, run, search, data_dir=DATADIR, presplit=presplit_)
    return exp


@app.route('/api/v2/status')
@requires_auth
def v2_status():
    return _json_response(
        dict(
            ok=True,
            service="bcmproteomics-legacy-api",
        )
    )


@app.route('/api/v2/schema/tables')
@requires_auth
def v2_schema_tables():
    table_type = (request.args.get("table_type") or "TABLE").strip()
    include_type = _truthy(request.args.get("include_type"))

    conn, error = _odbc_connect_from_request()
    if error:
        return _json_response(dict(ok=False, error=error), status=502)

    try:
        cursor = conn.cursor()
        rows = cursor.tables(tableType=table_type) if table_type else cursor.tables()
        items = []
        for row in rows:
            name = row[2]
            typ = row[3] if len(row) > 3 else None
            if include_type:
                items.append(dict(name=name, type=typ))
            else:
                items.append(name)

        if include_type:
            items.sort(key=lambda item: str(item.get("name", "")).lower())
        else:
            items = sorted({str(name) for name in items if name is not None}, key=str.lower)

        return _json_response(dict(ok=True, table_type=table_type, tables=items))
    finally:
        try:
            conn.close()
        except Exception:
            pass


@app.route('/api/v2/schema/tables/<tablename>/fields')
@requires_auth
def v2_schema_table_fields(tablename):
    include_meta = _truthy(request.args.get("include_meta"))

    conn, error = _odbc_connect_from_request()
    if error:
        return _json_response(dict(ok=False, error=error), status=502)

    try:
        cursor = conn.cursor()
        resolved = _resolve_table_name(cursor, tablename) or tablename

        try:
            query = f"SELECT * FROM {_quote_ident(resolved)} WHERE 1=0"
            cur = cursor.execute(query)
        except Exception:
            try:
                query = f"SELECT * FROM {resolved} WHERE 1=0"
                cur = cursor.execute(query)
            except Exception as exc:
                return _json_response(
                    dict(ok=False, error=str(exc), table=tablename),
                    status=400,
                )

        desc = cur.description or []
        fields = [col[0] for col in desc if col and col[0] is not None]
        payload = dict(ok=True, table=resolved, fields=fields)
        if include_meta:
            payload["columns"] = [
                dict(
                    name=col[0],
                    type=str(col[1]),
                    display_size=col[2],
                    internal_size=col[3],
                    precision=col[4],
                    scale=col[5],
                    nullable=col[6],
                )
                for col in desc
            ]
        return _json_response(payload)
    finally:
        try:
            conn.close()
        except Exception:
            pass


@app.route('/api/v2/experiments/<rec>_<run>_<search>/manifest')
@requires_auth
def v2_manifest(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"

    meta_paths = _as_list(ispec._find_file(target=f"{stem}*.json", path=DATADIR))
    e2g_paths = _as_list(ispec._find_file(target=f"{stem}*_e2g*t*", path=DATADIR))
    psm_paths = _as_list(ispec._find_file(target=f"{stem}*psm*t*", path=DATADIR))

    return _json_response(
        dict(
            ok=True,
            recno=recno,
            runno=runno,
            searchno=searchno,
            artifacts=dict(
                metadata=[fi for fi in (_file_info(p) for p in meta_paths) if fi],
                e2g=[fi for fi in (_file_info(p) for p in e2g_paths) if fi],
                psms=[fi for fi in (_file_info(p) for p in psm_paths) if fi],
            ),
        )
    )


@app.route('/api/v2/experiments/<rec>_<run>_<search>/metadata')
@requires_auth
def v2_metadata(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    only_local = _truthy(request.args.get("only_local", "1"))
    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"

    meta_path = ispec._find_file(target=f"{stem}*.json", path=DATADIR)
    if meta_path is None and not only_local:
        ispec.Experiment(recno, runno, searchno, data_dir=DATADIR, only_local=False)
        meta_path = ispec._find_file(target=f"{stem}*.json", path=DATADIR)

    if meta_path is None:
        return _json_response(
            dict(ok=False, error="metadata not found", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    payload = _read_metadata_file(meta_path)
    return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, metadata=payload))


@app.route('/api/v2/experiments/<rec>_<run>_<search>/e2g')
@requires_auth
def v2_e2g(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    fmt = (request.args.get("format") or "tsv").strip().lower()
    only_local = _truthy(request.args.get("only_local", "1"))

    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"
    e2g_path = ispec._find_file(target=f"{stem}*_e2g*t*", path=DATADIR)
    if e2g_path is None and only_local:
        return _json_response(
            dict(ok=False, error="e2g not found in local cache", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    exp = ispec.E2G(recno, runno, searchno, data_dir=DATADIR, only_local=only_local)
    df = exp.df

    if fmt in {"ispec2", "ispec", "ispec-full"}:
        rows = []
        for rec in df.to_dict(orient="records"):
            gene = str(rec.get("GeneID", "")).strip()
            if not gene or gene.lower() == "nan":
                continue
            rows.append(
                dict(
                    gene=gene,
                    geneidtype="GeneID",
                    label="0",
                    iBAQ_dstrAdj=_json_safe(rec.get("iBAQ_dstrAdj")),
                    peptideprint=_json_safe(rec.get("PeptidePrint")),
                )
            )
        return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, e2g=rows))

    if fmt == "json":
        return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, e2g=df.to_dict(orient="records")))

    if fmt != "tsv":
        return _json_response(dict(ok=False, error=f"Unsupported format: {fmt}"), status=400)

    return Response(_df_to_tsv_stream(df), mimetype='text/tab-separated-values')


@app.route('/api/v2/experiments/<rec>_<run>_<search>/psms')
@requires_auth
def v2_psms(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    fmt = (request.args.get("format") or "tsv").strip().lower()
    only_local = _truthy(request.args.get("only_local", "1"))
    presplit = _truthy(request.args.get("presplit"))

    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"
    psm_path = ispec._find_file(target=f"{stem}*psm*t*", path=DATADIR)
    if psm_path is None and only_local:
        return _json_response(
            dict(ok=False, error="psms not found in local cache", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    exp = ispec.PSMs(
        recno,
        runno,
        searchno,
        data_dir=DATADIR,
        presplit=presplit,
        only_local=only_local,
    )
    df = exp.df

    if fmt == "json":
        return _json_response(
            dict(ok=True, recno=recno, runno=runno, searchno=searchno, psms=df.to_dict(orient="records"))
        )

    if fmt != "tsv":
        return _json_response(dict(ok=False, error=f"Unsupported format: {fmt}"), status=400)

    return Response(_df_to_tsv_stream(df), mimetype='text/tab-separated-values')

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
    _toexclude = [x for x in exp.__dict__.keys() if x.startswith('_') or x == "metadata"]
    #for attr in [x for x in exp.__dict__.keys() if not x.startswith('_')]:
    if exp.metadata is not None:
        app.logger.info("sending metadata through json export")
        d['metadata'] = exp.metadata.to_json()
    for attr in [*set(exp.__dict__.keys()) - set(_toexclude)]:
        # logging.info(attr)
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
