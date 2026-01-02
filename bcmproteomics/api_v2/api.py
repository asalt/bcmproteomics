import json
import os
import re
import threading
import time
from datetime import datetime, timezone
from functools import wraps

from click import get_app_dir
from flask import Blueprint, Response, current_app, request

from bcmproteomics import ispec

bp = Blueprint("bcmproteomics_api_v2", __name__)


def _json_response(payload, status=200):
    return Response(
        json.dumps(payload, default=str),
        status=status,
        mimetype="application/json",
    )


_SCHEMA_CACHE_LOCK = threading.Lock()
_SCHEMA_CACHE = {}


def _schema_cache_ttl_seconds():
    try:
        ttl = int(float(current_app.config.get("BCMPROTEOMICS_API_V2_SCHEMA_CACHE_TTL_SECONDS", 300)))
    except (TypeError, ValueError):
        ttl = 300
    return max(0, ttl)


def _schema_cache_max_entries():
    try:
        max_entries = int(
            float(current_app.config.get("BCMPROTEOMICS_API_V2_SCHEMA_CACHE_MAX_ENTRIES", 256))
        )
    except (TypeError, ValueError):
        max_entries = 256
    return max(0, max_entries)


def _schema_cache_prune_locked(now):
    expired = [key for key, (expires_at, _) in _SCHEMA_CACHE.items() if expires_at <= now]
    for key in expired:
        _SCHEMA_CACHE.pop(key, None)

    max_entries = _schema_cache_max_entries()
    if max_entries <= 0:
        _SCHEMA_CACHE.clear()
        return

    while len(_SCHEMA_CACHE) > max_entries:
        _SCHEMA_CACHE.pop(next(iter(_SCHEMA_CACHE)), None)


def _schema_cache_get(key):
    ttl = _schema_cache_ttl_seconds()
    if ttl <= 0:
        return None

    now = time.monotonic()
    with _SCHEMA_CACHE_LOCK:
        entry = _SCHEMA_CACHE.get(key)
        if not entry:
            return None
        expires_at, payload = entry
        if expires_at <= now:
            _SCHEMA_CACHE.pop(key, None)
            return None
        return payload


def _schema_cache_set(key, payload):
    ttl = _schema_cache_ttl_seconds()
    if ttl <= 0:
        return

    now = time.monotonic()
    expires_at = now + ttl
    with _SCHEMA_CACHE_LOCK:
        _schema_cache_prune_locked(now)
        if _schema_cache_max_entries() <= 0:
            return
        _SCHEMA_CACHE[key] = (expires_at, payload)
        _schema_cache_prune_locked(now)


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
    yield "\t".join([str(c) for c in df.columns]) + "\n"
    ix = 0
    total_rows = len(df)
    while ix < total_rows:
        stop = ix + chunksize
        chunk = df.iloc[ix:stop]
        yield chunk.to_csv(header=False, sep="\t", index=False)
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


def _datadir():
    return current_app.config.get("BCMPROTEOMICS_DATADIR") or get_app_dir(
        "local-ispec", roaming=False, force_posix=True
    )


def check_auth(username, password):
    current_app.logger.info("%s is trying to login.", username)
    test_params = dict(
        user=username,
        pw=password,
        database=ispec.login_params.database,
        url=ispec.login_params.url,
    )

    conn = ispec.filedb_connect(params=test_params)
    if isinstance(conn, str):
        current_app.logger.info("%s failed with %s.", username, conn)
        return False
    try:
        conn.close()
    except Exception:
        pass
    current_app.logger.info("%s successfully logged in.", username)
    return True


def authenticate():
    """Sends a 401 response that enables basic auth."""
    return Response(
        "Could not verify your access level for that URL.\n"
        "You have to login with proper credentials",
        401,
        {"WWW-Authenticate": 'Basic realm="Login Required"'},
    )


def requires_auth(func):
    @wraps(func)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return func(*args, **kwargs)

    return decorated


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


_SAFE_IDENT_RE = re.compile(r"^[A-Za-z0-9_]+$")


def _safe_ident(identifier):
    ident = str(identifier)
    if not _SAFE_IDENT_RE.match(ident):
        raise ValueError(f"Unsafe identifier: {ident!r}")
    return ident


def _get_table_description(cursor, tablename):
    try:
        query = f"SELECT * FROM {_quote_ident(tablename)} WHERE 1=0"
        cur = cursor.execute(query)
    except Exception:
        query = f"SELECT * FROM {tablename} WHERE 1=0"
        cur = cursor.execute(query)
    return cur.description or []


def _get_table_columns(cursor, tablename):
    desc = _get_table_description(cursor, tablename)
    return [col[0] for col in desc if col and col[0] is not None]


def _resolve_fields(available, requested):
    mapping = {str(name).strip().lower(): name for name in available if name is not None}
    resolved = []
    for field in requested:
        key = str(field).strip().lower()
        if not key:
            continue
        actual = mapping.get(key)
        if actual is None:
            raise KeyError(field)
        resolved.append(actual)
    return resolved


def _resolve_table_name(cursor, requested_name):
    if not requested_name:
        return None
    try:
        tables = [row[2] for row in cursor.tables()]
    except Exception:
        return None
    mapping = {str(name).strip().lower(): name for name in tables if name is not None}
    return mapping.get(str(requested_name).strip().lower())


def _maybe_parse_datetime(value):
    if value is None:
        return None
    if not isinstance(value, str):
        return value
    text = value.strip()
    if not text:
        return value
    text = text.replace("Z", "+00:00") if text.endswith("Z") else text
    try:
        dt = datetime.fromisoformat(text)
    except ValueError:
        return value
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
    return dt

LEGACY_TABLE_DEFAULTS = {
    # Canonical defaults for key fields when FileMaker/ODBC doesn't report PKs reliably.
    # Keys are case-insensitive table names.
    "ispec_projects": dict(pk_field="prj_PRJRecNo", modified_field="prj_ModificationTS"),
    # ExperimentRuns are uniquely identified by (EXPRecNo, RunNo, SearchNo) in practice.
    "ispec_experimentruns": dict(
        pk_fields=["exprun_EXPRecNo", "exprun_EXPRunNo", "exprun_EXPSearchNo"],
        modified_field="exprun_ModificationTS",
    ),
}


def _table_defaults(table_name):
    if not table_name:
        return {}
    return LEGACY_TABLE_DEFAULTS.get(str(table_name).strip().lower(), {})


def _maybe_parse_number(value, default=0):
    if value is None:
        return default
    if isinstance(value, (int, float)):
        return int(value) if float(value).is_integer() else value
    text = str(value).strip()
    if not text:
        return default
    try:
        num = float(text)
    except ValueError:
        return text
    return int(num) if num.is_integer() else num


def _infer_pk_field(cursor, table_name, available_fields):
    defaults = _table_defaults(table_name)
    pk_default = defaults.get("pk_field")
    if pk_default:
        try:
            return _resolve_fields(available_fields, [pk_default])[0]
        except KeyError:
            return None
    try:
        rows = list(cursor.primaryKeys(table=table_name))
    except Exception:
        rows = []

    candidates = []
    for row in rows:
        name = None
        try:
            name = row.column_name
        except Exception:
            try:
                # best-effort fallback for tuple-like rows
                name = row[3]
            except Exception:
                name = None
        if name:
            candidates.append(name)

    if candidates:
        normalized = {str(name).strip().lower() for name in candidates}
        if len(normalized) == 1:
            return _resolve_fields(available_fields, [candidates[0]])[0]
        if len(candidates) == 1:
            return _resolve_fields(available_fields, [candidates[0]])[0]

    recno_candidates = [
        name for name in available_fields if str(name).strip().lower().endswith("recno")
    ]
    if len(recno_candidates) == 1:
        return recno_candidates[0]
    return None


def _infer_pk_fields(cursor, table_name, available_fields):
    defaults = _table_defaults(table_name)
    pk_fields_default = defaults.get("pk_fields")
    if pk_fields_default:
        try:
            return _resolve_fields(available_fields, pk_fields_default)
        except KeyError:
            return None

    pk_field = defaults.get("pk_field")
    if pk_field:
        try:
            return [_resolve_fields(available_fields, [pk_field])[0]]
        except KeyError:
            return None

    try:
        rows = list(cursor.primaryKeys(table=table_name))
    except Exception:
        rows = []

    candidates = []
    for row in rows:
        name = None
        try:
            name = row.column_name
        except Exception:
            try:
                name = row[3]
            except Exception:
                name = None
        if name and name not in candidates:
            candidates.append(name)

    if candidates:
        try:
            resolved = _resolve_fields(available_fields, candidates)
        except KeyError:
            resolved = None
        if resolved:
            return resolved

    recno_candidates = [
        name for name in available_fields if str(name).strip().lower().endswith("recno")
    ]
    if len(recno_candidates) == 1:
        return [recno_candidates[0]]

    return None

@bp.route("/status")
@requires_auth
def v2_status():
    return _json_response(dict(ok=True, service="bcmproteomics-legacy-api"))


@bp.route("/schema/tables")
@requires_auth
def v2_schema_tables():
    table_type = (request.args.get("table_type") or "TABLE").strip()
    include_type = _truthy(request.args.get("include_type"))

    auth = request.authorization
    username = auth.username if auth else None
    cache_key = (
        "schema_tables",
        username,
        ispec.login_params.url,
        ispec.login_params.database,
        table_type,
        include_type,
    )
    cached = _schema_cache_get(cache_key)
    if cached is not None:
        return _json_response(dict(ok=True, table_type=table_type, tables=cached))

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

        _schema_cache_set(cache_key, items)
        return _json_response(dict(ok=True, table_type=table_type, tables=items))
    finally:
        try:
            conn.close()
        except Exception:
            pass


@bp.route("/schema/tables/<tablename>/fields")
@requires_auth
def v2_schema_table_fields(tablename):
    include_meta = _truthy(request.args.get("include_meta"))

    auth = request.authorization
    username = auth.username if auth else None
    table_key = str(tablename).strip().lower()
    cache_key = (
        "schema_table_fields",
        username,
        ispec.login_params.url,
        ispec.login_params.database,
        table_key,
        include_meta,
    )
    cached = _schema_cache_get(cache_key)
    if cached is not None:
        return _json_response(cached)

    conn, error = _odbc_connect_from_request()
    if error:
        return _json_response(dict(ok=False, error=error), status=502)

    try:
        cursor = conn.cursor()
        resolved = _resolve_table_name(cursor, tablename) or tablename

        try:
            desc = _get_table_description(cursor, resolved)
        except Exception as exc:
            return _json_response(dict(ok=False, error=str(exc), table=tablename), status=400)

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

        _schema_cache_set(cache_key, payload)
        return _json_response(payload)
    finally:
        try:
            conn.close()
        except Exception:
            pass


@bp.route("/legacy/tables/<tablename>/rows")
@requires_auth
def v2_legacy_table_rows(tablename):
    fmt = (request.args.get("format") or "json").strip().lower()
    fields_raw = request.args.get("fields")  # comma-separated

    since_raw = request.args.get("since") or request.args.get("modified_since")
    since_pk_raw = request.args.get("since_pk")

    modified_field_raw = request.args.get("modified_field") or request.args.get("modified_ts")
    pk_fields_raw = request.args.get("pk_fields")
    pk_field_raw = request.args.get("pk_field") or request.args.get("pk")
    ids_raw = request.args.get("ids") or request.args.get("id")

    try:
        limit = int(float(request.args.get("limit", 1000)))
    except (TypeError, ValueError):
        return _json_response(dict(ok=False, error="Invalid limit"), status=400)
    if limit < 1 or limit > 20000:
        return _json_response(dict(ok=False, error="limit must be between 1 and 20000"), status=400)

    try:
        offset = int(float(request.args.get("offset", 0)))
    except (TypeError, ValueError):
        return _json_response(dict(ok=False, error="Invalid offset"), status=400)
    if offset < 0:
        return _json_response(dict(ok=False, error="offset must be >= 0"), status=400)

    if fmt not in {"json", "tsv"}:
        return _json_response(dict(ok=False, error=f"Unsupported format: {fmt}"), status=400)

    conn, error = _odbc_connect_from_request()
    if error:
        return _json_response(dict(ok=False, error=error), status=502)

    try:
        cursor = conn.cursor()
        resolved_table = _resolve_table_name(cursor, tablename) or tablename
        try:
            available_fields = _get_table_columns(cursor, resolved_table)
        except Exception as exc:
            return _json_response(dict(ok=False, error=str(exc), table=tablename), status=400)

        if fields_raw:
            requested_fields = [part.strip() for part in fields_raw.split(",") if part.strip()]
            if requested_fields == ["*"]:
                selected_fields = available_fields
            else:
                try:
                    selected_fields = _resolve_fields(available_fields, requested_fields)
                except KeyError as exc:
                    return _json_response(
                        dict(ok=False, error=f"Unknown field: {exc.args[0]}", table=resolved_table),
                        status=400,
                    )
        else:
            selected_fields = available_fields

        defaults = _table_defaults(resolved_table)
        modified_field = None
        pk_fields = None
        pk_field = None
        where_quoted = None
        where_parts = []
        params = []
        since_pk = None

        since = _maybe_parse_datetime(since_raw)
        if since_raw:
            modified_default = defaults.get("modified_field")
            modified_field_key = (modified_field_raw or modified_default or "").strip()
            if modified_field_key:
                try:
                    modified_field = _resolve_fields(available_fields, [modified_field_key])[0]
                except KeyError:
                    return _json_response(
                        dict(
                            ok=False,
                            error=f"Unknown modified_field: {modified_field_key}",
                            table=resolved_table,
                        ),
                        status=400,
                    )
            else:
                candidates = [
                    name
                    for name in available_fields
                    if str(name).strip().lower().endswith("modificationts")
                ]
                if len(candidates) == 1:
                    modified_field = candidates[0]
                elif not candidates:
                    return _json_response(
                        dict(
                            ok=False,
                            error="No ModificationTS-like field found; pass modified_field explicitly",
                            table=resolved_table,
                        ),
                        status=400,
                    )
                else:
                    return _json_response(
                        dict(
                            ok=False,
                            error="Multiple ModificationTS-like fields found; pass modified_field explicitly",
                            table=resolved_table,
                            candidates=candidates,
                        ),
                        status=400,
                    )

            if pk_fields_raw:
                requested = [part.strip() for part in str(pk_fields_raw).split(",") if part.strip()]
                if not requested:
                    return _json_response(dict(ok=False, error="No pk_fields provided"), status=400)
                try:
                    pk_fields = _resolve_fields(available_fields, requested)
                except KeyError as exc:
                    return _json_response(
                        dict(ok=False, error=f"Unknown pk_fields entry: {exc.args[0]}", table=resolved_table),
                        status=400,
                    )
            elif pk_field_raw:
                try:
                    pk_field = _resolve_fields(available_fields, [pk_field_raw])[0]
                except KeyError:
                    return _json_response(
                        dict(ok=False, error=f"Unknown pk_field: {pk_field_raw}", table=resolved_table),
                        status=400,
                    )
            else:
                inferred = _infer_pk_fields(cursor, resolved_table, available_fields)
                if not inferred:
                    return _json_response(
                        dict(
                            ok=False,
                            error="Unable to infer pk_field(s); pass pk_field or pk_fields explicitly",
                            table=resolved_table,
                        ),
                        status=400,
                    )
                if len(inferred) == 1:
                    pk_field = inferred[0]
                else:
                    pk_fields = inferred

            if pk_fields:
                raw = str(since_pk_raw or "").strip()
                if not raw:
                    since_pk = [0] * len(pk_fields)
                else:
                    parts = [part.strip() for part in raw.split(",")]
                    if len(parts) != len(pk_fields):
                        return _json_response(
                            dict(
                                ok=False,
                                error="since_pk must match pk_fields length",
                                table=resolved_table,
                                pk_fields=pk_fields,
                                since_pk=since_pk_raw,
                            ),
                            status=400,
                        )
                    since_pk = [_maybe_parse_number(part, default=0) for part in parts]

                clauses = []
                tie_params = []
                for ix, field in enumerate(pk_fields):
                    prefix = " AND ".join(f"{_quote_ident(pk_fields[jx])} = ?" for jx in range(ix))
                    gt = f"{_quote_ident(field)} > ?"
                    clause = gt if not prefix else f"{prefix} AND {gt}"
                    clauses.append(f"({clause})")
                    tie_params.extend(since_pk[:ix] + [since_pk[ix]])

                tie = "(" + " OR ".join(clauses) + ")"
                where_quoted = (
                    f"(({_quote_ident(modified_field)} > ?) OR "
                    f"({_quote_ident(modified_field)} = ? AND {tie}))"
                )
                params = [since, since] + tie_params
            else:
                since_pk = _maybe_parse_number(since_pk_raw, default=0)
                where_quoted = (
                    f"(({_quote_ident(modified_field)} > ?) OR "
                    f"({_quote_ident(modified_field)} = ? AND {_quote_ident(pk_field)} > ?))"
                )
                params = [since, since, since_pk]
        elif pk_field_raw:
            try:
                pk_field = _resolve_fields(available_fields, [pk_field_raw])[0]
            except KeyError:
                return _json_response(
                    dict(ok=False, error=f"Unknown pk_field: {pk_field_raw}", table=resolved_table),
                    status=400,
                )
        elif pk_fields_raw:
            requested = [part.strip() for part in str(pk_fields_raw).split(",") if part.strip()]
            if not requested:
                return _json_response(dict(ok=False, error="No pk_fields provided"), status=400)
            try:
                pk_fields = _resolve_fields(available_fields, requested)
            except KeyError as exc:
                return _json_response(
                    dict(ok=False, error=f"Unknown pk_fields entry: {exc.args[0]}", table=resolved_table),
                    status=400,
                )
            if len(pk_fields) == 1:
                pk_field = pk_fields[0]
                pk_fields = None

        ids = None
        if ids_raw:
            if pk_fields:
                return _json_response(
                    dict(ok=False, error="ids filter requires a single pk_field", table=resolved_table),
                    status=400,
                )
            parts = [p.strip() for p in str(ids_raw).split(",") if p.strip()]
            ids = []
            for part in parts:
                try:
                    ids.append(int(part))
                except ValueError:
                    ids.append(part)

            if not ids:
                return _json_response(dict(ok=False, error="No ids provided"), status=400)

            if pk_field is None:
                pk_field = _infer_pk_field(cursor, resolved_table, available_fields)
                if pk_field is None:
                    return _json_response(
                        dict(
                            ok=False,
                            error="Unable to infer pk_field; pass pk_field explicitly",
                            table=resolved_table,
                        ),
                        status=400,
                    )

            placeholders = ",".join(["?"] * len(ids))
            ids_clause = f"{_quote_ident(pk_field)} IN ({placeholders})"
            if where_quoted:
                where_quoted = f"({where_quoted}) AND ({ids_clause})"
            else:
                where_quoted = ids_clause
            params.extend(ids)

        # Ensure cursor fields are included when cursoring.
        if modified_field is not None and modified_field not in selected_fields:
            selected_fields = [modified_field] + selected_fields
        if pk_fields:
            missing = [name for name in pk_fields if name not in selected_fields]
            if missing:
                selected_fields = missing + selected_fields
        if pk_field is not None and pk_field not in selected_fields:
            selected_fields = [pk_field] + selected_fields

        order_by_raw = (request.args.get("order_by") or "").strip()
        order_by_parts = []
        if order_by_raw:
            parts = [p.strip() for p in order_by_raw.split(",") if p.strip()]
            for part in parts:
                direction = "ASC"
                name = part
                if part.startswith("-"):
                    direction = "DESC"
                    name = part[1:]
                try:
                    actual = _resolve_fields(available_fields, [name])[0]
                except KeyError:
                    return _json_response(
                        dict(ok=False, error=f"Unknown order_by field: {name}", table=resolved_table),
                        status=400,
                    )
                order_by_parts.append((actual, direction))
        elif modified_field is not None and pk_fields:
            order_by_parts.append((modified_field, "ASC"))
            for field in pk_fields:
                order_by_parts.append((field, "ASC"))
        elif modified_field is not None and pk_field is not None:
            order_by_parts.append((modified_field, "ASC"))
            order_by_parts.append((pk_field, "ASC"))
        elif pk_fields:
            for field in pk_fields:
                order_by_parts.append((field, "ASC"))
        elif pk_field is not None:
            order_by_parts.append((pk_field, "ASC"))

        order_by_quoted = [f"{_quote_ident(name)} {direction}" for name, direction in order_by_parts]

        select_cols = ", ".join(_quote_ident(name) for name in selected_fields)
        sql = f"SELECT {select_cols} FROM {_quote_ident(resolved_table)}"
        if where_quoted:
            sql += " WHERE " + where_quoted
        if order_by_quoted:
            sql += " ORDER BY " + ", ".join(order_by_quoted)

        fetch_limit = limit + 1 if fmt == "json" else limit

        try:
            cur = cursor.execute(sql, params)
        except Exception as exc:
            # fallback for drivers that don't accept quoted identifiers
            try:
                select_cols = ", ".join(_safe_ident(name) for name in selected_fields)
                sql = f"SELECT {select_cols} FROM {_safe_ident(resolved_table)}"
                if where_quoted:
                    where_unquoted = where_quoted
                    replacements = []
                    if modified_field is not None:
                        replacements.append(modified_field)
                    if pk_field is not None:
                        replacements.append(pk_field)
                    if pk_fields:
                        replacements.extend(pk_fields)
                    for name in replacements:
                        where_unquoted = where_unquoted.replace(_quote_ident(name), _safe_ident(name))
                    sql += " WHERE " + where_unquoted
                if order_by_parts:
                    order_by_unquoted = [f"{_safe_ident(name)} {direction}" for name, direction in order_by_parts]
                    sql += " ORDER BY " + ", ".join(order_by_unquoted)
            except ValueError as ident_exc:
                return _json_response(
                    dict(ok=False, error=str(exc), identifier_error=str(ident_exc), table=resolved_table),
                    status=400,
                )
            cur = cursor.execute(sql, params)

        skipped = 0
        while skipped < offset:
            batch = cur.fetchmany(min(1000, offset - skipped))
            if not batch:
                break
            skipped += len(batch)

        fetched = []
        while len(fetched) < fetch_limit:
            batch = cur.fetchmany(min(1000, fetch_limit - len(fetched)))
            if not batch:
                break
            fetched.extend(batch)

        if fmt == "tsv":
            header = "\t".join(str(name) for name in selected_fields)
            lines = [header]
            for row in fetched:
                lines.append("\t".join("" if row[ix] is None else str(row[ix]) for ix in range(len(selected_fields))))
            return Response("\n".join(lines) + "\n", mimetype="text/tab-separated-values")

        has_more = len(fetched) > limit
        fetched = fetched[:limit]
        items = [
            {selected_fields[ix]: row[ix] for ix in range(len(selected_fields))}
            for row in fetched
        ]

        next_since = None
        next_since_pk = None
        if items and modified_field is not None and (pk_field is not None or pk_fields):
            last = items[-1]
            next_since = last.get(modified_field)
            if pk_fields:
                next_since_pk = [last.get(name) for name in pk_fields]
            else:
                next_since_pk = last.get(pk_field)

        return _json_response(
            dict(
                ok=True,
                table=resolved_table,
                fields=selected_fields,
                items=items,
                rows=items,
                since=since_raw,
                since_pk=since_pk if since_raw else None,
                modified_field=modified_field,
                pk_field=pk_field,
                pk_fields=pk_fields,
                next_since=next_since,
                next_since_pk=next_since_pk,
                has_more=has_more,
                limit=limit,
                offset=offset,
            )
        )
    finally:
        try:
            conn.close()
        except Exception:
            pass


@bp.route("/experiments/<rec>_<run>_<search>/manifest")
@requires_auth
def v2_manifest(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"
    datadir = _datadir()

    meta_paths = _as_list(ispec._find_file(target=f"{stem}*.json", path=datadir))
    e2g_paths = _as_list(ispec._find_file(target=f"{stem}*_e2g*t*", path=datadir))
    psm_paths = _as_list(ispec._find_file(target=f"{stem}*psm*t*", path=datadir))

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


@bp.route("/experiments/<rec>_<run>_<search>/metadata")
@requires_auth
def v2_metadata(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    only_local = _truthy(request.args.get("only_local", "1"))
    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"
    datadir = _datadir()

    meta_path = ispec._find_file(target=f"{stem}*.json", path=datadir)
    if meta_path is None and not only_local:
        ispec.Experiment(recno, runno, searchno, data_dir=datadir, only_local=False)
        meta_path = ispec._find_file(target=f"{stem}*.json", path=datadir)

    if meta_path is None:
        return _json_response(
            dict(ok=False, error="metadata not found", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    payload = _read_metadata_file(meta_path)
    return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, metadata=payload))


@bp.route("/experiments/<rec>_<run>_<search>/e2g")
@requires_auth
def v2_e2g(rec=None, run=1, search=1):
    triplet = _normalize_triplet(rec, run, search)
    if triplet is None:
        return _json_response(dict(ok=False, error="Invalid rec/run/search triplet"), status=400)

    fmt = (request.args.get("format") or "tsv").strip().lower()
    only_local = _truthy(request.args.get("only_local", "1"))

    recno, runno, searchno = triplet
    stem = f"{recno}_{runno}_{searchno}"
    datadir = _datadir()

    e2g_path = ispec._find_file(target=f"{stem}*_e2g*t*", path=datadir)
    if e2g_path is None and only_local:
        return _json_response(
            dict(ok=False, error="e2g not found in local cache", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    exp = ispec.E2G(recno, runno, searchno, data_dir=datadir, only_local=only_local)
    df = exp.df

    if fmt in {"ispec2", "ispec", "ispec-full"}:
        rows = []
        for rec_ in df.to_dict(orient="records"):
            gene = str(rec_.get("GeneID", "")).strip()
            if not gene or gene.lower() == "nan":
                continue
            rows.append(
                dict(
                    gene=gene,
                    geneidtype="GeneID",
                    label="0",
                    iBAQ_dstrAdj=_json_safe(rec_.get("iBAQ_dstrAdj")),
                    peptideprint=_json_safe(rec_.get("PeptidePrint")),
                )
            )
        return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, e2g=rows))

    if fmt == "json":
        return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, e2g=df.to_dict(orient="records")))

    if fmt != "tsv":
        return _json_response(dict(ok=False, error=f"Unsupported format: {fmt}"), status=400)

    return Response(_df_to_tsv_stream(df), mimetype="text/tab-separated-values")


@bp.route("/experiments/<rec>_<run>_<search>/psms")
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
    datadir = _datadir()

    psm_path = ispec._find_file(target=f"{stem}*psm*t*", path=datadir)
    if psm_path is None and only_local:
        return _json_response(
            dict(ok=False, error="psms not found in local cache", recno=recno, runno=runno, searchno=searchno),
            status=404,
        )

    exp = ispec.PSMs(
        recno,
        runno,
        searchno,
        data_dir=datadir,
        presplit=presplit,
        only_local=only_local,
    )
    df = exp.df

    if fmt == "json":
        return _json_response(dict(ok=True, recno=recno, runno=runno, searchno=searchno, psms=df.to_dict(orient="records")))

    if fmt != "tsv":
        return _json_response(dict(ok=False, error=f"Unsupported format: {fmt}"), status=400)

    return Response(_df_to_tsv_stream(df), mimetype="text/tab-separated-values")
