import json
import os
import re
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


def _infer_pk_field(cursor, table_name, available_fields):
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


@bp.route("/status")
@requires_auth
def v2_status():
    return _json_response(dict(ok=True, service="bcmproteomics-legacy-api"))


@bp.route("/schema/tables")
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


@bp.route("/schema/tables/<tablename>/fields")
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
            fields = _get_table_columns(cursor, resolved)
        except Exception as exc:
            return _json_response(dict(ok=False, error=str(exc), table=tablename), status=400)

        payload = dict(ok=True, table=resolved, fields=fields)
        if include_meta:
            desc = _get_table_description(cursor, resolved)
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


@bp.route("/legacy/tables/<tablename>/rows")
@requires_auth
def v2_legacy_table_rows(tablename):
    fmt = (request.args.get("format") or "json").strip().lower()
    fields_raw = request.args.get("fields")  # comma-separated

    since_raw = request.args.get("since") or request.args.get("modified_since")
    since_pk_raw = request.args.get("since_pk")

    modified_field_raw = request.args.get("modified_field") or request.args.get("modified_ts")
    pk_field_raw = request.args.get("pk_field") or request.args.get("pk")

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

        modified_field = None
        pk_field = None
        where_quoted = None
        where_parts = []
        params = []

        since = _maybe_parse_datetime(since_raw)
        if since_raw:
            if modified_field_raw:
                try:
                    modified_field = _resolve_fields(available_fields, [modified_field_raw])[0]
                except KeyError:
                    return _json_response(
                        dict(
                            ok=False,
                            error=f"Unknown modified_field: {modified_field_raw}",
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

            if pk_field_raw:
                try:
                    pk_field = _resolve_fields(available_fields, [pk_field_raw])[0]
                except KeyError:
                    return _json_response(
                        dict(ok=False, error=f"Unknown pk_field: {pk_field_raw}", table=resolved_table),
                        status=400,
                    )
            else:
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

            try:
                since_pk = int(float(since_pk_raw or 0))
            except (TypeError, ValueError):
                return _json_response(dict(ok=False, error="Invalid since_pk"), status=400)

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

        # Ensure cursor fields are included when cursoring.
        if modified_field is not None and modified_field not in selected_fields:
            selected_fields = [modified_field] + selected_fields
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
        elif modified_field is not None and pk_field is not None:
            order_by_parts.append((modified_field, "ASC"))
            order_by_parts.append((pk_field, "ASC"))
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
                    for name in {modified_field, pk_field} - {None}:
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
        if items and modified_field is not None and pk_field is not None:
            last = items[-1]
            next_since = last.get(modified_field)
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
