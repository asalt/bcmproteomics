import base64
import json
import tempfile
import unittest
from datetime import datetime
from types import SimpleNamespace

from flask import Flask

from bcmproteomics.api_v2 import api_v2_blueprint
from bcmproteomics.api_v2 import api as api_v2_module


def _auth_headers(username="user", password="pass"):
    token = base64.b64encode(f"{username}:{password}".encode("utf-8")).decode("ascii")
    return {"Authorization": f"Basic {token}"}


def _col(name):
    # match pyodbc cursor.description (7-tuple)
    return (name, "TEXT", None, None, None, None, True)


class FakeCursor:
    def __init__(self, *, tables=None, table_descriptions=None, rows=None, pk_field=None, fail_quoted=False):
        self._tables = list(tables or [])
        self._table_descriptions = dict(table_descriptions or {})
        self._rows = list(rows or [])
        self._pk_field = pk_field
        self._fail_quoted = fail_quoted

        self.executed = []
        self.description = []
        self._pos = 0

    def tables(self, tableType=None):
        return self._tables

    def primaryKeys(self, table=None):
        if self._pk_field is None:
            return []
        return [SimpleNamespace(column_name=self._pk_field)]

    def execute(self, query, params=None):
        if self._fail_quoted and '"' in query:
            raise Exception("quoted identifiers not supported")

        self.executed.append((query, list(params or [])))
        self._pos = 0

        # Populate description for `SELECT * FROM <table> WHERE 1=0` lookups.
        if "WHERE 1=0" in query.upper():
            table = None
            parts = query.split()
            for ix, part in enumerate(parts):
                if part.upper() == "FROM" and ix + 1 < len(parts):
                    table = parts[ix + 1].strip('"')
                    break
            if table is not None:
                self.description = self._table_descriptions.get(table, [])

        return self

    def fetchmany(self, size):
        if self._pos >= len(self._rows):
            return []
        batch = self._rows[self._pos : self._pos + size]
        self._pos += len(batch)
        return batch


class FakeConn:
    def __init__(self, cursor):
        self._cursor = cursor

    def cursor(self):
        return self._cursor

    def close(self):
        return None


class TestAPIV2(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

        self.cursor = FakeCursor(
            tables=[(None, None, "iSPEC_Projects", "TABLE", None)],
            table_descriptions={
                "iSPEC_Projects": [_col("prj_PRJRecNo"), _col("prj_ModificationTS"), _col("prj_ProjectTitle")]
            },
            pk_field="prj_PRJRecNo",
        )
        self.conn = FakeConn(self.cursor)

        self._orig_login_params = api_v2_module.ispec.login_params
        self._orig_filedb_connect = api_v2_module.ispec.filedb_connect
        api_v2_module.ispec.login_params = SimpleNamespace(database="iSPEC_BCM", url="127.0.0.1")
        api_v2_module.ispec.filedb_connect = lambda params=None: self.conn

        app = Flask(__name__)
        app.config["TESTING"] = True
        app.config["BCMPROTEOMICS_DATADIR"] = self.tmp.name
        app.register_blueprint(api_v2_blueprint, url_prefix="/api/v2")
        self.client = app.test_client()

    def tearDown(self):
        api_v2_module.ispec.login_params = self._orig_login_params
        api_v2_module.ispec.filedb_connect = self._orig_filedb_connect
        self.tmp.cleanup()

    def test_status_requires_auth(self):
        resp = self.client.get("/api/v2/status")
        self.assertEqual(resp.status_code, 401)

    def test_status_ok_with_auth(self):
        resp = self.client.get("/api/v2/status", headers=_auth_headers())
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertTrue(payload["ok"])

    def test_schema_tables_lists_tables(self):
        resp = self.client.get("/api/v2/schema/tables", headers=_auth_headers())
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertIn("iSPEC_Projects", payload["tables"])

    def test_schema_table_fields_lists_fields(self):
        resp = self.client.get(
            "/api/v2/schema/tables/iSPEC_Projects/fields",
            headers=_auth_headers(),
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertEqual(payload["table"], "iSPEC_Projects")
        self.assertIn("prj_PRJRecNo", payload["fields"])
        self.assertIn("prj_ModificationTS", payload["fields"])

    def test_manifest_lists_cached_artifacts(self):
        stem = "123_1_1"
        meta_path = f"{self.tmp.name}/{stem}.json"
        e2g_path = f"{self.tmp.name}/{stem}_e2g.tsv"
        psm_path = f"{self.tmp.name}/{stem}_psms.tsv"

        with open(meta_path, "w") as handle:
            json.dump(json.dumps({"sample": "x"}), handle)
        open(e2g_path, "w").close()
        open(psm_path, "w").close()

        resp = self.client.get(f"/api/v2/experiments/{stem}/manifest", headers=_auth_headers())
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        meta_names = [item["name"] for item in payload["artifacts"]["metadata"]]
        self.assertIn(f"{stem}.json", meta_names)

    def test_metadata_reads_cached_json_string_payload(self):
        stem = "123_1_1"
        meta_path = f"{self.tmp.name}/{stem}.json"
        with open(meta_path, "w") as handle:
            json.dump(json.dumps({"sample": "x"}), handle)

        resp = self.client.get(f"/api/v2/experiments/{stem}/metadata", headers=_auth_headers())
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertEqual(payload["metadata"]["sample"], "x")

    def test_legacy_rows_builds_since_pk_query_and_paginates(self):
        self.cursor.executed.clear()
        self.cursor._rows = [
            (11, datetime(2025, 1, 1, 0, 0, 0), "A"),
            (12, datetime(2025, 1, 1, 0, 0, 0), "B"),
            (13, datetime(2025, 1, 2, 0, 0, 0), "C"),
        ]

        resp = self.client.get(
            "/api/v2/legacy/tables/iSPEC_Projects/rows"
            "?fields=prj_ProjectTitle"
            "&since=2025-01-01T00:00:00Z"
            "&since_pk=10"
            "&modified_field=prj_ModificationTS"
            "&pk_field=prj_PRJRecNo"
            "&limit=2",
            headers=_auth_headers(),
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()

        self.assertTrue(payload["has_more"])
        self.assertEqual(len(payload["items"]), 2)
        self.assertEqual(payload["items"][0]["prj_ProjectTitle"], "A")
        self.assertEqual(payload["next_since_pk"], 12)

        # Validate query shape and parameter binding for since/since_pk cursoring.
        executed_sql, executed_params = self.cursor.executed[-1]
        self.assertIn("WHERE", executed_sql.upper())
        self.assertIn("prj_ModificationTS", executed_sql)
        self.assertIn("prj_PRJRecNo", executed_sql)
        self.assertEqual(len(executed_params), 3)

    def test_legacy_rows_filters_by_ids(self):
        self.cursor.executed.clear()
        self.cursor._rows = [
            (12, "B"),
        ]

        resp = self.client.get(
            "/api/v2/legacy/tables/iSPEC_Projects/rows"
            "?fields=prj_ProjectTitle"
            "&ids=12"
            "&limit=10",
            headers=_auth_headers(),
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertEqual(payload["items"][0]["prj_ProjectTitle"], "B")

        executed_sql, executed_params = self.cursor.executed[-1]
        self.assertIn("IN", executed_sql.upper())
        self.assertIn("prj_PRJRecNo", executed_sql)
        self.assertIn(12, executed_params)


if __name__ == "__main__":
    unittest.main()
