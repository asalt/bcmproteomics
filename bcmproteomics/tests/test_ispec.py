"""Pytest suite for bcmproteomics ispec module."""

import pandas as pd

from bcmproteomics import ispec


class FakeConnection:
    """Minimal connection stub returned by patched ``filedb_connect``."""

    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


def test_get_funcats_uses_stubbed_connection(monkeypatch):
    """Ensure ``get_funcats`` relies on the patched connection and stays offline."""

    fake_conn = FakeConnection()
    captured = {}

    def fake_filedb_connect(params=None):
        captured["called"] = True
        return fake_conn

    def fake_read_sql(sql, conn, index_col=None):
        captured["sql"] = sql
        captured["conn"] = conn
        captured["index_col"] = index_col
        return pd.DataFrame(
            {
                "gene_u2gPeptiBAQAveCount": [1.0],
                "gene_GeneSymbol": ["ABC"],
                "gene_GeneDescription": ["Example"],
                "gene_FunCats": ["Category"],
            },
            index=pd.Index([123], name="gene_GeneID"),
        )

    monkeypatch.setattr(ispec, "filedb_connect", fake_filedb_connect)
    monkeypatch.setattr(ispec.pd, "read_sql", fake_read_sql)

    result = ispec.get_funcats([123])

    assert captured["called"] is True
    assert captured["conn"] is fake_conn
    assert captured["index_col"] == "gene_GeneID"
    assert captured["sql"].startswith("Select gene_GeneID")
    assert "123" in captured["sql"]
    assert "GeneSymbol" in result.columns
    assert result.index.name == "GeneID"
    assert fake_conn.closed is True
