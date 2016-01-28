"""Test for bcmproteomics
"""
import sys
import unittest
import pyodbc
sys.path.append('..')
import ispec

class TestiSPEC(unittest.TestCase):
    def test_connect_to_db(self):
        conn = ispec.filedb_connect()
        assert isinstance(conn, pyodbc.Connection)

