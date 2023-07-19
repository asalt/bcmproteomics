"""Functions and classes used for querying and manipulating data in
BCM Proteomics iSPEC.

"""
from __future__ import print_function
import json
import os
import re
import configparser
from getpass import getpass
from warnings import warn
from glob import glob

import numpy as np
import pandas as pd
from collections import OrderedDict
from functools import lru_cache
import click
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

from click import get_app_dir

CONF_DIR = get_app_dir("ispec", roaming=False, force_posix=True)
CONF_FILE = os.path.join(CONF_DIR, "ispec.conf")

if not os.path.exists(CONF_DIR):
    os.makedirs(CONF_DIR)


class Conf:

    helps = dict(
        url="iSPEC IP address",
        database="iSPEC database name",
        user="iSPEC username",
        pw="password",
    )

    # parameters we get out of config['iSPEC'] configparser
    attrs = ("url", "database", "user", "pw")

    CONF_FILE = CONF_FILE

    def __init__(self):
        self._user = None
        self._pw = None
        self._database = None
        self._url = None
        self.save_all = None
        self.config = self.read_conf(self.CONF_FILE)

        # self.orig_params = { attr: getattr(self, attr) for attr in self.attrs }

    @lru_cache()
    def read_conf(self, CONF_FILE):

        config = configparser.ConfigParser()
        config.optionxform = str

        if not os.path.exists(CONF_FILE):
            config["iSPEC"] = {
                "url": "",
                "database": "iSPEC_BCM",
                "user": "",
                "pw": "",
            }
            click.echo("Created ispec conf file at {}".format(CONF_FILE))
            with open(CONF_FILE, "w") as f:
                config.write(f)

        with open(CONF_FILE, "r") as f:
            config.read_file(f)

        if "iSPEC" not in config:
            click.echo("iSPEC section missing, config file invalid.")
            return
        self.orig_params = dict(config["iSPEC"].items())
        return config

    # @lru_cache()
    def get_item(self, item):
        # ispec = self.config['iSPEC']
        # for key in 'url', 'database', 'user', 'pw':
        if item not in self.config["iSPEC"]:
            self.config["iSPEC"][item] = ""
        value = self.config["iSPEC"][item]
        if value == "":
            hide_input = True if item == "pw" else False
            help_text = self.helps.get(item, "")
            new_value = click.prompt(
                "No entry for {} found\n{}".format(item, help_text),
                hide_input=hide_input,
            )
            self.config["iSPEC"][item] = new_value
            return new_value
        else:
            return value

    def update_config(self):
        # diff = [x for x in self.config['iSPEC']]
        diff = {
            k: self.config["iSPEC"][k]
            for k in self.config["iSPEC"]
            if k in self.orig_params and self.config["iSPEC"][k] != self.orig_params[k]
        }
        if not diff:  # nothing to do
            return

        if self.save_all is None:
            save_all = click.confirm(
                "Save user/password? Not recommended on shared computers!"
            )
            self.save_all = save_all
        if not self.save_all:  # remove
            user, pw = self.config["iSPEC"]["user"], self.config["iSPEC"]["pw"]
            self.config["iSPEC"]["user"] = ""
            self.config["iSPEC"]["pw"] = ""

        with open(CONF_FILE, "w") as f:
            self.config.write(f)

        if not self.save_all:  # restore
            self.config["iSPEC"]["user"], self.config["iSPEC"]["pw"] = user, pw

        # now update the  "orig params" so we don't keep saving the same thing
        # over and over
        for k, v in diff.items():
            self.orig_params[k] = v

    # do we really need all these?
    @property
    def url(self):
        return self.get_item("url")
        # if self._url is None:
        #     self._url =  self.get_item('url')
        # return self._url

    @property
    def database(self):
        return self.get_item("database")
        # if self._database is None:
        #     self._database =  self.get_item('database')
        # return self._database

    @property
    def user(self):
        return self.get_item("user")
        # if self._user is None:
        #     self._user =  self.get_item('user')
        # return self._user

    @property
    def pw(self):
        return self.get_item("pw")
        # if self._pw is None:
        #     self._pw =  self.get_item('pw')
        # return self._pw

    @property
    def login_params(self):
        params = {attr: getattr(self, attr) for attr in self.attrs}
        # now call save:
        self.update_config()
        return params


login_params = Conf()

pd.set_option(
    "display.width",
    None,
    "display.max_colwidth",
    100,
    "display.max_columns",
    500,
)

try:
    import pyodbc
except ImportError:
    pass
from bcmproteomics.e2gitems import e2gcolumns, psm_columns, tmt_columns, itraq_columns

try:
    from bcmproteomics.classify import score_experiments

    _classify = True
except ImportError:
    _classify = False


user = None
pw = None
url = "10.16.2.74"
database = "iSPEC_BCM"
params = {"user": user, "pw": pw, "url": url, "database": database}


def reset_index_if_not_unique(df):
    if not df.index.is_unique:
        if (
            not df.index.name in df.columns
        ):  # don't throw away the index if we don't have it
            df[df.index.name] = df.index
        df = df.reset_index(drop=True)
        # warn(
        #     """Returned dataframe not indexed on GeneID due to duplicate records.
        # If this is a labeled experiment you can safely ignore this warning.
        # Note that joining multiple experiments will NOT work correctly."""
        # )
    return df


def _find_file(target, path):
    """Try to find a file in a given path

    :param target: target file name
    :param path: path name
    :returns: path to file . ext
    :rtype: str

    """

    dtype = ""
    if target.endswith("json"):
        dtype = "metadata"
    else:
        dtype = "e2g" if "e2g" in target else "psm"

    # result = [x for x in os.listdir(path) if x == target]

    # result = glob(os.path.join(path, target))
    globstr = os.path.join(path, "**", target)
    result = glob(globstr, recursive=True)
    # import ipdb

    # ipdb.set_trace()
    # prioritize dtype_QUAL and dtype_QUANT if present
    # eventually this will become default
    if any("_QUAL" in x for x in result) and any("_QUANT" in x for x in result):
        result3 = [
            x
            for x in result
            if dtype in x and ("QUAL" in x or "QUANT" in x)
            # if any(y in x for y in
            # (f"{dtype}_QUANT", f"{dtype}_QUAL")
            #  ["{}_QUANT".format(dtype), "{}_QUAL".format(dtype)])
        ]
        result3 = [x for x in result3 if "labelTMT_TMT_" not in x]
        if len(result3) == 2:
            return result3

    if result and len(result) == 1:
        # return os.path.abspath(os.path.join(path, result[0]))
        return result[0]
    elif result and len(result) > 1:  # more than 1 file, try to guess which one
        dtype = ""
        if target.endswith("json"):
            dtype = "metadata"
        else:
            dtype = "e2g" if "e2g" in target else "psm"
        result2 = [
            x for x in result if re.search(r"(?<!\d)_{}s?_\w+.tsv".format(dtype), x)
        ]
        # if dtype =='psm':
        #     import ipdb; ipdb.set_trace()
        if len(result2) == 1:
            return result2[0]
        # result3 = [x for x in result if any(y in x for y in ['{}_QUANT'.format(dtype), '{}_QUAL'.format(dtype)])]
        result3 = [x for x in result if any(y in x for y in ["_QUANT", "_QUAL"])]
        result3 = [x for x in result3 if "labelTMT_TMT_" not in x]
        result3 = [x for x in result3 if "_all.txt" not in x]

        if len(result3) == 2:
            return result3

        ret = sorted(result, key=len)[-1]
        warn("More than 1 file found, {}\nUsing{}".format(result, ret))
        return ret
    return None


class Experiment:
    """Container for accessing metadata for a given experiment.

    :param recno: int/str for record number
    :param runno: int/str, default 1
    :param searchno: int/str, default 1
    :param auto_populate: Whether or not to try to auto populate the data, default True
    :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
    :returns: Experiment instance
    :rtype: ispec.Experiment

    .. note::
        If data_dir is specified, data will *automatically* be saved in the given directory
    """

    def __init__(
        self,
        recno=None,
        runno=None,
        searchno=None,
        auto_populate=True,
        data_dir=None,
        only_local=False,
    ):
        """Loads metadata for a given experiment.

        :param recno: int/str for record number
        :param runno: int/str, default 1
        :param searchno: int/str, default 1
        :param auto_populate: Whether or not to try to auto populate the data, default True
        :param data_dir: (optional) data directory for saving/loading data
                         If specified, will first look in data_dir for data before making network calls
        :returns: Experiment instance
        :rtype: ispec.Experiment
        .. note::
            If data_dir is specified, data will **automatically** be saved in the given directory
        """

        try:
            recno = int(float(recno))
        except ValueError:
            raise ValueError("`recno` could not be coerced to int")

        if runno:
            try:
                runno = int(float(runno))
            except ValueError:
                raise ValueError("`runno` could not be coerced to int")

        if searchno:
            try:
                searchno = int(float(searchno))
            except ValueError:
                raise ValueError("`searchno` could not be coerced to int")

        self.recno = recno
        self.runno = runno or 1
        self.searchno = searchno or 1
        # self.techrepno =
        self.sample = ""
        self.treatment = ""
        self.exptype = ""
        self.genotype = ""
        self.added_by = ""
        self.identifiers = ""
        self.taxon_ratios = dict()
        self.labeltype = ""
        self.pygrouper_version = ""
        if data_dir is not None:
            data_dir = os.path.abspath(data_dir)
        self.data_dir = data_dir
        self._df = pd.DataFrame()  # else set as empty dataframe
        if recno is not None and auto_populate:
            if self.data_dir is not None:
                self.load_local_metadata(self.data_dir, only_local=only_local)
            elif not only_local:
                self.get_metadata(
                    recno, runno, searchno
                )  # self._df gets set through here
            else:
                print(
                    "No data directory provided and `only_local=True`, not querying iSPEC"
                )
        self._joined = False
        self.ibaq_normalize = None

    def __repr__(self):
        """Returns record and run number"""
        return "{}_{}_{}".format(self.recno, self.runno, self.searchno)
        # return 'Record number {}, run number {}'.format(self.recno, self.runno)

    # def __getattr__(self, name):
    #     if name not in self.__dict__ and hasattr(pd.DataFrame, name):
    #         return getattr(self.df, name)
    #     else:
    #         raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__, name))

<<<<<<< HEAD

    @property
    def metadata(self):
        return self._metadata

    @property
    def df(self):
        "Placeholder"
	# now in metadata
        return self._df

    @property
    def _database(self):
        return params.get("database")

    def reload(self):
        """Reload the dataset"""
        if self._joined:
            raise NotImplementedError("Cannot reload data for a joined experiment.")
        self.get_exprun(self.recno, self.runno, self.searchno)

    def get_metadata(self, recno, runno=1, searchno=1, label=0):
        if not recno:
            raise ValueError("recno must be specified")
        if runno is None:
            runno = 1  # make sure default runno is 1
        if searchno is None:
            searchno = 1
        self.recno = recno
        self.runno = runno

        conn = filedb_connect()
        database = ispec.login_params.database # should be str
        if isinstance(conn, str):
            return # failed to make connection, user is informed in filedb_connect()

        #col_str = "xp_EXPLabelFLAG, exp_LabelType"
        COLUMNS = [
        "exprun_EXPRecNo", "exprun_EXPRunNo", "exprun_EXPSearchNo",
        "exp_EXPLabelFLAG",
        "exp_EXPLabelType",
        "exp_Extract_No",
        "exp_Extract_CellTissue",
        "exp_Extract_Genotype",
        "exp_ExpType",
        "exp_DisplayID",
        "exp_Exp_Experimenter",
        "exp_Exp_Date",
        "exp_Extract_Treatment",
        "exp_Extract_Amount",
        "exp_Extract_Fractions",
        "exp_Extract_TaxonID_lu",
        "exp_ProtocolNo",
        "exp_Washes",
        "exp_Separation_1",
        "exp_Separation_1Detail",
        "exp_Separation_2",
        "exp_Separation_2Detail",
        # "expSeparation_3",
        # "exp_Separation_3Detail",
        "exprun_TaxonID",
        "exprun_Search_EnzymeSetting",
        "exprun_MS_Instrument",
        "exprun_MS_Experimenter",
        "exprun_Grouper_RefDatabase",
        "exprun_search_RefDatabase",
        ]
        col_str = str.join(", ", COLUMNS)
        #col_str = "exprun_EXPRecNo, exprun_EXPRunNo, exprun_EXPSearchNo "

        sql_query =("SELECT",
            col_str,
            f"from {database}.iSPEC_ExperimentRuns INNER JOIN {database}.iSPEC_Experiments on",
            f"{database}.iSPEC_Experiments.exp_EXPRecNo = {database}.iSPEC_Experimentruns.exprun_EXPRecNo",
	    "WHERE ",
        f"exprun_EXPRecNo={recno} AND",
        f"exprun_EXPRunNo={runno} AND",
        f"exprun_EXPSearchNo={searchno}"
	    # f"{database}.iSPEC_ExperimentRuns.exprun_EXPRecNo={recno} ",
	    # f"AND {database}.iSPEC_ExperimentRuns.exprun_EXPRunNo={runno} ",
	    # f"AND {database}.iSPEC_ExperimentRuns.exprun_EXPSearchNo={searchno} ",
            )
            # f"{database}.iSPEC_Experiments.exprun_EXPRecNo={recno} ",
        # database = self._database
        # #  or database = ispec.login_params.database
        # testing 
        logger.info(f"sql_query: {sql_query}")
        df = pd.read_sql(str.join(" ",sql_query), conn)
        self._metadata = df
	# this is our "replacement" metadata to supersede what comes next
        #info = df.to_dict("list")
        # testing 
        #

        #sql_description = ("SELECT exp_EXPClass, exp_Extract_CellTissue, exp_Exp_Description,"
        #                   "exp_Extract_Treatment, exp_IDENTIFIER, exp_Extract_No, "
        #                   "exp_Extract_Genotype, exp_AddedBy, exp_Digest_Type, exp_Digest_Enzyme "
        #                   "from {database}.iSPEC_Experiments where exp_EXPRecNo={recno}").format(database=self._database,
        #                                                                                          recno=recno)

        sql_description = (f"SELECT exp_EXPClass, exp_Extract_CellTissue, exp_Exp_Description,"
                           f"exp_Extract_Treatment, exp_IDENTIFIER, exp_Extract_No, "
                           f"exp_Extract_Genotype, exp_AddedBy, exp_Digest_Type, exp_Digest_Enzyme "
                           f"from {database}.iSPEC_Experiments where exp_EXPRecNo={recno}")
        info = pd.read_sql(sql_description, conn).to_dict('list') # a 1 row dataframe

        #import ipdb; ipdb.set_trace()
        info = df.to_dict("list")
        self.sample = ''.join(item for item in info.get('exp_Extract_CellTissue', '') if item)
        self.treatment = ''.join(item for item in info.get('exp_Extract_Treatment', '') if item)
        self.exptype = ''.join(item for item in info.get('exp_EXPClass', '') if item)
        self.description = [item for items in info.get('exp_Exp_Description', '')
                            for item in (items.splitlines() if items else '')]
        self.genotype = ''.join(item for item in info.get('exp_Extract_Genotype', '') if item)
        self.added_by = ''.join(item for item in info.get('exp_AddedBy', '') if item)
        self.digest_type = ''.join(item for item in info.get('exp_Digest_Type', '') if item)
        self.digest_enzyme = ''.join(item for item in info.get('exp_Digest_Enzyme', '') if item)
        # import ipdb; ipdb.set_trace()
        _extract = info.get('exp_Extract_No', 0)
        if isinstance(_extract, list):
            _extract = set(_extract)
        # self.extract_no = ''.join(str(item) for item in set(info.get('exp_Extract_No', [0])) if item)
        #self.extract_no = ''.join(str(item) for item in _extract if item)
        self.extract_no = f"{_extract}",
        self.identifiers = ''.join(item for item in info.get('exp_IDENTIFIER', '') if item).splitlines()
            return  # failed to make connection, user is informed in filedb_connect()

        sql_description = (
            "SELECT exp_EXPClass, exp_Extract_CellTissue, exp_Exp_Description,"
            "exp_Extract_Treatment, exp_IDENTIFIER, exp_Extract_No, "
            "exp_Extract_Genotype, exp_AddedBy, exp_Digest_Type, exp_Digest_Enzyme "
            "from {database}.iSPEC_Experiments where exp_EXPRecNo={recno}"
        ).format(database=self._database, recno=recno)
        info = pd.read_sql(sql_description, conn).to_dict("list")  # a 1 row dataframe
        self.sample = "".join(
            item for item in info.get("exp_Extract_CellTissue", "") if item
        )
        self.treatment = "".join(
            item for item in info.get("exp_Extract_Treatment", "") if item
        )
        self.exptype = "".join(item for item in info.get("exp_EXPClass", "") if item)
        self.description = [
            item
            for items in info.get("exp_Exp_Description", "")
            for item in (items.splitlines() if items else "")
        ]
        self.genotype = "".join(
            item for item in info.get("exp_Extract_Genotype", "") if item
        )
        self.added_by = "".join(item for item in info.get("exp_AddedBy", "") if item)
        self.digest_type = "".join(
            item for item in info.get("exp_Digest_Type", "") if item
        )
        self.digest_enzyme = "".join(
            item for item in info.get("exp_Digest_Enzyme", "") if item
        )
        self.extract_no = "".join(
            str(item) for item in info.get("exp_Extract_No", 0) if item
        )
        self.identifiers = "".join(
            item for item in info.get("exp_IDENTIFIER", "") if item
        ).splitlines()
        self.recno = recno
        self.runno = runno
>>>>>>> 9af07b5b69a96860bf806bf7fa45c6ac13545381
        self.searchno = searchno
        if self.extract_no:
            extract_fractions, extract_protocol = self.get_extract_data(self.extract_no)
            self.extract_fractions = extract_fractions
            self.extract_protocol = extract_protocol

        additional_info = (
            "SELECT exprun_Fraction_9606, exprun_Fraction_10090, exprun_Fraction_9031, "
            "exprun_LabelType, exprun_Grouper_Version "
            "from {database}.iSPEC_ExperimentRuns where "
            "exprun_EXPRecNo={recno} "
            "AND exprun_EXPRunNo={runno} "
            "AND exprun_EXPSearchNo={searchno}"
        ).format(recno=recno, runno=runno, searchno=searchno, database=self._database)
        cursor = conn.cursor()
        cursor.execute(additional_info)
        additional_info_result = cursor.fetchone()
        if additional_info_result:
            hu, mou, gg, labeltype, pygrouper_version = additional_info_result
            self._add_taxon_ratios(hu, mou, gg)
            self._assign_labeltype(labeltype)
            self.pygrouper_version = pygrouper_version

        return self

    def _assign_labeltype(self, labeltype=None):
        self.labeltype = labeltype
        return self

    def populate_metadata(self, info: dict):
        """Populate metadata from a dictionary originating from json dump
        of original data"""
        self.sample = info.get("sample")
        self.treatment = info.get("treatment")
        self.exptype = info.get("exptype")
        self.description = info.get("description", [])
        self.genotype = info.get("genotype")
        self.added_by = info.get("added_by")
        self.digest_type = info.get("digest_type")
        self.digest_enzyme = info.get("digest_enzyme")
        self.extract_no = info.get("extract_no", 0)
        self.identifiers = info.get("identifiers")
        self.extract_fractions = info.get("extract_fractions")
        self.extract_protocol = info.get("extract_protocol")
        self.taxon_ratios = info.get("taxon_ratios", dict())
        return self

    def _add_taxon_ratios(self, hu, mou, gg):
        for taxa, id in ((hu, "9606"), (mou, "10090"), (gg, "9031")):
            self.taxon_ratios[id] = taxa
        return self

    def get_extract_data(self, extractno):
        """Extract"""
        conn = filedb_connect()
<<<<<<< HEAD
        sql_extractdata = ("SELECT ext_Fractions, ext_Protocol FROM "
                           "{database}.iSPEC_Extracts "
                           "WHERE ext_ExtRecNo={extract_no}").format(database=self._database, extract_no=self.extract_no)
        try:
            extract_info = pd.read_sql(sql_extractdata, conn).to_dict('list') # a 1 row dataframe
        except pd.io.sql.DatabaseError:
            return '', ''
=======
        sql_extractdata = (
            "SELECT ext_Fractions, ext_Protocol FROM "
            "{database}.iSPEC_Extracts "
            "WHERE ext_ExtRecNo={extract_no}"
        ).format(database=self._database, extract_no=self.extract_no)
        extract_info = pd.read_sql(sql_extractdata, conn).to_dict(
            "list"
        )  # a 1 row dataframe
>>>>>>> 9af07b5b69a96860bf806bf7fa45c6ac13545381
        try:
            extract_fractions = "".join(extract_info.get("ext_Fractions", ""))
        except TypeError:
            extract_fractions = "None"
        try:
            extract_protocol = "".join(extract_info.get("ext_Protocol", "")).replace(
                "\r", ", "
            )
        except TypeError:
            extract_protocol = "None"
        conn.close()
        return (extract_fractions, extract_protocol)

    @property
    def _database(self):
        return params.get("database")

    @property
    def json_metadata(self):
        """Metadata to json

        :returns: json
        :rtype: str

        """
        json_data = json.dumps(
            {
                attr: self.__getattribute__(attr)
                for attr in self.__dict__.keys()
                if not attr.startswith("_")
            }
        )
        return json_data

    def save(self, data_dir=None):
        """Save data to given directory. Called automatically if data_dir is provided
        and data not present locally

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: None
        :rtype: NoneType

        """
        if data_dir is None:
            data_dir = self.data_dir
        with open(os.path.join(data_dir, "{!r}.json".format(self)), "w") as metaf:
            json.dump(self.json_metadata, metaf)
        return

    def load_local_metadata(self, data_dir=None, only_local=False):
        """Try to load the data from disk rather than over network.
        This method is usually invoked automatically and does not usually
        need to be manually invoked.

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: self
        :rtype: ispec.Experiment

        """
        if self.runno is None:
            raise ValueError("recno must be specified")
        if data_dir is None:
            data_dir = self.data_dir
        target = _find_file(target="{!r}*.json".format(self), path=data_dir)
        if target is None and not only_local:
            self.get_metadata(self.recno, self.runno, self.searchno)
            self.save(data_dir)
            return self
        elif target is None and only_local:
            return self
        # try:
        #     metadata = json.load(open(target, 'r'))
        # except:
        try:
            metadata = json.loads(
                json.load(open(target, "r"))
            )  # not sure why this is necessary sometimes
        except TypeError:
            metadata = json.load(open(target))
            # metadata = json.loads(open(target).read())
        self.populate_metadata(metadata)
        return self


class PSMs(Experiment):
    """Container for accessing psms data for a given experiment

    :param recno: int/str for record number
    :param runno: int/str, default 1
    :param searchno: int/str, default 1
    :param auto_populate: Whether or not to try to auto populate the data, default True
    :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
    :returns: PSMs instance
    :rtype: ispec.PSMs

    .. seealso:: ispec.Experiment
    .. todo:: add presplit support

    """

    def __init__(
        self,
        recno=None,
        runno=None,
        searchno=None,
        presplit=False,
        auto_populate=True,
        data_dir=None,
        only_local=False,
    ):
        """Container for accessing psms data for a given experiment

        :param recno: int/str for record number
        :param runno: int/str, default 1
        :param searchno: int/str, default 1
        :param auto_populate: Whether or not to try to auto populate the data, default True
        :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
        :returns: PSMs instance
        :rtype: ispec.PSMs

        .. seealso:: ispec.Experiment
        .. todo:: add presplit support

        """
        super().__init__(
            recno=recno,
            runno=runno,
            searchno=searchno,
            auto_populate=auto_populate,
            data_dir=data_dir,
            only_local=only_local,
        )
        self.presplit = presplit
        ##
        ##
        if recno is not None and auto_populate:
            if self.data_dir is not None:
                self.load_local(self.data_dir, only_local=only_local)
            elif not only_local:
                self.get_exprun(
                    recno, self.runno, self.searchno
                )  # self._df gets set through here
            else:
                print(
                    "No data directory provided and `only_local=True`, not querying iSPEC"
                )
        if self._df is None or len(self.df) == 0:
            self._df = pd.DataFrame(
                columns=[x.split("_")[1] for x in psm_columns]
            )  # else set as empty dataframe
        self._joined = False

    def get_exprun(self, recno, runno, searchno):
        """queries iSPEC database and grabs the gene product table which is stored in self.df"""
        conn = filedb_connect()
        if isinstance(conn, str):
            return  # failed to make connection, user is informed in filedb_connect()
        _psm_columns = psm_columns.copy()
        if "tmt" in self.labeltype.lower():
            _psm_columns += tmt_columns
        elif "itraq" in self.labeltype.lower():
            _psm_columns += itraq_columns
        print(self.labeltype)
        print("\n", _psm_columns, "\n")
        sql = (
            "SELECT {psm_cols} from {database}.iSPEC_data where "
            "psm_EXPRecNo={recno} "
            "AND psm_EXPRunNo={runno} "
            "AND psm_EXPSearchNo={searchno}"
        ).format(
            psm_cols=", ".join(_psm_columns),
            recno=recno,
            runno=runno,
            searchno=searchno,
            database=self._database,
        )
        if self.presplit:
            sql += "\nAND psm_oriFLAG=1"

        self._df = self._construct_df(sql, conn)
        conn.close()
        return self

    @staticmethod
    def _construct_df(sql, conn):
        """Construct a pandas dataframe from data in the iSPEC."""
        df = pd.read_sql(
            sql,
            conn,
        )
        df.rename(
            columns={
                k: k.split("psm_")[1]
                for k in [col for col in df.columns if col.startswith("psm_")]
            },
            inplace=True,
        )

        return df

    def reload(self):
        """Reload the dataset"""
        if self._joined:
            raise NotImplementedError("Cannot reload data for a joined experiment.")
        self.get_exprun(self.recno, self.runno, self.searchno)

    def save(self, data_dir=None):
        """Save data to given directory. Called automatically if data_dir is provided
        and data not present locally

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: None
        :rtype: NoneType

        """
        if data_dir is None:
            data_dir = self.data_dir
        super().save(data_dir=data_dir)
        if len(self.df) == 0:  # don't save if no data!
            return
        with open(os.path.join(data_dir, "{!r}_psms.tab".format(self)), "w") as data:
            self.df.to_csv(data, sep="\t", index=False)

    def load_local(self, data_dir=None, only_local=False):
        """Try to load the data from disk rather than over network.

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: self
        :rtype: ispec.PSMs

        """
        self.load_local_metadata(data_dir=data_dir, only_local=only_local)
        # psmsfile = _find_file(target='{!r}*psms*t[sa][vb]'.format(self), path=data_dir)
        psmsfile = _find_file(target="{!r}*psm*t*".format(self), path=data_dir)

        # TODO fix this for QUAL and QUANT
        ## this is now done within the `_find_files` func.
        # if len(psmsfile) > 2:
        #     psmsfile = [x for x in psmsfile if "psms_all.txt" not in x]
        #     # hack to remove "psms_all" combined file when "psms_QUAL" and "psms_QUANT" exist
        if psmsfile is None and not only_local:

            self.get_exprun(self.recno, self.runno, self.searchno)
            self.save(data_dir)
        elif psmsfile is None and only_local:
            return self  # don't query iSPEC
        else:
            if len(psmsfile) == 2:

                quant = [x for x in psmsfile if "QUANT" in x]
                qual = [x for x in psmsfile if "QUAL" in x]
                assert len(quant) == len(qual) == 1
                _quant_df = pd.read_table(quant[0], dtype={"GeneID": "str"})
                _quant_df = _quant_df[~_quant_df.PeptRank.isna()]
                _qual_df = pd.read_table(qual[0], dtype={"GeneID": "str"})
                # how does this work for SILAC?
                # not_sra = [x for x in _qual_df.columns if x !='SRA' and x != 'LabelFLAG']
                _qual_cols = [
                    "RTmin",
                    "Charge",
                    "FirstScan",
                    "IonScore",
                    "Modifications",
                    "SequenceModi",
                    "EXPRecNo",
                    "EXPRunNo",
                    "EXPSearchNo",
                    "GeneID",
                    "SpectrumFile",
                    "PSM_UseFLAG",
                    "Sequence",
                ]

                self._df = _quant_df.merge(
                    _qual_df[_qual_cols],
                    on=[
                        "EXPRecNo",
                        "EXPRunNo",
                        "EXPSearchNo",
                        "SpectrumFile",
                        "GeneID",
                        "Charge",
                        "FirstScan",
                        "SequenceModi",
                    ],
                    indicator=True,
                    how="left",
                )

            # alternative way, not as clear
            # if len(psmsfile) == 1:
            #     self._df = pd.read_table(psmsfile)
            # elif len(psmsfile) == 2:
            #     _qual = list(filter(lambda x: "QUAL" in x, psmsfile))[0]
            #     _quant = list(filter(lambda x: "QUANT" in x, psmsfile))[0]
            #     _df1 = pd.read_table(_qual)
            #     _df2 = pd.read_table(_quant)
            #     _toignore = (
            #         "LabelFLAG",
            #         "SequenceModi",
            #         "PrecursorArea",
            #         "PrecursorArea_dstrAdj",
            #         "PrecursorArea_split",
            #         "SequenceArea",  # this is ultimately what we want
            #         # "SequenceModi",
            #     )
            #     _df = _df2.merge(
            #         # _df1[[x for x in _df1 if x not in ("LabelFLAG", "SequenceModi")]],
            #         _df1[[x for x in _df1 if x not in _toignore]],
                del _quant_df
                del _qual_df

            else:
                self._df = pd.read_table(psmsfile)
        return self


class E2G(Experiment):
    """Container for accessing gene product data for a given experiment

    :param recno: int/str for record number
    :param runno: int/str, default 1
    :param searchno: int/str, default 1
    :param auto_populate: Whether or not to try to auto populate the data, default True
    :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
    :returns: E2G instance
    :rtype: ispec.E2G

    .. seealso:: ispec.Experiment

    """

    def __init__(
        self,
        recno=None,
        runno=None,
        searchno=None,
        auto_populate=True,
        data_dir=None,
        only_local=False,
    ):
        """Container for accessing gene product data for a given experiment

        :param recno: int/str for record number
        :param runno: int/str, default 1
        :param searchno: int/str, default 1
        :param auto_populate: Whether or not to try to auto populate the data, default True
        :param data_dir: (optional) data directory for saving/loading data
                     If specified, will first look in data_dir for data before making network calls
        :returns: E2G instance
        :rtype: ispec.E2G

        .. seealso:: ispec.Experiment

        """
        # """Different metadata as well as data"""
        super().__init__(
            recno=recno,
            runno=runno,
            searchno=searchno,
            auto_populate=auto_populate,
            data_dir=data_dir,
            only_local=only_local,
        )
        if recno is not None and auto_populate:
            if data_dir is not None:
                self.load_local(self.data_dir, only_local=only_local)
            elif not only_local:
                self.get_exprun(
                    recno, self.runno, self.searchno
                )  # self._df gets set through here
            else:
                print(
                    "No data directory provided and `only_local=True`, not querying iSPEC"
                )
        if self._df is None or len(self.df) == 0:
            self._df = pd.DataFrame(
                columns=[x.split("_")[1] for x in e2gcolumns]
            )  # else set as empty dataframe

        self._joined = False
        self.ibaq_normalize = None
        if self._df.empty:
            return
        if "SRA" not in self.df:
            try:
                self.assign_sra(self.df)
            except ValueError:
                print("Could not assign SRA for {}".format(self))

    def __add__(self, other):
        return join_exps(self, other)

    def __len__(self):
        return len(self.df)

    @property
    def df(self):
        """Data for the experiment

        :returns: df
        :rtype: pandas.DataFrame

        """
        return self._df

    def reload(self):
        """Reload the dataset"""
        if self._joined:
            raise NotImplementedError("Cannot reload data for a joined experiment.")
        self.get_exprun(self.recno, self.runno, self.searchno)

    def get_exprun(self, recno=None, runno=1, searchno=1):
        """queries iSPEC database and grabs the gene product table which is stored in self.df"""
        conn = filedb_connect()
        if isinstance(conn, str):
            return  # failed to make connection, user is informed in filedb_connect()

        sql = (
            "SELECT {e2gcols} from {database}.iSPEC_exp2gene where "
            "e2g_EXPRecNo={recno} "
            "AND e2g_EXPRunNo={runno} "
            "AND e2g_EXPSearchNo={searchno}"
        ).format(
            e2gcols=", ".join(e2gcolumns),
            recno=recno,
            runno=runno,
            searchno=searchno,
            database=self._database,
        )

        self._df = self._construct_df(sql, conn)
        conn.close()
        return self

    @staticmethod
    def assign_sra(df):
        df["SRA"] = "A"
        df["SRA"] = df["SRA"].astype("category", categories=("S", "R", "A"))

        df.loc[(df["IDSet"] == 1) & (df["IDGroup_u2g"] <= 3), "SRA"] = "S"

        df.loc[(df["IDSet"] == 2) & (df["IDGroup"] <= 3), "SRA"] = "S"

        df.loc[
            (df["IDSet"] == 1) & (df["SRA"] != "S") & (df["IDGroup_u2g"] <= 5), "SRA"
        ] = "R"

        df.loc[
            (df["IDSet"] == 2) & (df["SRA"] != "S") & (df["IDGroup"] <= 5), "SRA"
        ] = "R"

    @staticmethod
    def _construct_df(sql, conn):
        """Construct a pandas dataframe from data in the iSPEC."""
        df = pd.read_sql(sql, conn, index_col="e2g_GeneID")
        df.index.rename("GeneID", inplace=True)
        df["GeneID"] = df.index.astype("object")
        df.index = df.index.astype("object")
        df.rename(
            columns={
                k: k.split("e2g_")[1]
                for k in [e2gcol for e2gcol in df.columns if e2gcol.startswith("e2g_")]
            },
            inplace=True,
        )
        df.rename(
            columns={
                "n_iBAQ_dstrAdj": "iBAQ_dstrAdj",
            },
            inplace=True,
        )

        # df.rename(columns={'e2g_GeneID':'GeneID'}, inplace=True)
        geneidlist = [str(x) for x in df.GeneID.tolist() if not np.isnan(x)]
        genesql = (
            "Select gene_GeneID, "
            "gene_GeneSymbol, gene_GeneDescription, "
            "gene_FunCats "
            "from iSPEC_BCM.iSPEC_Genes "
            "where gene_GeneID in ({})".format(", ".join(geneidlist))
        )
        if geneidlist:
            genedf = pd.read_sql(
                genesql, conn, index_col="gene_GeneID"
            )  # all headers start with gene_
            generename = {c: c.split("gene_")[1] for c in genedf.columns}
            genedf.rename(columns=generename, inplace=True)
            genedf.index.rename("GeneID", inplace=True)
            df.index = df.index.astype("object")
            # genedf.rename(columns={'u2gPeptIBAQAveCount':'GeneCapacity'}, inplace=True)
            # df['e2g_GeneCapacity'] = df['e2g_nGPArea_Sum_dstrAdj']/df['e2g_nGPArea_Sum_dstrAdj']
            # funcat_table = pd.read_table('E:/reference_databases/iSPEC_genes_Aug_2015.tab')
            # df = pd.merge(df, funcat_table, how='left', on='GeneID')
            # df = pd.merge(df, genedf, on='GeneID')
            df = df.join(genedf)
            df["FunCats"].fillna("", inplace=True)
        if "GeneID" not in df.columns:
            df = reset_index_if_not_unique(df)
        return df

    def save(self, data_dir=None):
        """Save data to given directory. Called automatically if data_dir is provided
        and data not present locally

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: None
        :rtype: NoneType

        """
        if data_dir is None:
            data_dir = self.data_dir
        super().save(data_dir=data_dir)
        if len(self.df) == 0:  # don't save if no data!
            return
        with open(os.path.join(data_dir, "{!r}_e2g.tab".format(self)), "w") as data:
            # self.df.to_csv(data, sep='\t', index=True if 'GeneID' == self.df.index.name else False)
            self.df.to_csv(data, sep="\t", index=False)

    def load_local(self, data_dir=None, only_local=False):
        """Try to load the data from disk rather than over network.
        This method is usually invoked automatically and does not usually
        need to be manually invoked.

        :param data_dir: (Optional) Directory to save data.
                         Uses self.data_dir if not specified here.
        :returns: self
        :rtype: ispec.E2G

        """
        self.load_local_metadata(data_dir=data_dir, only_local=only_local)
        # e2gfile = _find_file(target='{!r}*_e2g*t[sa]|[vb]'.format(self), path=data_dir)
        e2gfile = _find_file(target="{!r}*_e2g*t*".format(self), path=data_dir)
        if e2gfile is None and not only_local:
            self.get_exprun(self.recno, self.runno, self.searchno)
            self.save(data_dir)
        elif e2gfile is None and only_local:
            return self  # don't query iSPEC
        else:
            if len(e2gfile) == 2:

                quant = [x for x in e2gfile if "QUANT" in x]
                qual = [x for x in e2gfile if "QUAL" in x]
                assert len(quant) == len(qual) == 1

                _quant_df = pd.read_table(quant[0], dtype={"GeneID": "str"})
                _qual_df = pd.read_table(qual[0], dtype={"GeneID": "str"})
                # how does this work for SILAC?
                not_sra = [
                    x for x in _qual_df.columns if x != "SRA" and x != "LabelFLAG"
                ]

                self._df = _quant_df.merge(
                    _qual_df[not_sra],
                    on=["EXPRecNo", "EXPRunNo", "EXPSearchNo", "GeneID"],
                )
                del _quant_df
                del _qual_df

            else:
                _df = pd.read_table(e2gfile, index_col="GeneID", engine="c")
                _df = _df.rename(columns={c: c.split("e2g_")[-1] for c in _df})
                # self._df = pd.read_table(e2gfile, index_col='GeneID', engine='c')
                self._df = _df
                self._df = reset_index_if_not_unique(self.df)
            if "FunCats" not in self._df:
                self._df["FunCats"] = ""
            self._df["FunCats"] = self.df["FunCats"].fillna("")
            if self.df.index.name == "GeneID" and "GeneID" not in self.df.columns:
                self.df["GeneID"] = self.df.index
        return self

    def strict(self, df=None, set1=False):
        """Returns a strict selection of gene products from the dataframe
        Attributes
        ----------
        df : desired dataframe to be filtered (optional, default is E2G.df)

        set1 : filter even more stringently with only IDSet 1 (optional, default False)
        """

        if df is None:
            df = self.df
        if not set1:
            cutoff = 3
        elif set1:
            cutoff = 2

        if self._joined:
            return df[
                (
                    (df["IDSet_x"] < cutoff) & (df["IDGroup_x"] <= 3)
                    | (df["IDSet_y"] < cutoff) & (df["IDGroup_y"] <= 3)
                )
            ]

        return df[(df["IDSet"] < cutoff) & (df["IDGroup"] <= 3)]

    def relaxed(self, df=None, set1=False):
        """Returns a 'relaxed' selection of gene products from the dataframe
        Attributes
        ----------
        df : desired dataframe to be filtered (optional, default is E2G.df)

        set1 : filter even more stringently with only IDSet 1 (optional, default False)
        """

        if df is None:
            df = self.df
        if not set1:
            cutoff = 3
        elif set1:
            cutoff = 2

        if self._joined:
            return df[
                (
                    (df["IDSet_x"] < cutoff) & (df["IDGroup_x"] <= 5)
                    | (df["IDSet_y"] < cutoff) & (df["IDGroup_y"] <= 5)
                )
            ]

        return df[(df["IDSet"] < cutoff) & (df["IDGroup"] <= 5)]

    @property
    def tfs_coRs(self):
        """Return gene products annotated as a DBTF or CoRegulator"""
        return self.df[self.df["FunCats"].str.contains("DBTF|TC|DBTC|CBTC")]

    @property
    def tfs(self):
        """Return gene products annotated as a DBTF"""
        return self.df[self.df["FunCats"].str.contains("DBTF")]

    @property
    def kinases(self):
        """Return gene products annotated as a kinase"""
        return self.df[self.df["FunCats"].str.contains("KI")]

    def USD_selector(self, df=None):
        """Returns U, S, D DataFrames for a joined E2G instance"""
        if self._joined is False:
            raise ValueError("Not a joined E2G instance.")
        if df is None:
            df = self.df
        ret = list()
        for cat in list("USD"):
            ret.append(df[df.USD == cat])
        return ret

    def summary(self):
        """Print out a quick summary of the types of gene products observed in the experiment."""

        if len(self._df) == 0:
            raise ValueError("DataFrame is empty!")

        searchstr = OrderedDict(
            {
                "TFs": "DBTF",
                "Tfs and CoRs": "DBTF|TC|DBTC|CBTC",
                "Kinases": "KI",
            }
        )
        print(
            "\nSummary for record {}, run {}.\n".format(
                self.recno,
                self.runno,
            )
        )

        if self._joined:
            dfU, dfS, dfD = self.USD_selector()
            print(
                "Total gene products found : "
                "{} same, {} up, {} down".format(len(dfS), len(dfU), len(dfD))
            )
            print(
                "Total strict gene products found : "
                "{} same, {} up, {} down".format(
                    len(self.strict(dfS)), len(self.strict(dfU)), len(self.strict(dfD))
                )
            )
            print("-" * 15, "\n")
            for s in searchstr:
                df_subsection = self.df[(self.df["FunCats"].str.contains(searchstr[s]))]
                dfU, dfS, dfD = self.USD_selector(df_subsection)
                print(
                    "Total {} gene products found : "
                    "{} same, {} up, {} down".format(s, len(dfS), len(dfU), len(dfD))
                )
                print(
                    "Total {} strict gene products found : "
                    "{} same, {} up, {} down".format(
                        s,
                        len(self.strict(dfS)),
                        len(self.strict(dfU)),
                        len(self.strict(dfD)),
                    )
                )
                print("-" * 15, "\n")
            return

        print(
            "Total gene products found : {}"
            "\n\tiBAQ total : {:4f}.".format(len(self._df), self._df.iBAQ_dstrAdj.sum())
        )

        df_strict = self.strict()
        print(
            "Total strict gene products found : {}"
            "\n\tiBAQ total : {:4f}.".format(
                len(df_strict), df_strict.iBAQ_dstrAdj.sum()
            )
        )
        print("-" * 15, "\n")
        for s in searchstr:
            df_subsection = self.df[(self.df["FunCats"].str.contains(searchstr[s]))]
            print(
                "Total {} gene products found : {}"
                "\n\tiBAQ total : {:4f}".format(
                    s, len(df_subsection), df_subsection.iBAQ_dstrAdj.sum()
                )
            )

            df_sub_strict = self.strict(df_subsection)
            print(
                "Total {} strict gene products found : {}"
                "\n\tiBAQ total : {:4f}.".format(
                    s, len(df_sub_strict), df_sub_strict.iBAQ_dstrAdj.sum()
                )
            )
            print("-" * 15, "\n")


# def _display(iterable):
#     """Display an iterable of things"""
#     mapping = dict()
#     for ix, element in enumerate(iterable):
#         mapping[ix+1] = element
#         print("({}) -- {}".format(ix+1,
#                                   element))
#     return mapping
# def _make_selection(iterable):
#     selected = None
#     value = None
#     mapping = _display(iterable)
#     while not selected or not value:
#         selected = click.prompt('Make a selection from above', type=int)
#         try:
#             value = mapping[selected]
#         except (KeyError, ValueError):
#             print("Invalid Selection\n")
#     return value

# def _getlogin():
#     """Checks if the username and password have been established for the current python session.
#     If username or password is undefined, will prompt the user.

#     Defaults to bcmproteomics.
#     """
#     servers = OrderedDict({'bcmproteomics': '10.16.2.74',
#                            'jun lab': '10.13.14.171',
#     })
#     databases = {'10.16.2.74': ['iSPEC_BCM', 'iSPEC_BCM_psms', 'iSPEC_BCM_IDG'],
#                  '10.13.14.171': ['iSPEC'],
#                  }
#     if params.get('user') is None:
#         print('Username is not set')
#         user = click.prompt('Enter your iSPEC username')
#         params['user'] = user
#     if params.get('pw') is None:
#         print('Password is not set')
#         pw = click.prompt('Enter your password', hide_input=True)
#         params['pw'] = pw
#     if params.get('url') is None:
#         print('iSPEC is not set, the options are:')
#         server = _make_selection(servers)
#         # print(*[server for server in servers], sep='\n')
#         # server = input('Select an iSPEC : ').strip()
#         params['url'] = servers.get(server, '10.16.2.74')
#     # elif 'url' in params:
#     #     params['url'] = servers.get(params['url'], '10.16.2.74')

#     if params.get('database') is None:
#         server_url = params['url']
#         if len(databases.get(server_url, [])) == 1:
#             params['database'] = databases[server_url][0]
#         else:
#             print('iSPEC database is not set, the options are:')
#             db = _make_selection(databases[server_url])
#             # print(*[database for database in databases[server_url]], sep='\t')
#             # db = input('Select an iSPEC database: ').strip()
#             params['database'] = databases.get(db, 'iSPEC_BCM')
#     return params


def filedb_connect(params=None):
    """Establish a connection to the BCM Proteomics FileMaker Server.

    Returns a database connection, which can then be used to execute sql statements.

    If the password has not been established yet, user will be prompted.
    Password is only stored in memory for the duration of the python session.
    """

    # params = _getlogin()
    # config = Conf()
    # params = config.login_params
    if not params:
        params = login_params.login_params

    driver = "DRIVER={FileMaker ODBC};"

    login_info = "UID={user};PWD={pw}".format(**params)

    server_info = "SERVER={url};".format(**params)

    database_info = "DATABASE={database};".format(**params)

    conn_string = "{0}{1}{2}{3}".format(driver, server_info, database_info, login_info)
    try:
        conn = pyodbc.connect(conn_string)
    except pyodbc.Error as e:  # any ODBC error is returned to python as "Error"
        error = e.__str__()
        if "password incorrect" in error.lower():
            print("Invalid password.")
        elif "account" in error.lower():
            print("Incorrect account.")
        elif "failed to connect" in error.lower():
            print("Error making connection, check the address.")
        else:
            print("Error with username or password")  # maybe raise an error instead?
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
        raise TypeError("Input must be a list with at least 1 element")
    # need to make sure every element in the geneidlist is a string
    conn = filedb_connect()
    if isinstance(conn, str):
        return conn  # failed to make connection, user is informed in filedb_connect()

    genesql = (
        "Select gene_GeneID, gene_u2gPeptiBAQAveCount, "
        "gene_GeneSymbol, gene_GeneDescription, "
        "gene_FunCats "
        "from iSPEC_BCM.iSPEC_Genes "
        "where gene_GeneID in ({})".format(", ".join([str(x) for x in geneidlist]))
    )

    genedf = pd.read_sql(
        genesql, conn, index_col="gene_GeneID"
    )  # all headers start with gene_
    generename = {c: c.split("gene_")[1] for c in genedf.columns}
    genedf.rename(columns=generename, inplace=True)
    genedf.index.rename("GeneID", inplace=True)
    genedf.rename(columns={"u2gPeptIBAQAveCount": "GeneCapacity"}, inplace=True)
    genedf["GeneID"] = [str(int(x)) for x in genedf.index.tolist()]
    genedf["FunCats"].fillna("", inplace=True)
    conn.close()
    return genedf


def get_geneids(taxonid):
    """Get all gene ids for a given taxon"""
    if not isinstance(taxonid, int):
        try:
            taxonid = int(taxonid.strip())
        except:
            raise TypeError("Input must be a taxon id")

    conn = filedb_connect()
    if isinstance(conn, str):
        return conn  # failed to make connection, user is informed in filedb_connect()
    sql = (
        "SELECT gene_GeneID from iSPEC_BCM.iSPEC_Genes "
        "WHERE gene_TaxonID={}".format(taxonid)
    )
    cursor = conn.execute(sql)
    gids = cursor.fetchall()
    genes = [int(g) for gid in gids for g in gid]
    return genes


def get_all_funcats(taxonid, data_dir=None):
    """Get all the gene information for a given taxonid.

    Optional data_dir argument will first look for output file
    in given directory.
    If not found, will make network call and save results there."""
    fname = "geneinfo_{}.tab".format(taxonid)

    if data_dir:
        f = _find_file(fname, data_dir)
        if f is not None:
            return pd.read_table(f, index_col="GeneID")

    gids = get_geneids(taxonid)
    funcats = get_funcats(gids)
    if data_dir:
        funcats.to_csv(os.path.join(data_dir, fname), sep="\t")
    return funcats


def align(exps, metadata=None):
    """
    exps is a collection of ispec.E2Gs or ispec.PSMs

    metadata is an optional mapping between the repr of each exp
    (within exps) and the metadata for display.
    If not specified, uses the repr of each exp
    """
    if not metadata:
        d = {repr(exp): exp.df for exp in exps}
    else:
        d = {metadata.get(repr(exp)): exp.df for exp in exps}
    df = pd.concat(d.values(), axis=1, keys=d.keys())
    return df


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
    if not any(isinstance(exp, E2G) for exp in [exp1, exp2]):
        raise TypeError("Incorrect input type")

    for exp in [exp1, exp2]:
        if "iBAQ_dstrAdj_raw" in exp.df.columns:  # revert columns
            exp.df.rename(columns={"iBAQ_dstrAdj": "ibaq_norm"}, inplace=True)
            exp.df.rename(columns={"iBAQ_dstrAdj_raw": "iBAQ_dstrAdj"}, inplace=True)

    if normalize is not None:
        for exp in [exp1, exp2]:
            if normalize.strip() == "mean":
                exp.df["ibaq_norm"] = exp.df.iBAQ_dstrAdj / exp.df.iBAQ_dstrAdj.mean()
            elif normalize.strip() == "median":
                exp.df["ibaq_norm"] = exp.df.iBAQ_dstrAdj / exp.df.iBAQ_dstrAdj.median()
            elif normalize.strip() == "sum":
                exp.df["ibaq_norm"] = exp.df.iBAQ_dstrAdj / exp.df.iBAQ_dstrAdj.sum()
            exp.df.rename(columns={"iBAQ_dstrAdj": "iBAQ_dstrAdj_raw"}, inplace=True)
            exp.df.rename(columns={"ibaq_norm": "iBAQ_dstrAdj"}, inplace=True)

    funcat_cols = ["GeneCapacity", "GeneSymbol", "GeneDescription", "FunCats"]
    funcat1 = exp1.df[funcat_cols]
    funcat2 = exp2.df[funcat_cols]
    funcats = pd.concat([funcat1, funcat2]).drop_duplicates()

    joinexp = E2G()
    for attr in [key for key in joinexp.__dict__.keys() if not key.startswith("_")]:

        value = (
            str(exp1.__getattribute__(attr)) + " vs " + str(exp2.__getattribute__(attr))
        )
        joinexp.__setattr__(attr, value)

    nofuncats = [
        col for col in exp1.df.columns if col not in funcat_cols + ["GeneID.1"]
    ]
    joinexp._df = exp1.df[nofuncats].join(
        exp2.df[nofuncats],
        lsuffix="_x",
        rsuffix="_y",
        how="outer",
    )
    joinexp._df = joinexp._df.join(funcats, how="left")
    joinexp._df["GeneID"] = [
        str(int(x)) if not np.isnan(x) else x for x in joinexp.df.index
    ]  # for convienence
    joinexp._joined = True
    joinexp.ibaq_normalize = normalize
    if _classify == False:
        return joinexp
    try:
        score_experiments(joinexp, seed=seed)
    except Exception as e:
        print("Error scoring experiments")
        raise (e)
        print(e)
    return joinexp


from .data_loading import *
