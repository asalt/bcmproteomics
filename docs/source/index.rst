.. bcmproteomics documentation master file, created by
   sphinx-quickstart on Sun Aug 28 19:34:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to bcmproteomics's documentation!
=========================================

.. warning:: Make sure that both bcmproteomics and bcmproteomics_ext is installed.
             bcmproteomics requires additional drivers for connection to database that
             are not available on all systems.
             You will need to know the IP address for data access via bcmproteomics_ext [#f1]_


Quickstart
##########
Load some data:
  >>> from bcmproteomics import ispec
  >>> experiment = ispec.E2G(recno=12345, runno=1, searchno=1)
  >>> experiment.df
  pandas.DataFrame with available values.

A few examples of convenient data access:
  >>> experiment.added_by
  alex
  >>> experiment.tfs
  subselection of gene products marked as transcription factors
  >>> experiment.strict
  subselection of gene products identified with high confidence


Save loaded data:
  >>> experiment.save('./data/directory/')

Automatically save and load data from a directory::
  >>> from bcmproteomics import ispec
  >>> experiment = ispec.E2G(recno=12345, runno=1, searchno=1, data_dir='./data/directory')

.. note::
   Whenever ``data_dir`` is specified, the experiment will first attempt to be loaded from the given directory.
   If the data is not present, then the data is accessed over the network and automatically saved.
   To reload data over the network after it has been saved locally, simply remove the local files.

.. rubric:: Footnotes
.. [#f1] Current IP address is 10.16.3.226:5000 but is subject to change.

Contents:
##########

.. toctree::
   :maxdepth: 2
.. automodule:: bcmproteomics.ispec
    :members:
.. automodule:: bcmproteomics.multicomparison
    :members:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
