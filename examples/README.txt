cmaqsatproc Example Gallery
===========================

This gallery houses examples on how to use cmaqsatproc with various satellite
products and using different data acquisition methods.

Satellite Products
==================

* Configure
  * config_netrc.py : setup to use OpenDAP connection.
* TropOMI:
  * tropomi_cmr.py: create custom L3 file from TropOMI via OpenDAP server.
  * cmaq_tropomi.py: create custom L3 file from CMAQ and TropOMI with CMAQ AMF
* VIIRS
  * viirs_local.py: create custom L3 file from VIIRS via local files.


All examples require `cmaqsatproc` which can be installed via `shell` or in a
Jupyter Notebook. In `shell`, the command is:

.. code-block:: bash

   python -m pip install --user cmaqsatproc

In notebooks (e.g., on Google Colab), this can be done live with
the command below (may require kernel restart):

.. code-block:: python

   %pip install --user cmaqsatproc

To run an example in a notebook:

1. Copy the command above into a new cell and run it.
2. Click on example below to see code.
3. Copy the code from the example into a new cell and run it.

