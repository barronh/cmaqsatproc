.. cmaqsatproc documentation master file, created by
   sphinx-quickstart on Fri Mar 17 21:46:43 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

User's Guide
============

Satellite Processors designed for simple CMAQ comparisons.

What can you do?
----------------

* convert L2 or L3 satellite products to L3 on CMAQ grids
  * 2-d species like total or tropospheric columns
  * n-d vairables like averaging kernels and scattering weights.
* convert CMAQ concentrations to L3-like products
  * Apply satellite averaging kernels to CMAQ concentrations to make satellite-like CMAQ
  * Apply CMAQ to create alternative air mass factors to make CMAQ-like satellite products.

What makes it simple?
---------------------

`cmaqsatproc` has an easy full suite approach

1. Operates on local files or dynamically finds remote files
  1. User can specify input files from their disk.
  2. Queries NASA's Common Metadata Repository (CMR) or NOAA AWS
2. Allows for spatial subsetting based on a simple box.
  1. User can specify the box based on lat/lon
  2. The CMAQ grid can be used to automatically define the box.
3. Provides L2 access as a dataframe or makes Level 3 data as a dataset
4. `Simple instructions <https://github.com/barronh/cmaqsatproc/blob/main/COLABINSTALL.md>`_ are provided to configure Google Colab.


Getting Started
---------------

The best way to get started is to install (see below) and then explore the
examples gallery.


Installation
------------

`cmaqsatproc` is avalable through pypi.org, but is still in rapid development. You
can get the latest release from pypi via the command below.

.. code-block::

    pip install cmaqsatproc

Or, you can get the latest version with this command.

.. code-block::

    pip install git+https://github.com/barronh/cmaqsatproc.git


Issues
------

If you're having any problems, open an issue on github.

https://github.com/barronh/cmaqsatproc/issues


Quick Links
-----------

* :doc:`auto_examples/index`
* :doc:`cmaqsatproc`

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Table of Contents

   self
   auto_examples/index
   cmaqsatproc
