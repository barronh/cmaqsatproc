# cmaqsatproc

Satellite Processors designed for simple CMAQ comparisons.

## What can you do?

* convert L2 or L3 satellite products to L3 on CMAQ grids
* convert CMAQ concentrations to L3

## What makes it simple?

* Converting L2 satellite products takes:
  * raw L2 inputs in HDF-EOS5 (.he5), or
  * opendap urls (requires user opendap configuration), or
  * Common Metadata Repository short_name and daterange
* Converting CMAQ uses raw inputs and outputs.
  * Output: 3D CONC file
  * Inputs:
    * METCRO2D (only CFRAC, PRSFC)
    * METCRO3D (only PV)
* It makes some simplifying assumptions.

## What assumptions are being made?

* Spatial matching is pretty good
  * For satellite products with pixel corners, fractional area weighting is used.
  * For satellite products like MODIS, satellite pixel centers are matched to CMAQ projected grids. No attempt to apply area-fractions is made.
* Satellite AveragingKernels are averaged
  * within a single day
  * within grid cells
* CMAQ stratosphere is removed according to potential vorticity.

## How to?

An example notebook and python script is provided in the examples folder.

* ./examples/basic.ipynb
* ./examples/modis.py
* ./examples/omi_no2.py

## Prerequisites

* numpy
* xarray
* pyproj
* pandas
* geopandas
* PseudoNetCDF (http://github.com/barronh/PseudoNetCDF) and associated prereqs
