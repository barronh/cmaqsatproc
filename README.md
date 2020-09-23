# cmaqsatproc

Satellite Processors designed for simple CMAQ comparisons.

## What can you do?

* convert L2 satellite products to L3 on CMAQ grids
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
  * Satellite pixel centers are matched to CMAQ projected grids.
  * No attempt to apply area-fractions is made.
* Satellite AveragingKernels are averaged
  * within a single day
  * within grid cells
* CMAQ stratosphere is removed according to potential vorticity.

## How to?

Edit the scripts below and run them

* ./convertsat.sh
* ./convertcmaq.sh

## Prerequisites

* PseudoNetCDF (http://github.com/barronh/PseudoNetCDF) and associated prereqs

