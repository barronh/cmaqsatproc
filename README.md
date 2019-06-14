# cmaqsatproc

Satellite Processors designed for simple CMAQ comparisons.

## What can you do?

* convert L2 satellite products to L3 on CMAQ grids
* convert CMAQ concentrations to L3

## What make it simple?

* Converting L2 satellite products takes raw L2 inputs
* Converting CMAQ uses raw inputs and outputs.
  * Output: 3D CONC file
  * Inputs:
    * METCRO2D (CFRAC, PRSFC)
    * METCRO3D (PV)
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

* ./convertsat.sh
* ./convertcmaq.sh

## Prerequisites

* PseudoNetCDF (http://github.com/barronh/PseudoNetCDF) and associated prereqs

