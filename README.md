# cmaqsatproc

Satellite Processors designed for simple CMAQ comparisons.

## What can you do?

* convert L2 satellite products to L3 on CMAQ grids
* convert CMAQ concentrations to L3

## What make it simple?

* It takes raw CMAQ inputs and outputs
* It makes some simplifying assumptions.

## What assumptions are being made?

* Satellite pixel centers are matched when they are within CMAQ projected
  grid cells. No attempt to apply area-fractions is made.
* CMAQ stratosphere is removed according to 

## How to?

* ./convertsat.sh
* ./convertcmaq.sh

