# cmaqsatproc

Satellite Processors designed for simple CMAQ comparisons.

*Note*: Updating to a version 2 that will look very different, but make working with TropOMI and geostationary data easier. Work is currently in the cmaqsatproc2 branch.

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

## Annotated Configuration Options

The configuration files are Python files with a dictionary definition.
This is functionally like a json file but, unlike json, can be annotated with
comments. Below is an example config.py file with extensive annotation

```
{
    # HDF-EOS5 has variables in groups, so we define a data group
    # and a geographic group, which each has keeys to process.
    "datagrp": "HDFEOS/SWATHS/OMI Vertical Ozone Profile/Data Fields",
    "datakeys": [
        "O3TroposphericColumn", "O3AveragingKernel", "TropopauseIndex",
        "ProfileLevelPressure", "RMS", "ExitStatus", "AverageResiduals",
        "EffectiveCloudFraction"
    ],
    "geogrp": "HDFEOS/SWATHS/OMI Vertical Ozone Profile/Geolocation Fields",
    "geokeys": ["Time", "SolarZenithAngle", "Latitude", "Longitude"],

    # HDF-EOS5 has a structured meta-data, but netCDF4 does not process it
    # instead, all the dimensions have names like "phony_dim_1". To make
    # processing easier, the configuration file has renaming. The dimensions
    # must include nTimes, nXtrack, nLevels. nLevelEdges is optional
    # Dimensions are separately named for data and geo groups
    "datadims": {
        "phony_dim_0": "nTimes",
        "phony_dim_1": "nXtrack",
        "phony_dim_4": "nLevels",
        "phony_dim_8": "nLevelEdges"
    },
    "geodims": {
        "phony_dim_9": "nTimes",
        "phony_dim_10": "nXtrack"
    },
    # Optional, when pressure is available, it should be specified. This allows
    # the vertical coordinate in IOAPI to be calculated
    "pressurekey": "ProfileLevelPressure",

    # Optional, some files (e.g., OMPROFOZ) have the surface level at the top
    # of the file (e.g., (0 hPa, 10 hPa, ..., 1013 hPa). IOAPI expects the order
    # to be reversed, so flipdims identifies files that need to be reordered.
    "flipdims": ["nLevels", "nLevelEdges"],

    # Data filtering is needs to be designed on a source-basis, and is
    # generally composed of a ground-level (i.e., 2D) and a data (sometimes
    # 3D) component. Here, I am showing OMPROFOZ filters
    "grndfilterexpr": "(SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(ExitStatus[:] <= 0) | (ExitStatus[:] >= 10) | " +
        "(RMS[:].max(-1) > 3) | (AverageResiduals[:].max(-1) >= 3) | "
        "(EffectiveCloudFraction >= 0.3)"
    ),

    # IOAPI can only have names that are 16 characters or less.
    # Less is better. To be safe, we use 15 as an upper limit.
    # Each variable also has a count variabel (i.e., N<name>)
    # so base names should be 14 in length.
    "renamevars": {
        "O3TroposphericColumn": "VCD_O3_Trop",
        "O3APrioriProfile": "O3APriori",
        "O3RetrievedProfile": "O3Retrieved",
        "O3AveragingKernel": "AvgKernel",
        "TropopauseIndex": "TropopauseIdx",
        "ProfileLevelPressure": "ProfilePress",
        "AverageResiduals": "AvgResiduals",
        "EffectiveCloudFraction": "EffectCldFrac"
    }
}
```
