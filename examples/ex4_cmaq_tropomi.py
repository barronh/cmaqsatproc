"""
Create Custom L3 from VIIRS
===========================

This script is designed to make custom L3 files from SNPP VIIRS Deep Blue
using files that have already been downloaded. It can be edited to work with
TropOMI HCHO, VIIRS_AERDB, OMNO2, OMHCHO, ...

This example assumes you have downloaded satellite files. The code is largely
the same as the previous. Instead of `cmr_to_level3`, it the method uses `glob`
to make a list of files that it passes to `paths_to_level3`.
"""
# %%
# Import Library and Configure
# ----------------------------
from glob import glob
import cmaqsatproc as csp

GDNAM = '12US1'
date='2019-07-24'
readername = 'VIIRS_AERDB' # or TropOMIHCHO, VIIRS_AERDB, ...
outpath = f'{readername}_{date}_{GDNAM}.nc'

cg = csp.open_griddesc(GDNAM)
satreader = csp.reader_dict[readername]

paths = sorted(glob('AERDB_L2_VIIRS_SNPP*.nc'))
l3 = satreader.paths_to_level3(
    paths, bbox=cg.csp.bbox(), grid=cg.csp.geodf, verbose=9
)
l3.to_netcdf(outpath)