"""
Create Custom L3 from TropOMI
=============================

This script is designed to create a custom L3 file from TropOMI NO2, but can
be edited to work with TropOMI HCHO, VIIRS_AERDB, OMNO2, OMHCHO, ...
"""
# %%
# Import Library and Configure
# ----------------------------
import cmaqsatproc as csp

# Using common EPA Lambert Conic Conformal 12-km grid
GDNAM = '12US1'
# Doing just one day
date='2019-07-24'
# Using TropOMINO2 reader
readername = 'TropOMINO2' # or TropOMIHCHO, VIIRS_AERDB, ...
# Defining output path
outpath = f'{readername}_{date}_{GDNAM}.nc'

# %%
# Use NASA's CMR to find and process files
# ----------------------------------------

# Get a CMAQ grid definition
cg = csp.open_griddesc(GDNAM)
# Get the Satellite Reader object
satreader = csp.reader_dict[readername]
# Use CMAQ grid definition and date to drive cmr query
l3 = satreader.cmr_to_level3(
    temporal=f'{date}T00:00:00Z/{date}T23:59:59Z',
    bbox=cg.csp.bbox(), grid=cg.csp.geodf, verbose=9
)
# Save file to disk
l3.to_netcdf(outpath)
