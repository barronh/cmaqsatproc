"""
TropOMI L3 OpenDAP
==================

This script is designed to create a custom L3 file from TropOMI NO2 files
streamed remotely from an OpenDAP server. OpenDAP access to TropOMI data
requires that you:

- Configure `.netrc`, `.urs_cookies`, and `.dodsrc` (see example)
- accept the "NASA GESDISC DATA ARCHIVE" and "Sentinel" End User License
  Agreements (EULAs) on urs.earthdata.nasa.gov using the username and password.

This example can be edited to work with TropOMI HCHO, VIIRS_AERDB, OMI OMNO2,
OMI OMHCHO, ..., etc. The non-TropOMI products do not require the Sentinal
EULA."""
# %%
# Import Library and Configure
# ----------------------------
import cmaqsatproc as csp

# Using common EPA Lambert Conic Conformal 12-km grid
GDNAM = '12US1'
# Doing just one day
date='2019-07-24'
# Using TropOMINO2 reader
readername = 'TropOMINO2' # or TropOMIHCHO, AERDB_L2_VIIRS_SNPP, ...
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
