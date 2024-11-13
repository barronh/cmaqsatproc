"""
TropOMI/CMAQ L3
===============

This script is designed to create a CMAQ L3 file including alternative TropOMI
NO2 using the CMAQ AMF.

It assumes that you have a 3D CMAQ CONC file, the METCRO3D file, and a TropOMI
NO2 custom L3 file (see tropomi_cmr.py)
"""
# %%
# Import Library and Configure
# ----------------------------
import cmaqsatproc as csp
import xarray as xr

# Using common EPA Lambert Conic Conformal 12-km grid
GDNAM = '12US1'
# Doing just one day
date='2019-07-24'
# Using TropOMINO2 reader
readername = 'TropOMINO2'
# Defining output path
outpath = f'{readername}_{date}_{GDNAM}_CMAQ.nc'

# %%
# Gather Data and Processors
# --------------------------

# Get satellite reader
satreader = csp.reader_dict[readername]
# Get custom L3 file made earlier
l3 = xr.open_dataset(f'{readername}_{date}_{GDNAM}.nc')

# Open CMAQ CONC and METCRO3D files and combine met variables with CONC
qf = csp.open_ioapi(f'CCTM_CONC_{date}_{GDNAM}.nc')[['NO2']]
mf = csp.open_ioapi(f'METCRO3D_{date}_{GDNAM}.nc')
qf['DENS'] = mf['DENS']
qf['ZF'] = mf['ZF']
qf['PRES'] = mf['PRES']

# Create satellite according to CMAQ, and CMAQ according to satellite
overf = satreader.cmaq_process(qf, l3)
overf.to_netcdf(outpath)
