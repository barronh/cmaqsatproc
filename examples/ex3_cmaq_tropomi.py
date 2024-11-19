"""
TropOMI/CMAQ L3
===============

This script is designed to create a CMAQ L3 file including alternative TropOMI
NO2 using the CMAQ AMF.

It assumes that you have a 3D CMAQ CONC file, the METCRO3D file, and a TropOMI
NO2 custom L3 file. The TropOMI custom L3 file is created in the TropOMI L3
OpenDAP example. For the 3D CONC or METCRO3D, there are two options.

For CONC and METCRO3D, you can either use your own data or download data from
the EPA's Air QUAlity TimE Series (`EQUATES`_) Project. The EQUATES project is
available through EPA's Remote Sensing Information Gateway (RSIG) data archive.
RSIG data is served one variable at a time, so you will have four files: (1) 3D
CONC file and (3) METCRO3D files -- DENS, ZF, and PRES. The links below are each
approximately 1 GB in size and must be saved using the names suggested.

* Download `CMAQ CONC3D NO2 <NO2_3D>`_ and save as `CCTM_CONC_2019-07-14_12US1.nc`
* Download `METCRO3D PRES <PRES_3D>`_ and save as `METCRO3D_PRES_2019-07-14_12US1.nc`
* Download `METCRO3D DENS <DENS_3D>`_ and save as `METCRO3D_DENS_2019-07-14_12US1.nc`
* Download `METCRO3D ZF <ZF_3D>`_ and save as `METCRO3D_ZF_2019-07-14_12US1.nc`

.. _EQUATES: https://www.epa.gov/cmaq/equates
.. _NO2_3D: https://ofmpub.epa.gov/rsig/rsigserver?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=netcdf-ioapi&TIME=2019-07-24T00:00:00Z/2019-07-24T23:59:59Z&BBOX=-135.000000,15.000000,-55.000000,70.000000&COVERAGE=cmaq.equates.conus.conc.NO2&COMPRESS=0&NOLONLATS=1
.. _PRES_3D: https://ofmpub.epa.gov/rsig/rsigserver?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=netcdf-ioapi&TIME=2019-07-24T00:00:00Z/2019-07-24T23:59:59Z&BBOX=-135.000000,15.000000,-55.000000,70.000000&COVERAGE=cmaq.equates.conus.metcro3d.PRES&COMPRESS=0&NOLONLATS=1
.. _DENS_3D: https://ofmpub.epa.gov/rsig/rsigserver?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=netcdf-ioapi&TIME=2019-07-24T00:00:00Z/2019-07-24T23:59:59Z&BBOX=-135.000000,15.000000,-55.000000,70.000000&COVERAGE=cmaq.equates.conus.metcro3d.DENS&COMPRESS=0&NOLONLATS=1
.. _ZF_3D: https://ofmpub.epa.gov/rsig/rsigserver?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=netcdf-ioapi&TIME=2019-07-24T00:00:00Z/2019-07-24T23:59:59Z&BBOX=-135.000000,15.000000,-55.000000,70.000000&COVERAGE=cmaq.equates.conus.metcro3d.ZF&COMPRESS=0&NOLONLATS=1
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

# Open a CMAQ 3D CONC file that has the NO2 variable.
qf = csp.open_ioapi(f'CCTM_CONC_{date}_{GDNAM}.nc')[['NO2']]

# Open a single METCRO3D file or separate files if downloaded from RSIG.
# Defaulting to separate files
presf = csp.open_ioapi(f'METCRO3D_PRES_{date}_{GDNAM}.nc')
densf = csp.open_ioapi(f'METCRO3D_DENS_{date}_{GDNAM}.nc')
zff = csp.open_ioapi(f'METCRO3D_ZF_{date}_{GDNAM}.nc')

# Combine DENS, ZF and PRES from METCRO3D files (or file) with CONC
qf['DENS'] = densf['DENS']
qf['ZF'] = zff['ZF']
qf['PRES'] = presf['PRES']

# Create satellite according to CMAQ, and CMAQ according to satellite
overf = satreader.cmaq_process(qf, l3)
overf.to_netcdf(outpath)
