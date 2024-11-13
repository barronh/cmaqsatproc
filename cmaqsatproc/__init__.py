__all__ = [
    'utils', 'readers', 'cmaq', 'drivers', 'open_ioapi', 'open_griddesc',
    'reader_dict', 'print_reader_list'
]

__doc__ = """
Overview
========

cmaqsatproc provides satellite data processing for CMAQ. This has four
basic steps:
  1. Create a custom level 3 gridded satellite product on CMAQ grid,
  2. Subset CMAQ consistent with satellite overpass and valid samples,
  3. integrate CMAQ mixing ratios to number density (or other metric), and
  4. apply CMAQ air mass factor to the satellite.

Core Objects
============

  * open_griddesc : define a CMAQ grid by name or GRIDDESC.
  * reader_dict : dictionary of satellite readers.
  * open_ioapi : Used to open CMAQ and MCIP data.
  * print_reader_list : useful to find out what can be done.

By TropOMI NO2 Example
======================

    # Import Libraries
    import cmaqsatproc as csp

    # Define grid, satellite data (by reader), and date to process
    GDNAM = '12US1'
    date='2019-07-24'
    readername = 'TropOMINO2' # or TropOMIHCHO, VIIRS_AERDB, ...

    cg = csp.open_griddesc(GDNAM)
    satreader = csp.reader_dict[readername]

    # Custom L3 output path
    outl3path = f'{readername}_{date}_{GDNAM}.nc'
    # CMAQ satellite output
    outcsppath = f'{readername}_{date}_{GDNAM}_CMAQ.nc'

    # Use NASA Common Metadata Repository to find satellite data
    # and remote access methods. Grid results on CMAQ grid.
    l3 = satreader.cmr_to_level3(
        temporal=f'{date}T00:00:00Z/{date}T23:59:59Z',
        bbox=cg.csp.bbox(), grid=cg.csp.geodf, verbose=9
    )
    # Optionally, store custom l3 result
    l3.to_netcdf(outl3path)

    # Collect CMAQ 3D CONC and MCIP data
    qf = csp.open_ioapi(f'CCTM_CONC_{date}_{GDNAM}.nc')[['NO2']]
    mf = csp.open_ioapi(f'METCRO3D_{date}_{GDNAM}.nc')
    qf['DENS'] = mf['DENS']
    qf['ZF'] = mf['ZF']
    qf['PRES'] = mf['PRES']

    # cmaq_process applies temporal subsetting consistent with the satellite
    # overpass, integrates NO2 to mole density, and applies CMAQ vertical
    # to create an alternative AMF and "adjusted" satellite column
    overf = satreader.cmaq_process(qf, l3)
    overf.to_netcdf(outcsppath)

"""

from . import utils
from . import readers
from . import cmaq
from . import drivers

__version__ = '0.4.0'

reader_dict = readers.reader_dict
open_ioapi = cmaq.open_ioapi
open_griddesc = cmaq.open_griddesc
print_reader_list = readers.print_reader_list
