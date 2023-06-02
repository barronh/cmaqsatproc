2023-02-07:
* Only updates to open_griddesc and ioapi functionality
* Updated 1US1, 4US2, 1US2 to use correct cell sizes in open_griddesc
* Added robustness for IOAPI dates. SDATE -365 is converted to 1970001
* Added a to_ioapi cleanup for output.

2022-11-21:
* Updated TropOMI reader for compatibility with newer versions of xarray.

2022-11-03:
* Updated TropOMI CO, NO2, HCHO, and CH4 to automatically use _HiR products for dates on or after 2019-08-06. Before, you had to manually set the short_name. This is an ease of use update that I am treating as a bug, so the version incremented to 0.2.2.

2022-10-25:
* Updated TropOMI CO to use meter-based averaging kernel.
* Updated to_level3 and wrappers to add cmaqsatproc_version to outputs
* Updated version to 0.2.1

2022-10-14:
* Updated qa minimum values for GOES and TropOMI
* Added GOES AOD, OMPS, and IASI

2022-09-19:
* VIIRS Dark Target and Deep Blue functionality added.
* OpenDAP is currently not working for VIIRS. Links exist for .nc.html, but provides error `OSError: [Errno -70] NetCDF: DAP server error: b'https://ladsweb.modaps.eosdis.nasa.gov/opendap/allData/5110/AERDB_L2_VIIRS_SNPP/2019/262/AERDB_L2_VIIRS_SNPP.A2019262.0948.001.2019263132953.nc'`
