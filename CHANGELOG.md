2026-02-25 v0.5.1:
* Added support for Python3.12 https://github.com/pandas-dev/pandas/issues/57500

2026-02-15 v0.5.0:
* Added cmr_download method to retrieve files based on cmr keywords.
* Added TEMPO support for NO2 and HCHO
* Added CAMx support.
* Added filterfunc named keyword support

2024-11-19 v0.4.1:
* Dramatically improved documentation.
* Fixed CMAQ date issue when just one time is available.

2024-11-10 v0.4.0:
* Use xarray.Dataset.sizes instead of dims to avoid warning.
* Remove pygeos dependency since geopandas is moving to shapely 2.0
* Add CMAQ SDATE as date for time independent files.
* TropOMI extending warning and date ranges (issue #9)

2024-01-25 v0.3.0:
* Added methane_mixing_ratio_bias_corrected to the default key for S5P_L2__CH4___
* Changed all cmr_links commands to use concept_id instead of short_name.
* Altered grouped_weighted_avg to output weight_mean and weight_sum instead of weighted weight.
* Updated open_griddesc to read GRIDDESC with commas, double precision (D or d), and comments
* Updated documentation

2023-10-12 v0.2.5:
* Minor updates to cmaq.py for date robustness
* Updated OMI readers to correctly handle CloudFraction with scaling metadata in local he5 files
  * OpenDAP and NetCDF files correctly scale cloud fractions to (0, 1)
  * he5 has ScaleFactor and Offset metadata, but is not being scaled correctly.
  * cloudleq has been updated to check

2023-02-07 v0.2.4:
* Only updates to open_griddesc and ioapi functionality
* Updated 1US1, 4US2, 1US2 to use correct cell sizes in open_griddesc
* Added robustness for IOAPI dates. SDATE -365 is converted to 1970001
* Added a to_ioapi cleanup for output.

2022-11-21 v0.2.3:
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
