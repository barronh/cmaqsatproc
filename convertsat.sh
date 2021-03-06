#!/bin/bash
mkdir -p outsat

# Make OMNO2 L3 output from downloaded files
OUTPATH=outsat/OMI-Aura_L3_108NHEMI2-OMNO2_2016m0101_v003.nc
CFGPATH=config/OMNO2v003.py
INPATHS=/work/ROMO/users/bhenders/obs/OMNO2d/aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMNO2.003/2016/001/*
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

# Make OMNO2 L3 output from prequeried OpenDAP paths
# Requires .dodsrc or .daprc configured to allow
# NASA OpenDAQ access
OUTPATH=outsat/OMI-Aura_L3_108NHEMI2-OMNO2_dap_2016m0101_v003.nc
INPATHS="https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2016/001/OMI-Aura_L2-OMNO2_2016m0101t1410-o60978_v003-2019m0819t154317.he5 https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2016/001/OMI-Aura_L2-OMNO2_2016m0101t1549-o60979_v003-2019m0819t154325.he5 https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2016/001/OMI-Aura_L2-OMNO2_2016m0101t1727-o60980_v003-2019m0819t154324.he5 https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2016/001/OMI-Aura_L2-OMNO2_2016m0101t1906-o60981_v003-2019m0819t154319.he5 https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2016/001/OMI-Aura_L2-OMNO2_2016m0101t2045-o60982_v003-2019m0819t154319.he5"
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

# Make OMNO2 L3 output from CMR query for OpenDAP paths
# Requires .dodsrc or .daprc configured to allow
# NASA OpenDAQ access
OUTPATH=outsat/OMI-Aura_L3_108NHEMI2-OMNO2_dapq_2016m0101_v003.nc
INPATHS="{'short_name': 'OMNO2', 'daterange': '2016-01-01T00:00:00Z/P01D'}"
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

# Make OMPROFOZ L3 output from downloaded files
# I am not aware of a OMPROFOZ OpenDAP path
OUTPATH=outsat/OMI-Aura_L3_108NHEMI2-OMPROFOZ_2016m0101_v003.nc
CFGPATH=config/OMPROFOZv003.py
INPATHS=/work/ROMO/users/bhenders/obs/OMPROFOZ/avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2/OMPROFOZ/2016/01/01/*.he5
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

CFGPATH=config/MLO3v005he5.py
INPATHS=/work/ROMO/satellite/ML2O3/acdisc.gesdisc.eosdis.nasa.gov/data/Aura_MLS_Level2/ML2O3.005/2020/MLS-Aura_L2GP-O3_v05-01-c01_2020d267.he5
OUTPATH=outsat/MLS-Aura_L3_108NHEMI2-O3_v05-01-c01_2020d267.nc
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

CFGPATH=config/MLO3v005opendap.py
INPATHS="{'short_name': 'ML2O3', 'version': '005', 'daterange': '2020-01-01T00:00:00Z/P01D'}"
OUTPATH=outsat/MLS-Aura_L3_108NHEMI2-O3_dapq_v05-01-c01_2020d001.nc
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}

# Make OMNHCHO L3 output from downloaded files
OUTPATH=outsat/OMI-Aura_L3_108NHEMI2-OMHCHO_2016m0101_v003.nc
CFGPATH=config/OMHCHOv003.py
INPATHS=/work/ROMO/satellite/OMHCHO/aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMHCHO.003/2016/001/*he5
python scripts/sat2cmaq.py GRIDDESC 108NHEMI2 ${OUTPATH} ${CFGPATH} ${INPATHS}
