#!/bin/bash
mkdir -p outsat
python scripts/sat2cmaq.py \
  GRIDDESC \
  HEMIS \
  outsat/OMI-Aura_L2-OMNO2_2016m0101_v003.nc \
  config/OMNO2v003.py \
  /work/ROMO/users/bhenders/obs/OMNO2d/aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMNO2.003/2016/001/*
