mkdir -p figs
python scripts/plot_compare.py \
  outsat/OMI-Aura_L2-OMNO2_2016m0101_v003.nc "ColumnAmountNO2[:]/2.69e16; VCD.units = 'dobsons'" \
  outcmaq/CCTM_CONC_v521_intel17.0_HEMIS_cb6_20160101.NO2.nc NO2 \
  config/plotopt.py figs/NO2.png
