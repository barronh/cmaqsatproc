KEY=NO2
YYYYhMMhDD=2016-01-01
YYYYmMMDD=$(date -ud "${YYYYhMMhDD}" +%Ym%m%d)
YYYYMMDD=$(date -ud "${YYYYhMMhDD}" +%Y%m%d)
YYMMDD=$(date -ud "${YYYYhMMhDD}" +%y%m%d)
OUTPATH=CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.${KEY}.nc
AKPATH=OMI-Aura_L2-OMNO2_2016m0101_v003.nc
AKKEY="ScatteringWeight[:]/AmfTrop[:]"
echo ${OUTPATH}
if [ ! -e ${OUTPATH} ]; then
python scripts/cmaq2sat.py \
  /work/ROMO/global/CMAQv5.2.1/2016fe_hemi_cb6_16jh/108km/basecase/extr/O3NO2SO2/CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.O3NO2SO2FORM \
  /work/ROMO/global/CMAQv5.2.1/2016fe_hemi_cb6_16jh/108km/basecase/input/mcip_extr/SAT/METCRO2D.108NHEMI2.44L.${YYMMDD}.SAT \
  /work/ROMO/global/CMAQv5.2.1/2016fe_hemi_cb6_16jh/108km/basecase/input/mcip_extr/SAT/METCRO3D.108NHEMI2.44L.${YYMMDD}.SAT \
  ${OUTPATH} ${KEY} \
  ${AKPATH} ${AKKEY}
fi
