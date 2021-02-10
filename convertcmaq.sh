KEY=NO2
AKPATH=outsat/OMI-Aura_L3_108NHEMI2-OMNO2_2016m0101_v003.nc
AKKEY="ScatWgt[:]/AmfTrop[:, 0]"
# KEY=O3
# AKPATH=outsat/OMI-Aura_L3_108NHEMI2-OMPROFOZ_2016m0101_v003.nc
# AKKEY="O3AveragingKernel[:]"
YYYYhMMhDD=2016-01-01
YYYYmMMDD=$(date -ud "${YYYYhMMhDD}" +%Ym%m%d)
YYYYMMDD=$(date -ud "${YYYYhMMhDD}" +%Y%m%d)
YYMMDD=$(date -ud "${YYYYhMMhDD}" +%y%m%d)
OUTPATH=outcmaq/CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.${KEY}.nc
DATAROOT=/work/ROMO/global/CMAQv5.2.1/2016fe_hemi_cb6_16jh/108km/basecase/

echo ${OUTPATH}
if [ ! -e ${OUTPATH} ]; then
  mkdir -p $(dirname ${OUTPATH})

  CONCPATH=${DATAROOT}/extr/O3NO2SO2/CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.O3NO2SO2FORM
  M2DPATH=${DATAROOT}/input/mcip_extr/SAT/METCRO2D.108NHEMI2.44L.${YYMMDD}.SAT
  M3DPATH=${DATAROOT}/input/mcip_extr/SAT/METCRO3D.108NHEMI2.44L.${YYMMDD}.SAT
  M3DPATH=${DATAROOT}/../ZEAS/input/mcip/METCRO3D.108NHEMI2.44L.${YYMMDD}

python scripts/cmaq2sat.py \
  -v --satpath ${AKPATH} --akexpr "${AKKEY}" \
  --tppexpr "TropopausePres[:] * 100" \
  ${CONCPATH} ${M2DPATH} ${M3DPATH} \
  ${OUTPATH} ${KEY}

fi
