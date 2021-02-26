#!/bin/bash

function setpaths(){
  YYYYmMMDD=$(date -ud "${YYYYhMMhDD}" +%Ym%m%d)
  YYYYMMDD=$(date -ud "${YYYYhMMhDD}" +%Y%m%d)
  YYMMDD=$(date -ud "${YYYYhMMhDD}" +%y%m%d)
  CONCPATH=${DATAROOT}/extr/O3NO2SO2/CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.O3NO2SO2FORM
  M2DPATH=${DATAROOT}/input/mcip_extr/SAT/METCRO2D.108NHEMI2.44L.${YYMMDD}.SAT
  M3DPATH=${DATAROOT}/input/mcip_extr/SAT/METCRO3D.108NHEMI2.44L.${YYMMDD}.SAT
  M3DPATH=${DATAROOT}/../ZEAS/input/mcip/METCRO3D.108NHEMI2.44L.${YYMMDD}
  OUTPATH=outcmaq/CCTM_CONC_v521_intel17.0_HEMIS_cb6_${YYYYMMDD}.${KEY}.nc
  mkdir -p outcmaq
}


DATAROOT=/work/ROMO/global/CMAQv5.2.1/2016fe_hemi_cb6_16jh/108km/basecase/
DATES="2016-01-01"

KEY=NO2
for YYYYhMMhDD in ${DATES}
do
setpaths
python scripts/cmaq2sat.py \
  -v -v -v -v \
  --satpath outsat/OMI-Aura_L3_108NHEMI2-OMNO2_${YYYYmMMDD}_v003.nc \
  --akexpr "ScatWgt[:]/AmfTrop[:, 0]" --swexpr "ScatWgt[:]" \
  --tppexpr "TropopausePres[:] * 100" \
  ${CONCPATH} ${M2DPATH} ${M3DPATH} \
  ${OUTPATH} ${KEY}
done

KEY=FORM
for YYYYhMMhDD in ${DATES}
do
setpaths
python scripts/cmaq2sat.py \
  -v -v -v -v \
  --satpath outsat/OMI-Aura_L3_108NHEMI2-OMHCHO_${YYYYmMMDD}_v003.nc \
  --akexpr "ScatWgt[:]/AirMassFactor[:, 0]" --swexpr "ScatWgt[:]" \
  ${CONCPATH} ${M2DPATH} ${M3DPATH} \
  ${OUTPATH} ${KEY}
done

KEY=O3
for YYYYhMMhDD in ${DATES}
do
# Not compatible right now
#  --akexpr "AvgKernel[:]" --prexpr "O3APriori[:]" \
setpaths
python scripts/cmaq2sat.py \
  -v -v -v -v \
  --satpath outsat/OMI-Aura_L3_108NHEMI2-OMPROFOZ_2016m0101_v003.nc \
  --isvalidexpr "(VCD_O3_Trop[:].mask == False)" \
  --maxcld 1 \
  ${CONCPATH} ${M2DPATH} ${M3DPATH} \
  ${OUTPATH} ${KEY}
done
