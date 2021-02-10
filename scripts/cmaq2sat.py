import os
import sys
import time
import numpy as np
import PseudoNetCDF as pnc
import argparse
from warnings import warn

parser = argparse.ArgumentParser()
parser.description = """
Convert 3-D IOAPI mixing ratios (ppm) to dobson units within the
clear tropopsphere at overpass times.

* The troposphere is diagnosed either by:
  1. Satellite Tropopause Pressure.
  2. Potential vorticity < 2 (see Itahashi 2018)
  3. Lapse rate (See WMA 1957 as in doi: 10.1029/2003GL018240

* Clear is defined as cloud fraction (METCRO2D CFRAC) less than
  maxcld (default 0.2)

* Overpass times are defined in local hour between minlsthour and
  maxlsthour local hour is defined as CMAQ UTC + lon / 15.

If provided, the averaging kernel is applied to the CMAQ profile.

If either the averaging kernel or the tropopause pressure are used
from the satellite file, then the output will be restricted to
valid satellite retrievals.
"""

parser.epilog = """
Example No Averaging Kernel:
    python cmaq2sat.py CMAQ.ACONC_20160101.nc \
           MECRO2D_160101 METCRO3D_160101 CMAQ_L3_20160101.nc O3

Example No Averaging Kernel:
    python cmaq2sat.py CMAQ.ACONC_20160101.nc \
           MECRO2D_160101 METCRO3D_160101 CMAQ_L3_20160101.nc O3 AVGK.nc

Where AVGK.nc was made by gridding the AveragingKernel from an L2 file and
averaging to 1-day
"""
aa = parser.add_argument
aa(
    '--minlsthour', default=12, type=float,
    help=(
        'Used with --maxlsthour to restrict times to overpass time. ' +
        'For example, Aura Overpass is 13:40 solar time at the equator ' +
        'and 12:00 solar time at 70 N (minlsthour <= lsthour <= maxlsthour).' +
        'See Figure 2.4 ' +
        'https://aura.gsfc.nasa.gov/images/aura_validation_v1.0.pdf'
    )
)
aa(
    '--maxlsthour', default=14, type=float,
    help='See --minlsthour'
)
aa(
    '--maxcld', default=0.2, type=float,
    help='CFRAC < maxcld'
)
aa(
    '--satpath', default=None,
    help='Averaging Kernel path; contains akexpr and or tppexpr'
)
aa(
    '--tppexpr', default=None, type=str,
    help='expr for TropoPausePressure (Pa) of tropopause in satellite product'
)
aa(
    '--isvalidexpr', default=None, type=str,
    help=(
        'Is valid expr from the satellite product'
        '(e.g., VCD_NO2_Trop.mask == False)'
    )
)
aa('--akexpr', default=None, help='Averaging Kernel expression')
aa('cpath', help='CMAQ CONC path; must contain key')
aa('m2path', help='METCRO2D path; contains CFRAC and PRSFC')
aa('m3path', help='METCRO3D path; contains PV')
aa('outpath', help='output path')
aa('key', help='CMAQ key to process')


def findtroposphere(m3f, rawsatf, myargs):
    tppexpr = myargs['tppexpr']
    if tppexpr is not None and 'PRES' in m3f.variables:
        myargs['istropmethod'] = 'Satellite Tropopause Pressure'
        tppf = rawsatf.eval(
            f'TropoPausePressure = {tppexpr}'
        ).slice(LAY=0)
        # Satellite TropoPause Pressure surface values only
        stpp = tppf.variables['TropoPausePressure'][:]
        istrop = m3f.variables['PRES'][:] > stpp
    elif 'PV' in m3f.variables:
        myargs['istropmethod'] = 'Potential Vorticity < 2 (Itahashi 2018)'
        # Any PV < 2 starts the troposphere. Itahashi 2018 submitted
        PV = m3f.variables['PV'][:]
        istrop = np.cumsum(PV[:, ::-1] < 2, axis=1)[:, ::-1] > 0
    elif 'TA' in m3f.variables and 'ZH' in m3f.variables:
        warn(
            'PV was not available. Using approximation of WMO 1957 ' +
            'lapse-rate definition. see ' +
            'https://dx.doi.org/10.1029/2003GL018240'
        )
        myargs['istropmethod'] = (
            'Lapse-rate (similar to doi: 10.1029/2003GL018240)'
        )
        # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GL018240
        # According to the WMO [1957], the tropopause is defined as "the lowest
        # level at which the lapse-rate decreases to 2 C/km or less, provided
        # that the average lapse-rate between this level and all higher levels
        # within 2 km does not exceed 2 C/km"
        ZH = m3f.variables['ZH'][:] * 1
        dT = m3f.variables['TA'][:] * 1
        dZ = ZH * 1
        # dZ = ZH(l) - ZH(l-1) [=] km
        dZ[:, 1:-1] = (ZH[:, 2:] - ZH[:, :-2]) / 1000
        # dT = TA(l) - TA(l-1) [=] K or C
        dT[:, 1:-1] = (dT[:, 2:] - dT[:, :-2])
        # Set lowest values to layer above them
        dT[:, 0] = dT[:, 1]
        dT[:, -1] = dT[:, -2]
        dZ[:, 0] = dZ[:, 1]
        dZ[:, -1] = dZ[:, -2]
        dTdZ = - dT / dZ
        # assuming upper trop layers are ~1km, using next layer
        # as a surrogate for higher levels within 2 km.
        dTdZ2 = dTdZ * 1
        dTdZ2[:, :-1] = - (dT[:, :-1] + dT[:, 1:]) / (dZ[:, :-1] + dZ[:, 1:])

        # Adding minimum 5km tropopause by not allowing a below min flag
        # below 5km.
        belowmin = np.ma.masked_where(
            ZH[:] < 5000,
            np.ma.logical_and(
                dTdZ[:] < 2,
                dTdZ2[:] < 2
            )
        ).filled(False)
        istrop = np.cumsum(belowmin, axis=1) == 0
    else:
        raise KeyError('METCRO3D requires PV or TA and ZH, but had neither')

    return istrop


def pinterp(m2fd, m3f, aksatf):
    prsfc = m2fd.variables['PRSFC'][:]
    # The output shape will be 3d matching the metcro3d file
    akshape = list(prsfc.shape)
    akshape[1] = len(m3f.dimensions['LAY'])

    ak = np.ma.masked_all(akshape, dtype='f')

    # The partial pressure represented by the satellite
    akptop = aksatf.VGTOP
    akdp = (aksatf.VGLVLS[0] - akptop)

    # Calculate mid point pressures
    akpmid = (aksatf.VGLVLS[:-1] + aksatf.VGLVLS[1:]) / 2

    # Create a sigma approximation
    aksigma = (akpmid - akptop) / akdp

    akin = aksatf.variables['AveragingKernel']

    # Out sigma is mid points
    outsigma = (m3f.VGLVLS[:-1] + m3f.VGLVLS[1:]) / 2
    outvgtop = m3f.VGTOP

    # Iterate over indices (k=0, while t, j, and i change
    print('Start inteprolate AK', flush=True)
    for t, k, j, i in np.ndindex(*prsfc.shape):
        # Hold a vertical column
        psfc = prsfc[t, k, j, i]
        if np.ma.getmaskarray(psfc).all():
            # Skip any columns that are 100% masked
            continue

        akv = akin[t, :, j, i]
        if np.ma.getmaskarray(akv).all():
            # Skip any columns that are 100% masked
            continue
        x = outsigma * (psfc - outvgtop) + outvgtop
        xp = akpmid
        if not pressure:
            x = (x - akptop) / akdp
            xp = aksigma

        ak[t, :, j, i] = np.interp(x, xp[::-1], akv[::-1])

    return ak


def ppm2du(
    cpath, m2path, m3path, outpath, key='O3',
    satpath=None, akexpr=None, tppexpr=None, minlsthour=12, maxlsthour=14,
    isvalidexpr=None, maxcld=0.2
):
    myargs = locals().copy()
    myargs['script'] = sys.argv[0]
    myargs['script_mtime'] = time.ctime(os.stat(sys.argv[0]).st_mtime)

    # outpath = os.path.join('du', os.path.basename(cpath))
    if os.path.exists(outpath):
        print('Using cached:', outpath)
        return

    # Check for weird options
    # Satellite expressions with no satellite path
    # Or satellite path with no satellite options
    if (
        (isvalidexpr is None and akexpr is None and tppexpr is None) and
        (satpath is not None)
    ):
        warn(
            'Satellite file unused; specify akexpr or tropoexpr or isvalidexpr'
        )
    elif (
        (akexpr is not None or tppexpr is not None) and
        (satpath is None)
    ):
        raise ValueError('akexpr or tppexpr provided without satellite')

    inf = pnc.pncopen(cpath, format='ioapi')
    m2f = pnc.pncopen(m2path, format='ioapi')
    m3f = pnc.pncopen(m3path, format='ioapi')
    cmaqunit = getattr(inf.variables[key], 'units', 'unknown').strip().lower()
    if cmaqunit in ('ppm', 'ppmv'):
        unitfactor = 1
    elif cmaqunit in ('ppb', 'ppbv'):
        unitfactor = 1e-3
    else:
        warn(f'Unit {cmaqunit} is unknown; assuming ppm')
        unitfactor = 1

    if satpath is not None:
        rawsatf = pnc.pncopen(satpath, format='ioapi')
    else:
        rawsatf = None
    j = np.arange(inf.NROWS)
    i = np.arange(inf.NCOLS)
    I, J = np.meshgrid(i, j)
    lon, lat = inf.ij2ll(I, J)

    # Define UTC hour using TFLAG where available
    if 'TFLAG' in inf.variables:
        utch = inf.variables['TFLAG'][:, 0, 1][:, None, None, None] / 10000
    else:
        if 'TSTEP' in inf.dimensions:
            nh = len(inf.dimensions['TSTEP'])
        else:
            nh = inf.variables[key].shape[0]
        utch = np.arange(0, nh)[:, None, None, None]

    # Define local soalr hour using 15 degree approximation
    # of time zone offset as an integer to get lst
    tzoff = (lon[None, None, :, :] // 15).round(0).astype('i')
    lsth = utch + tzoff

    # Is overpass based on local lsthour
    isoverpass = (lsth >= minlsthour) & (lsth <= maxlsthour)
    # Cloud free is based on CFRAC and a maximum value
    cf = m2f.variables['CFRAC'][:]
    iscldfree = cf < maxcld

    # Find layers in the troposphere
    istrop = findtroposphere(m3f, rawsatf, myargs)
    # Create a 2D and 3D mask for use later
    is2valid = (isoverpass & iscldfree)

    # Only select data with valid satellite retrieved overpasses
    if isvalidexpr is not None:
        myargs['isvalidmethod'] = f'Paired with satellite {isvalidexpr}'
        isvalidf = rawsatf.eval(
            f'ISVALID = {isvalidexpr}'
        ).slice(LAY=0)
        issatvalid = isvalidf.variables['ISVALID'][:]
        is2valid = is2valid & issatvalid
    else:
        myargs['isvalidmethod'] = 'Model filter only'

    is3valid = is2valid.repeat(inf.NLAYS, 1) & istrop
    is3strat = is2valid.repeat(inf.NLAYS, 1) & (~istrop)

    # subset variables to increase speed
    # create copies of m2f and inf with invalid pixels masked
    infm = inf.subsetVariables(
        [key]
    ).mask(
        ~is3valid, dims=('TSTEP', 'LAY', 'ROW', 'COL')
    )
    m2fm = m2f.subsetVariables(
        ['CFRAC', 'PRSFC']
    ).mask(
        ~is2valid, dims=('TSTEP', 'LAY', 'ROW', 'COL')
    )

    # Average only valid pixels for a daily file
    infd = infm.apply(TSTEP='mean')
    m2fd = m2fm.apply(TSTEP='mean')

    # calculate edge pressures from surface, sigma, and top pressure
    pedges = (
        (m2fd.variables['PRSFC'][:] - infd.VGTOP) *
        infd.VGLVLS[None, :, None, None] + infd.VGTOP
    )
    # https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMO3PR.003/doc/README.OMO3PR.pdf
    # http://www.temis.nl/data/conversions.pdf
    # assumes mixing ratio in PPM and dp in hPa
    hPa_to_du = (
        10 * 1.3807e-23 * 6.022e23 / 0.02894 * 273.15 / 9.80665 / 101325.
    )

    dp = -np.diff(pedges, axis=1) / 100
    ppm = infd.variables[key][:] * unitfactor

    # Using O3 gradient as tropopause
    # commented out to use PV instead
    # dppm = np.concatenate([ppm[:,[0]] * 0, np.diff(ppm, axis=1)], axis=1)
    # mask = (dppm > 0.050).astype('i')
    # cmask = np.cumsum(mask, axis=1) > 0
    # ppm = np.ma.masked_where(cmask, ppm)
    if akexpr is not None:
        # Calculate the averaging Kernel using variables in the raw satellite
        # file
        aksatf = rawsatf.eval(
            f'AveragingKernel = {akexpr}'
        )
        # Now interpolate the Averaging Kernel to the model levels
        # Currently, this supports
        if (
            aksatf.VGTYP == infd.VGTYP
            and np.allclose(aksatf.VGTYP, infd.VGLVLS)
        ):
            # The vertical coordinate types are the same and the vertical
            # level edges are the same, so no interpolation is necessary.
            akf = aksatf
        elif aksatf.VGTYP == infd.VGTYP and aksatf.VGTYP == 7:
            # The vertical coordinates are dry sigma pressures, which support
            # the use of the efficient interpSigma method for vertical
            # interpolation
            myargs['ak_interp'] = 'sigma2sigma'
            akf = aksatf.interpSigma(
                infd.VGLVLS, vgtop=infd.VGTOP, interptype='linear',
                extrapolate=False, fill_value='extrapolate', verbose=0
            )
            # Ensure that missing values are treated as missing
            ak = np.ma.masked_values(akf.variables['AveragingKernel'][:], -999)
        elif aksatf.VGTYP == 4:
            # This case is designed to support OMNO2d-like pressure coordinates
            # In this case, the VGLVLS will be in decreasing value order in
            # Pascals (e.g., 102500.f, 101500.f, ... , 115.f, 45.f)
            myargs['ak_interp'] = 'pressure2sigma'
            ak = pinterp(m2fd, m3f, aksatf, pressure=True)
        else:
            raise TypeError(
                f'Unknown VGTYP={aksatf.VGTYP}; cannot interpolate'
            )
        # Use averaging kernel if available
        pp = (ppm * dp)
        dus = hPa_to_du * pp
        du = (dus * ak).sum(1, keepdims=True)
        myargs['akmethod'] = f'AveragingKernel = {akexpr}'
    else:
        pp = (ppm * dp).sum(1)  # 1e6 hPa
        du = hPa_to_du * pp
        myargs['akmethod'] = 'None'

    # Create a template file from the input
    # for output results
    outf = infd.slice(LAY=0).subsetVariables([key])

    # copy in mean valid cloud fraction as a diagnostic
    outf.copyVariable(m2fd.variables['CFRAC'][:], key='CFRAC')
    if 'PRES' in m3f.variables:
        presf = m3f.subsetVariables(
            ['PRES']
        )
        troppf = presf.mask(
            ~is3valid, dims=('TSTEP', 'LAY', 'ROW', 'COL')
        ).apply(TSTEP='mean', LAY='min')
        stratpmaxf = presf.mask(
            is3valid | istrop, dims=('TSTEP', 'LAY', 'ROW', 'COL')
        ).apply(TSTEP='mean', LAY='max')
        troppf.variables['PRES'][:] += stratpmaxf.variables['PRES'][:]
        troppf.variables['PRES'][:] /= 2
        outf.copyVariable(troppf.variables['PRES'][:], key='TROPOPAUSEPRES')

    # Create a variable to store LAYER count in troposphere
    lays = outf.copyVariable(m2fd.variables['CFRAC'][:], key='LAYERS')
    lays.units = 'count'
    lays.long_name = 'LAYERS'
    lays.var_desc = 'LAYERS'
    lays[:] = 0
    lays[:] = (~ppm.mask).sum(1, keepdims=True)

    # Store column in variable of same name as concentration
    # with dobson units
    ovar = outf.variables[key]
    ovar.units = 'dobson units'
    ovar.var_desc = f'{key} dobson units (molec/cm2 = DU * 2.69e16)'.ljust(80)
    ovar[:, 0] = du
    outf.HISTORY = f'ppm2du(**{myargs})'

    # Save output
    outf.save(outpath, verbose=0)


if __name__ == '__main__':
    args = parser.parse_args()
    ppm2du(
        cpath=args.cpath, m2path=args.m2path, m3path=args.m3path,
        outpath=args.outpath, key=args.key, maxcld=args.maxcld,
        minlsthour=args.minlsthour, maxlsthour=args.maxlsthour,
        satpath=args.satpath, akexpr=args.akexpr, tppexpr=args.tppexpr,
        isvalidexpr=args.isvalidexpr,
    )
