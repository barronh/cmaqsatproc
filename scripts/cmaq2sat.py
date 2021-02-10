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
aa('-v', '--verbose', default=0, action='count', help='Increase verbosity')
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


def findtroposphere(m3f, rawsatf, myargs, verbose=0):
    """
    Identify cells that are within the troposphere

    Arguments
    ---------
    m3f : PseudoNetCDFFile
        Must have either PRES, PV or (TA and ZH)
        * if tppexpr available, then troposphere where greater than PRES
        * elif PV in m3f, then Itahashi 2018 method is used
        * elif TA and ZH in m3f, then using WMO 1957 approach
    rawsatf : PseudoNetCDFFile
        * used to calculate the satellite troposphere pressure using tppexpr
    myargs : dictionary
        must have tppexpr (str or None)
    verbose : int
        count verbosity level

    Returns
    -------
    istrop : array
        boolean is troposphere (True) or not (False)
    """
    tppexpr = myargs['tppexpr']
    if tppexpr is not None and 'PRES' in m3f.variables:
        myargs['istropmethod'] = 'Satellite Tropopause Pressure'
        if verbose > 0:
            print(' * Using:', myargs['istropmethod'])
        tppf = rawsatf.eval(
            f'TropoPausePressure = {tppexpr}'
        ).slice(LAY=0)
        # Satellite TropoPause Pressure surface values only
        stpp = tppf.variables['TropoPausePressure'][:]
        if verbose > 0:
            print(' * TropoPausePressure mean', stpp.mean())

        istrop = m3f.variables['PRES'][:] > stpp
    elif 'PV' in m3f.variables:
        myargs['istropmethod'] = 'Potential Vorticity < 2 (Itahashi 2018)'
        if verbose > 0:
            print(' * Using:', myargs['istropmethod'])
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
        if verbose > 0:
            print(' * Using:', myargs['istropmethod'])
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

    if verbose > 1:
        print(' * IsTrop average layers', istrop.sum(1).mean())

    return istrop


def pinterp(aksatf, PRES=None, PRSFC=None, VGLVLS=None, VGTOP=None, verbose=0):
    """
    Interpolates AveragingKernel from aksatf on satellite layers to pressure
    values per grid cell either defined using PRES or the approximation:

        PRES ~ (PRSFC - VGTOP) * VGLVLS + VGTOP

    Requires either PRES or all of PRSFC, VGLVLS, VGTOP

    Arguments
    ---------
    aksatf : NetCDFFile-like
        has variables dict with AveragingKernel dim(1, NLAYS, NROWS, NCOLS)
    PRES : array or None
        if array, dim(1, NLAYS, NROWS, NCOLS) and units Pascals
    PRSFC : array or None
        if array, dim(1, 1, NROWS, NCOLS) and units Pascals
    VGLVLS: array or None
        if array, shape NLAYS + 1 and units of (P - Ptop) / (Psfc - Ptop)
    VGTOP : array or None
        top of model in Pascals

    Returns
    -------
    ak : array
        Interpolated Averaging Kernel dim(1, NLAYS, NROWS, NCOLS) where
        NLAYS = VGLVLS.size - 1
    """
    approxpres = PRES is None
    if approxpres:
        if verbose > 0:
            print(' * Using: pressure approximated from PRSFC')
        if (VGLVLS is None or VGTOP is None or PRSFC is None):
            errstr = (
                'PRES or PRSFC/VGLVLS/VGTOP are required; got\n'
                + f'PRES={PRES},PRSFC={PRSFC},VGLVLS={VGLVLS},VGTOP={VGTOP}'
            )
            raise ValueError(errstr)

    # The averaging kernel used as y-value (np.interp yp)
    akin = aksatf.variables['AveragingKernel']

    # The output shape will be 3d matching the metcro3d file
    inshape = akin.shape
    nt = inshape[0]
    nj = inshape[2]
    ni = inshape[3]
    outshape = list(inshape)
    if approxpres:
        outshape[1] = VGLVLS.size - 1
        # Out sigma is mid points
        outsigma = (VGLVLS[:-1] + VGLVLS[1:]) / 2
        outvgtop = VGTOP
    else:
        outshape[1] = PRES.shape[1]

    ak = np.ma.masked_all(outshape, dtype='f')

    # Calculate mid point pressures, used as xcoordinate (np.interp xp)
    akpmid = (aksatf.VGLVLS[:-1] + aksatf.VGLVLS[1:]) / 2

    if verbose > 0:
        print(' * Start interpolate AK', flush=True)

    # Iterate over indices (k=0, while t, j, and i change
    for t, j, i in np.ndindex(nt, nj, ni):
        # Hold a vertical column
        akv = akin[t, :, j, i]
        if np.ma.getmaskarray(akv).all():
            # Skip any columns that are 100% masked
            continue

        if approxpres:
            psfc = PRSFC[t, 0, j, i]
            # Calculate pressures at mid-points (np.interp x)
            x = outsigma * (psfc - outvgtop) + outvgtop
        else:
            x = PRES[t, :, j, i]

        if np.ma.getmaskarray(x).all():
            # Skip any columns that are 100% masked
            continue

        ak[t, :, j, i] = np.interp(x, akpmid[::-1], akv[::-1])

    return ak


def ppm2du(
    cpath, m2path, m3path, outpath, key='O3',
    satpath=None, akexpr=None, tppexpr=None, minlsthour=12, maxlsthour=14,
    isvalidexpr=None, maxcld=0.2, verbose=0
):
    """
    Create a satellite-like file from CMAQ with retrieval processing.

    Arguments
    ---------
    cpath : str
        path to IOAPI 3D concentration file containing key
    m2path : str
        path to IOAPI 2D MCIP METCRO2D file containing key
    m3path : str
        path to IOAPI 3D MCIP METCRO3D file containing key
    outpath : str
        path for IOAPI 2D output path
    key : str
        key to use from cpath as a concentration file. Expected to have units
        of either ppm or ppb
    satpath : str or None
        optional, IOAPI 3D file with variables to derive akexpr and/or tppexpr
    akexpr : str or None
        optional, definition of Averaging Kernel which can be a single variable
        or can be dervied from Scattering Weights (w_z) and air mass factor (M)
        e.g., akexpr="ScatWgt / AmfTrop[:, 0]"
    tppexpr : str or None
        optional, definition of Tropopause Pressure in Pascals, which can be a
        single variable name or an expression.
        e.g., tppexpr='TropopausePres[:] * 100'
    minlsthour : int
        lst hour must be >= minlsthour (default: 12)
    maxlsthour : int
        lst hour must be <= maxlsthour (default: 14)
        isoverpass = (lsthour >= minlsthour) & (lsthour <= maxlsthour)
    isvalidexpr : str or None
        optional, expr to use from satpath to select valid cells
    maxcld : float
        Maximum cloud fraction (isclear = CFRAC < maxcld)
    verbose : int
        Count of verbosity level

    Returns
    -------
    outf : NetCDFFile
        File with vertical column density (named `key`), CFRAC,
        TROPOPAUSEPRES, and LAYERS.
        * TFLAG : <YYYYJJJ>, <HHMMSS>
            standard IOAPI variable
        * KEY : dobson units
            vertical column density V_z * w_z / M_z
        * CFRAC : fraction
            fraction of cloudiness in the model
        * TROPOPAUSEPRES : Pascals
            tropopause pressure level (ie. min pressure in troposphere)
        * LAYERS : count
            number of layers used in vertical column
    """
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

    if verbose > 0:
        print(f' * Scaling {cmaqunit} by {unitfactor} to get ppm')

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
    istrop = findtroposphere(m3f, rawsatf, myargs, verbose=verbose)

    # Create a 2D and 3D mask for use later
    is2valid = (isoverpass & iscldfree)

    # Only select data with valid satellite retrieved overpasses
    if isvalidexpr is not None:
        myargs['isvalidmethod'] = f'Paired with satellite {isvalidexpr}'
        if verbose > 0:
            print(' *', myargs['isvalidmethod'])

        isvalidf = rawsatf.eval(
            f'ISVALID = {isvalidexpr}'
        ).slice(LAY=0)
        issatvalid = isvalidf.variables['ISVALID'][:]
        is2valid = is2valid & issatvalid
    else:
        myargs['isvalidmethod'] = 'Model filter only'
        if verbose > 0:
            print(' *', myargs['isvalidmethod'])

    is3valid = is2valid.repeat(inf.NLAYS, 1) & istrop
    is3strat = is2valid.repeat(inf.NLAYS, 1) & (~istrop)

    if verbose > 0:
        print(' * Masking files to remove non-overpass or cloudy pixels')

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

    if verbose > 0:
        print(' * Averaging remaining times.')

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
        if verbose > 0:
            print(f' * Calculating AveragingKernel = {akexpr}')

        # Calculate the averaging Kernel using variables in the raw satellite
        # file
        aksatf = rawsatf.eval(
            f'AveragingKernel = {akexpr}'
        )
        if verbose > 0:
            aktmp = aksatf.variables['AveragingKernel']
            print(' * AveragingKernel shape', aktmp.shape)
            if verbose > 1:
                print(' * AveragingKernel mean profile')
                print(aktmp[:].mean((0, 2, 3)))
            del aktmp
        # Now interpolate the Averaging Kernel to the model levels
        # Currently, this supports
        if (
            aksatf.VGTYP == infd.VGTYP
            and np.allclose(aksatf.VGTYP, infd.VGLVLS)
        ):
            if verbose > 0:
                print(
                    ' * Satellite and CMAQ share vertical coordinate: '
                    + 'no interpolation is necessary.'
                )

            # The vertical coordinate types are the same and the vertical
            # level edges are the same, so no interpolation is necessary.
            akf = aksatf
        elif aksatf.VGTYP == infd.VGTYP and aksatf.VGTYP == 7:
            if verbose > 0:
                print(
                    ' * Both satellite and CMAQ use IOAPI sigma coordinate: '
                    + 'using PseudoNetCDF interpSigma method for interpolation'
                )

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
            if verbose > 0:
                print(
                    ' * Satellite uses pressure coordinate: using pinterp '
                    + 'which assumes CMAQ is on a sigma coordinate (VGTYP=7),'
                    + ' got VGTYP={m3f.VGYP}. This may be okay for the hybrid'
                    + ' coordinate.'
                )

            # This case is designed to support OMNO2d-like pressure coordinates
            # In this case, the VGLVLS will be in decreasing value order in
            # Pascals (e.g., 102500.f, 101500.f, ... , 115.f, 45.f)
            myargs['ak_interp'] = 'pressure2sigma'
            ak = pinterp(
                aksatf, VGLVLS=m3f.VGLVLS, VGTOP=m3f.VGTOP,
                PRSFC=m2fd.variables['PRSFC'], verbose=verbose
            )
            if verbose > 0:
                print(' * AveragingKernel shape', ak.shape)
                if verbose > 1:
                    print(' * AveragingKernel mean profile')
                    print(ak.mean((0, 2, 3)))
        else:
            raise TypeError(
                f'Unknown VGTYP={aksatf.VGTYP}; cannot interpolate'
            )

        if verbose > 0:
            print(
                ' * Calculating vertical column density:\n'
                + f' * {key} = {hPa_to_du} \\sum_z (K_z * ppm_z * dP_z)'
            )

        # Use averaging kernel if available
        pp = (ppm * dp)
        dus = hPa_to_du * pp
        du = (dus * ak).sum(1, keepdims=True)
        myargs['akmethod'] = f'AveragingKernel = {akexpr}'
    else:
        if verbose > 0:
            print(
                ' * Calculating vertical column density:\n'
                + f' * {key} = {hPa_to_du} \\sum_z (ppm_z * dP_z)'
            )

        pp = (ppm * dp).sum(1)  # 1e6 hPa
        du = hPa_to_du * pp
        myargs['akmethod'] = 'None'

    if verbose > 0:
        print(' * Creating output file')

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

    if verbose > 0:
        print(f' * Storing file to disk: {outpath}')

    # Save output
    return outf.save(outpath, verbose=max(0, verbose - 1))


if __name__ == '__main__':
    args = parser.parse_args()
    outf = ppm2du(
        cpath=args.cpath, m2path=args.m2path, m3path=args.m3path,
        outpath=args.outpath, key=args.key, maxcld=args.maxcld,
        minlsthour=args.minlsthour, maxlsthour=args.maxlsthour,
        satpath=args.satpath, akexpr=args.akexpr, tppexpr=args.tppexpr,
        isvalidexpr=args.isvalidexpr, verbose=args.verbose
    )
