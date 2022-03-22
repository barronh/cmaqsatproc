import os
import sys
import time
import numpy as np
import PseudoNetCDF as pnc
import argparse
from warnings import warn

# Constants Copied adopted from The NIST Reference on Constants, Units,
# and Uncertainty. US National Institute of Standards and Technology.
# May 2019. Retrieved 2021-12-08. They are compared to constants from
# CMAQ CONST.EXT
# https://github.com/USEPA/CMAQ/blob/
# c7518fb9449334cb2f8e0226b7ee0f6969059519/CCTM/src/ICL/fixed/const/CONST.EXT

# CMAQ
# mean gravitational acceleration [ m/sec**2 ]
# FSB: Value is mean of polar and equatorial values.
# Source: CRC Handbook (76th Ed) pp. 14-6
# GRAV = 9.80622
#
# Adopted
# https://physics.nist.gov/cgi-bin/cuu/Value?gn Retrieved 2021-12-08
GRAV = 9.80665

# CMAQ
# http://physics.nist.gov/cgi-bin/cuu/Value?na Retrieved 2017-04-21.
# DAVO = 6.02214085774e23
#
# Adopted
# http://physics.nist.gov/cgi-bin/cuu/Value?na Retrieved 2021-12-08
DAVO = 6.02214076e23

# CMAQ
# http://physics.nist.gov/cgi-bin/cuu/Value?r Retrieved 2017-04-21.
# DRGASUNIV = 8.314459848e0
#
# Adopted
# http://physics.nist.gov/cgi-bin/cuu/Value?r Retrieved 2021-12-08.
DRGASUNIV = 8.314462618e0

# standard atmosphere  [ Pa ]
STDATMPA = 101325.0
# Standard Temperature [ K ]
STDTEMP = 273.15


# CMAQ
# mean molecular weight for dry air [ g/mol ]
# FSB: 78.06% N2, 21% O2, and 0.943% A on a mole
# fraction basis ( Source : Hobbs, 1995) pp. 69-70
# MWAIR = 28.9628 / 1e3
#
# Adopted in kg/mole
MWAIR = 28.9628 / 1e3

# https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMO3PR.003/doc/README.OMO3PR.pdf
# http://www.temis.nl/data/conversions.pdf
# assumes mixing ratio in PPM and dp in hPa
hPa_to_du = (
    10 * DRGASUNIV / MWAIR * STDTEMP / GRAV / STDATMPA
)

# Conversion from Dobson units to molecules/m2
DU2MM2 = STDATMPA * DAVO / DRGASUNIV / STDTEMP * 0.01e-3
# Conversion from Dobson units to molecules/cm2
DU2MCM2 = DU2MM2 / 1e4

denseqnstr = f'DU = ppm / 1e6 * DENS / {MWAIR} * DZ * {DAVO} / {DU2MM2:.6e}'
prestaeqnstr = (
    f'DU = ppm / 1e6 * PRES / TA / {DRGASUNIV} * DZ * {DAVO} / {DU2MM2:.6e}'
)
sigmaeqnstr = f'DU = ppm * dp * {hPa_to_du}'


def getparser():
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
        '--todu', choices=['DENS', 'PRESTA', 'SIGMA'], default=None,
        help=(
            'CMAQ units of ppm are converted to dobson using DENS (preferred),'
            + ' PRESTA (pressure and temperature), or sigma. For WRF sigma'
            + ' coordinates, all three should be nearly identical. For WRF'
            + ' hybrid,it is important to use DENS or PRESTA. '
            + f'DENS: {denseqnstr}; PRES: {prestaeqnstr}; SIGMA: {sigmaeqnstr}'
        )
    )
    aa(
        '--minlsthour', default=12, type=float,
        help=(
            'Used with --maxlsthour to restrict times to overpass time. For'
            + ' example, Aura Overpass is 13:40 solar time at the equator and '
            + ' 12:00 solar time at 70 N (minlsthour <= lsthour <= maxlsthour)'
            + '. See Figure 2.4 '
            + 'https://aura.gsfc.nasa.gov/images/aura_validation_v1.0.pdf'
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
        help=(
            'expr for TropoPausePressure (Pa) of tropopause in satellite'
            + ' product'
        )
    )
    aa(
        '--isvalidexpr', default=None, type=str,
        help=(
            'Is valid expr from the satellite product'
            '(e.g., VCD_NO2_Trop.mask == False)'
        )
    )
    aa('--prexpr', default=None, help='Prior expression')
    aa('--swexpr', default=None, help='Scattering weight expression')
    aa('--akexpr', default=None, help='Averaging Kernel expression')
    aa('cpath', help='CMAQ CONC path; must contain key')
    aa('m2path', help='METCRO2D path; contains CFRAC and PRSFC')
    aa('m3path', help='METCRO3D path; contains PV')
    aa('outpath', help='output path')
    aa('key', help='CMAQ key to process')
    return parser


def findtroposphere(
    m3f, rawsatf, myargs, usesat=None, usepv=None, usethermal=None, verbose=0
):
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
        Used to calculate the satellite troposphere pressure using tppexpr.
    usesat : bool or None
        Use the satellite's definition of the tropopause using
        myargs['tppexpr']. If None, usesat is diagnosed by whether tppexpr
        is not None
    usepv : bool or None
        Use the PV approach
        If None, usepv is diagnosed by whether PV is in m3f
    usethermal : bool or None
        Use the thermal-lapse rate approach (aka WMO1957)
        If None, usethermal is diagnosed by whether TA and ZH is in m3f, and is
        only used if usepv is False.
    myargs : dictionary
        must have tppexpr (str or None)
    verbose : int
        count verbosity level

    Returns
    -------
    istrop : array
        boolean is troposphere (True) or not (False)

    Notes
    -----
    Right now, there are four approaches:
     * Find pressures below tropopause pressure from satellite (tppexpr)
     * PV approach based on description in Itahashi (10.5194/acp-20-3373-2020)
       which does not cite, but is likely similar to Holton (10.1029/95RG02097)
     * WMO1957 approach as described in doi:10.1029/2003GL018240
     * Hybrid: either PV or WMO1957

    If usepv and usethermal are set to true, then a Hybrid Approach is used.
    The hybrid approach uses either PV or thermal approach.

    The pros and cons of PV and WMO1957 are discussed in Reichler
    (10.1029/2003GL018240) and Knowland (10.1002/2017GL074532) chooses to use
    the maximum altitude of the two approaches, but used 3 PVU instead of the
    2 PVU typically used and used here.
    """
    tppexpr = myargs['tppexpr']

    if usesat is None:
        usesat = tppexpr is not None and 'PRES' in m3f.variables
    if usepv is None:
        usepv = 'PV' in m3f.variables
    if usethermal is None:
        usethermal = (
            not usepv and 'TA' in m3f.variables and 'ZH' in m3f.variables
        )

    if usesat:
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
    elif usepv and usethermal:
        istroppv = findtroposphere(
            m3f, rawsatf, myargs, usesat=False, usepv=True, usethermal=False,
            verbose=verbose
        )
        istropth = findtroposphere(
            m3f, rawsatf, myargs, usesat=False, usepv=False, usethermal=True,
            verbose=verbose
        )
        myargs['istropmethod'] = 'PV < 2 or Thermal Lapse-rate'
        istrop = istroppv | istropth
    elif usepv:
        myargs['istropmethod'] = 'Potential Vorticity < 2 (Itahashi 2018)'
        if verbose > 0:
            print(' * Using:', myargs['istropmethod'])
        # Any PV < 2 starts the troposphere. Itahashi 2018 submitted
        PV = m3f.variables['PV'][:]
        istrop = np.cumsum(PV[:, ::-1] < 2, axis=1)[:, ::-1] > 0
    elif usethermal:
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


def sinterp(satf, *args, **kwds):
    """
    Thin wrapper to PseudoNetCDF.cmaq.ioapi_base.interpSigma

    Wrapper subsets the file for non-masked values and applies
    the interpolation only over valid data. Because interpSigma
    is computationally expensive, it saves significant amount
    of time.

    Arguments
    ---------
    satf : PseudoNetCDFFile
        expected to have AveragingKernel, ScatteringWeights and/or APriori

    Returns
    -------
    outf : PseudoNetCDFFile
        Same as satf, but on new vertical coordiane system
    """
    # Expected variable
    expected_keys = ['AveragingKernel', 'ScatteringWeights', 'APriori']
    # Find a variable
    for tmpkey in expected_keys:
        if tmpkey in satf.variables:
            break

    # Variable will have dimensions (TSTEP, ..., ROW, COL)
    # and will be masked
    tmpv = satf.variables[tmpkey]

    # Create a nd mask
    mask = tmpv[:].mask
    # Reduce the mask from nd to 2d where if any value is unmasked, the whole
    # row/col combination is processed
    while mask.ndim > 2:
        mask = mask.all(0)

    # Indentify where the mask is not true; these are the indicies to process
    j, i = np.where(~mask)

    # Subset the file for the rows/columns. Name vector ROW
    tmpsatf = satf.slice(ROW=j, COL=i, newdims=('ROW',))

    # Interpolate those values on th sigma
    tmpoutf = tmpsatf.interpSigma(*args, **kwds)

    # Create an output file that is like satf, but with no variables
    outf = satf.subset([])
    # Change the layer dimension to be like the output
    outf.createDimension('LAY', len(tmpoutf.dimensions['LAY']))

    # add the IOAPI vertical coordinate definition
    outf.VGTOP = tmpoutf.VGTOP
    outf.VGLVLS = tmpoutf.VGLVLS

    # Copy the vector values into the ROW/COL locations
    for outkey in expected_keys:
        if outkey in tmpoutf.variables:
            inv = satf.variables[outkey]
            outak = outf.copyVariable(
                inv, key=outkey, withdata=False, fill_value=-999.
            )
            outak[:] = np.ma.masked
            outak[..., j, i] = tmpoutf.variables[outkey][:]

    return outf.mask(values=-999.)


def pinterp(
    aksatf, interpkey, PRES=None, PRSFC=None, VGLVLS=None, VGTOP=None,
    verbose=0
):
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
    interpkey : str
        variable to interpolate

    Returns
    -------
    vals : array
        Interpolated interpkey values dim(1, NLAYS, NROWS, NCOLS) where
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
    akin = aksatf.variables[interpkey]

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


def cmaq2tropvcd(
    cpath, m2path, m3path, outpath, key='O3', todu=None,
    satpath=None, akexpr=None, swexpr=None, prexpr=None, tppexpr=None,
    minlsthour=12, maxlsthour=14, isvalidexpr=None, maxcld=0.2, verbose=0
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
    todu : str
        CMAQ units of ppm are converted to dobson using DENS (preferred),
        PRESTA (pressure and temperature), or sigma. For WRF sigma
        coordinates, all three should be nearly identical. For WRF hybrid,
        it is important to use DENS or PRESTA.
        * DENS requires DENS and ZF variables in m3path
        * PRESTA requires PRES, TA and ZF variables in m3path
        * SIGMA does not require anything in m3path, can have big differences
          from the other methods under specific conditions. Specifically, when
          CMAQ VGTYP is not 7 and predictions cover high altitude (e.g.,
          mountainous) terrains
        * None will select based on variable availability and VGTYP
    satpath : str or None
        optional, IOAPI 3D file with variables to derive akexpr and/or tppexpr
    akexpr : str or None
        optional, definition of Averaging Kernel which can be a single variable
        or can be dervied from Scattering Weights (w_z) and air mass factor (M)
        e.g., akexpr="ScatWgt / AmfTrop[:, 0]"
    swexpr : str or None
        optional, definition of Scattering Weights which can be a single
        variable or can be dervied from the Averaging Kernel (K_z) and air mass
        factor (M) e.g., akexpr="AvgKernel * AmfTrop[:, 0]". When provided,
        cmaq2sat derives the CMAQ-based Amf, which is output as CmaqAmf
    swexpr : str or None
        like swexpr, but for the satellite a priori
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
        (
            isvalidexpr is None and akexpr is None and tppexpr is None
            and swexpr is None and prexpr
        ) and (satpath is not None)
    ):
        warn(
            'Satellite file unused; specify akexpr, swexpr, prexpr, or'
            + ' isvalidexpr'
        )
    elif (
        (
            akexpr is not None or tppexpr is not None or swexpr is not None
            or prexpr is not None
        ) and (satpath is None)
    ):
        raise ValueError(
            'akexpr, swexpr, prexpr, tppexpr or isvalidexpr provided without'
            + ' satellite'
        )

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

    # Define local solar hour using 15 degree approximation
    # of time zone offset as an integer to get lst
    # BHH 2021-12-09: updated to round fractional value rather than truncated
    tzoff = (lon[None, None, :, :] / 15).round(0).astype('i')
    lsth = (utch + tzoff) % 24

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

    dp = -np.diff(pedges, axis=1) / 100
    ppm = infd.variables[key][:] * unitfactor

    if prexpr is not None or swexpr is not None or akexpr is not None:
        aksatf = rawsatf.subset([])
        if prexpr is not None:
            if verbose > 0:
                print(f' * Calculating APriori = {prexpr}')
            prsatf = rawsatf.eval(f'APriori = {prexpr}')
            aksatf.copyVariable(
                prsatf.variables['APriori'], key='APriori'
            )
        if swexpr is not None:
            if verbose > 0:
                print(f' * Calculating ScatteringWeights = {swexpr}')
            swsatf = rawsatf.eval(f'ScatteringWeights = {swexpr}')
            aksatf.copyVariable(
                swsatf.variables['ScatteringWeights'], key='ScatteringWeights'
            )
        if akexpr is not None:
            if verbose > 0:
                print(f' * Calculating AveragingKernel = {akexpr}')
            krsatf = rawsatf.eval(f'AveragingKernel = {akexpr}')
            aksatf.copyVariable(
                krsatf.variables['AveragingKernel'], key='AveragingKernel'
            )
        # Calculate the averaging Kernel using variables in the raw satellite
        # file
        if verbose > 0 and prexpr is not None:
            prtmp = aksatf.variables['APriori']
            print(' * APriori shape', prtmp.shape)
            if verbose > 1:
                print(' * APriori mean profile')
                print(prtmp[:].mean((0, 2, 3)))
            del prtmp
        if verbose > 0 and swexpr is not None:
            swtmp = aksatf.variables['ScatteringWeights']
            print(' * ScatteringWeights shape', swtmp.shape)
            if verbose > 1:
                print(' * ScatteringWeights mean profile')
                print(swtmp[:].mean((0, 2, 3)))
            del swtmp
        if verbose > 0 and akexpr is not None:
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
            # Ensure that missing values are treated as missing
            if prexpr is not None:
                pr = np.ma.masked_values(
                    akf.variables['APriori'][:], -999
                )
            if swexpr is not None:
                sw = np.ma.masked_values(
                    akf.variables['ScatteringWeights'][:], -999
                )
            if akexpr is not None:
                ak = np.ma.masked_values(
                      akf.variables['AveragingKernel'][:], -999
                )
        elif aksatf.VGTYP == infd.VGTYP and aksatf.VGTYP == 7:
            if verbose > 0:
                print(
                    ' * Both satellite and CMAQ use IOAPI sigma coordinate: '
                    + 'using PseudoNetCDF interpSigma method for interpolation'
                )

            # The vertical coordinates are dry sigma pressures, which support
            # the use of the efficient interpSigma method for vertical
            # interpolation
            akf = sinterp(
                aksatf,
                vglvls=infd.VGLVLS, vgtop=infd.VGTOP, interptype='linear',
                extrapolate=False, fill_value='extrapolate', verbose=0
            )
            # Ensure that missing values are treated as missing
            if prexpr is not None:
                myargs['pr_interp'] = 'sigma2sigma'
                pr = akf.variables['APriori'][:]
            if swexpr is not None:
                myargs['sw_interp'] = 'sigma2sigma'
                sw = akf.variables['ScatteringWeights'][:]
            if akexpr is not None:
                myargs['ak_interp'] = 'sigma2sigma'
                ak = akf.variables['AveragingKernel'][:]
        elif aksatf.VGTYP == 4:
            if verbose > 0:
                print(
                    ' * Satellite uses pressure coordinate: using pinterp '
                    + 'which assumes CMAQ is on a sigma coordinate (VGTYP=7),'
                    + ' got VGTYP={m3f.VGYP}. This may be okay for the hybrid'
                    + ' coordinate.'
                )

            pkwds = dict(
                VGLVLS=m3f.VGLVLS, VGTOP=m3f.VGTOP,
                PRSFC=m2fd.variables['PRSFC'], verbose=verbose
            )
            # This case is designed to support OMNO2d-like pressure coordinates
            # In this case, the VGLVLS will be in decreasing value order in
            # Pascals (e.g., 102500.f, 101500.f, ... , 115.f, 45.f)
            if akexpr is not None:
                myargs['ak_interp'] = 'pressure2sigma'
                ak = pinterp(aksatf, 'AveragingKernel', **pkwds)
                if verbose > 0:
                    print(' * AveragingKernel shape', ak.shape)
                if verbose > 1:
                    print(' * AveragingKernel mean profile')
                    print(ak.mean((0, 2, 3)))
            if prexpr is not None:
                myargs['pr_interp'] = 'pressure2sigma'
                pr = pinterp(aksatf, 'APriori', **pkwds)
                if verbose > 0:
                    print(' * APriori shape', pr.shape)
                if verbose > 1:
                    print(' * APriori mean profile')
                    print(pr.mean((0, 2, 3)))
            if swexpr is not None:
                myargs['sw_interp'] = 'pressure2sigma'
                sw = pinterp(aksatf, 'ScatteringWeights', **pkwds)
                if verbose > 0:
                    print(' * ScatteringWeights shape', sw.shape)
                if verbose > 1:
                    print(' * ScatteringWeights mean profile')
                    print(sw.mean((0, 2, 3)))
        else:
            raise TypeError(
                f'Unknown VGTYP={aksatf.VGTYP}; cannot interpolate'
            )

    if myargs['todu'] is None:
        if 'DENS' in m3f.variables and 'ZF' in m3f.variables:
            myargs['todu'] = 'DENS'
        elif (
            inf.VGTYP != 7 and 'PRES' in m3f.variables
            and 'TA' in m3f.variables and 'ZF' in m3f.variables
        ):
            myargs['todu'] = 'PRESTA'
        else:
            myargs['todu'] = 'SIGMA'

    if myargs['todu'] == 'DENS':
        if not ('DENS' in m3f.variables and 'ZF' in m3f.variables):
            raise KeyError(
                'todu option "DENS" requires DENS and ZF variables; use a'
                + ' different option or use a METCRO3D file with DENS and'
                + ' ZF'
            )
        myargs['ppm2du'] = denseqnstr
        densf = m3f.subsetVariables(
            ['DENS', 'ZF']
        )
        densf = densf.mask(
            ~is3valid, dims=('TSTEP', 'LAY', 'ROW', 'COL')
        ).apply(TSTEP='mean')
        DENS = densf.variables['DENS'][:]
        ZF = densf.variables['ZF'][:]
        DZ = ZF.copy()
        DZ[:, 1:] -= ZF[:, :-1]
        cbar = DENS / MWAIR * DZ
        dus = cbar * ppm / 1e6 * DAVO / DU2MM2
    elif myargs['todu'] == 'PRESTA':
        if not (
            'PRES' in m3f.variables and 'TA' in m3f.variables
            and 'ZF' in m3f.variables
        ):
            raise KeyError(
                'todu option "PRESTA" requires PRES, TA and ZF variables; use'
                + ' a different option or use a METCRO3D file with PRES, TA,'
                + ' and ZF'
            )
        myargs['ppm2du'] = prestaeqnstr
        densf = m3f.subsetVariables(
            ['PRES', 'TA', 'ZF']
        )
        densf = densf.mask(
            ~is3valid, dims=('TSTEP', 'LAY', 'ROW', 'COL')
        ).apply(TSTEP='mean')
        ZF = densf.variables['ZF'][:]
        DZ = ZF.copy()
        DZ[:, 1:] -= ZF[:, :-1]
        PRES = densf.variables['PRES'][:]
        TA = densf.variables['TA'][:]
        cbar = PRES[:] / TA[:] / DRGASUNIV * DZ
        dus = cbar * ppm / 1e6 * DAVO / DU2MM2
    elif myargs['todu'] == 'SIGMA':
        myargs['ppm2du'] = sigmaeqnstr
        if verbose > 0:
            print(
                f' * Calculating vertical column density:\n{sigmaeqnstr}'
            )

        if inf.VGTYP != 7:
            print(
                '## Warn: VGTYP is not WRF sigma, but using sigma to'
                + ' approximate delta P'
            )
        # 1e6 hPa
        pp = (ppm * dp)
        dus = hPa_to_du * pp
    else:
        raise KeyError(
            myargs['todu'] + ' todu option unknown, use DENS, PRESTA or SIGMA'
        )

    # Simple VCD no AK
    du = dus.sum(1, keepdims=True)

    if swexpr is not None:
        cmaqAmf = (dus * sw).sum(1, keepdims=True) / du

    if akexpr is not None:
        # Use averaging kernel if available
        # For 2D averaging kernels, the result will be (TSTEP, LAY, ROW, COL)
        #  * typical for NO2, HCHO, SO2 and other optically thin gases
        # For 3D averaging kernels, the prior is required and will be used
        # as a part of the calculation
        #  * typical for O3, and potentially others
        if ak.ndim == 5:
            du = (pr + (ak * (dus - pr)).sum(1)).sum(1, keepdims=True)
            vcd_desc = (
                'Vertical column density using satellite Averaging Kernel'
                + ' (VCD = \\sum_{z2}{\\sum_{z1}{pr_z2 + ak_{z1,z2}'
                + ' * (cmaq_du_{z2} - pr_{z1,z2})}})'
            )
        else:
            du = (dus * ak).sum(1, keepdims=True)
            vcd_desc = (
                'Vertical column density using satellite Averaging Kernel'
                + ' (VCD = \\sum_z{cmaq_du * ak})'
            )
        # This will only be the case for 3D averaging kernels
        myargs['akmethod'] = f'AveragingKernel = {akexpr}'
    else:
        if verbose > 0:
            print(
                ' * Calculating vertical column density:\n'
                + f' * {key} = {hPa_to_du} \\sum_z (ppm_z * dP_z)'
            )
        vcd_desc = (
            'Vertical column density using satellite Averaging Kernel'
            + f' (VCD = \\sum_z({myargs["ppm2du"]}))'
        )
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

    if swexpr is not None:
        # Create a variable to store CmaqAmf
        camfv = outf.copyVariable(
            m2fd.variables['CFRAC'][:], key='CmaqAmfTrop'
        )
        camfv.long_name = 'CmaqAmf'
        camfv.units = 'none'
        camfv.var_desc = (
            'AmfTrop = SlantTropColumnDensity/VerticalTropColumnDensity'
        )
        camfv.description = '\\sum_z(cmaq_du_z * sw_z)'
        camfv[:] = cmaqAmf

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
    ovar.var_desc = (
        f'{key} dobson units (molec/cm2 = DU * {DU2MCM2:.6e})'.ljust(80)
    )
    ovar.description = vcd_desc
    ovar[:, 0] = du
    outf.HISTORY = f'cmaq2tropvcd(**{myargs})'

    if verbose > 0:
        print(f' * Storing file to disk: {outpath}')

    # Save output
    return outf.save(outpath, verbose=max(0, verbose - 1))


if __name__ == '__main__':
    parser = getparser()
    args = parser.parse_args()
    outf = cmaq2tropvcd(**vars(args))
