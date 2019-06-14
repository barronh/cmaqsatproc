import os
import sys
import numpy as np
import PseudoNetCDF as pnc


if len(sys.argv) < 6:
    print("""Usage: {} CONCPATH METCRO2D METCRO3D OUTPATH KEY [AverKern]

Description:
    Convert 3-D IOAPI mixing ratios (ppm) to dobson units within the
    tropopsphere as diagnosed by potential vorticity < 2 (see Itahashi 2018
    submitted).
    
    * Cloud Fraction (METCRO2D CFRAC) is used to screen out cloudy pixels (CFRAC > 0.2).
    * CMAQ times between 12 and 14 LST are eligible for comparison the satellite.

Arguments:
    CONCPATH : IOAPI CONC file must contain KEY
    METCRO2D : IOAPI 2 dimensional met file must contain CFRAC and PRSFC
    METCRO3D : IOAPI 3 dimensional met file must contain PV
    OUTPATH  : path for output (e.g., test.nc)
    KEY      : variable from CONCPATH for translation
    AKPATH   : IOAPI Averaging Kernel file (optional, must have AveragingKernel)
    
Example No Averaging Kernel:
    python cmaq2sat.py CMAQ.ACONC_20160101.nc MECRO2D_160101 METCRO3D_160101 CMAQ_L3_20160101.nc O3

Example No Averaging Kernel:
    python cmaq2sat.py CMAQ.ACONC_20160101.nc MECRO2D_160101 METCRO3D_160101 CMAQ_L3_20160101.nc O3 AVGK.nc

Where AVGK.nc was made by gridding the AveragingKernel from an L2 file and averaging to 1-day
""")
    exit()

cpath = sys.argv[1]
m2path = sys.argv[2]
m3path = sys.argv[3]
outpath = sys.argv[4]
key = sys.argv[5]
if len(sys.argv) > 6:
    akpath = sys.argv[6]
else:
    akpath = None

if len(sys.argv) > 7:
    akkey = sys.argv[7]
else:
    akkey = None


def ppm2du(cpath, m2path, m3path, outpath, key='O3', akpath=None, akkey=None):
    # outpath = os.path.join('du', os.path.basename(cpath))
    if os.path.exists(outpath):
        print('Using cached:', outpath)
        return

    inf = pnc.pncopen(cpath, format='ioapi')
    m2f = pnc.pncopen(m2path, format='ioapi')
    m3f = pnc.pncopen(m3path, format='ioapi')
    if akpath is not None:
        if akkey is None:
            akkey = 'AveragingKernel'
        akf = pnc.pncopen(akpath, format='ioapi').eval(
            'AveragingKernel = {}'.format(akkey)
        ).interpSigma(
            inf.VGLVLS, vgtop=inf.VGTOP, interptype='linear',
            extrapolate=False, fill_value='extrapolate', verbose=0
        )
        ak = np.ma.masked_values(akf.variables['AveragingKernel'][:], -999)
    cf = m2f.variables['CFRAC'][:]
    PV = m3f.variables['PV'][:]
    j = np.arange(inf.NROWS)
    i = np.arange(inf.NCOLS)
    I, J = np.meshgrid(i, j)
    lon, lat = inf.ij2ll(I, J)
    tzoff = (lon[None, None, :, :] // 15).round(0).astype('i')
    t = np.arange(0, 25)[:, None, None, None] + tzoff

    isoverpass = (t >= 12) & (t <= 14)
    iscldfree = cf < 0.2

    # Any PV < 2 starts the troposphere. Itahashi 2018 submitted
    istrop = np.cumsum(PV[:, ::-1] < 2, axis=1)[:, ::-1] > 0
    
    # Create a 2D and 3D mask for use later
    is3valid = (isoverpass & iscldfree).repeat(inf.NLAYS, 1) & istrop
    is2valid = (isoverpass & iscldfree)

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

    # Average only valid pixels
    inf = infm.apply(TSTEP='mean')
    m2f = m2fm.apply(TSTEP='mean')

    # calculate edge pressures from surface, sigma, and top pressure
    pedges =  (
        (m2f.variables['PRSFC'][:] - inf.VGTOP) *
        inf.VGLVLS[None, :, None, None] + inf.VGTOP
    )
    # https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMO3PR.003/doc/README.OMO3PR.pdf
    # http://www.temis.nl/data/conversions.pdf
    # assumes mixing ratio in PPM and dp in hPa
    hPa_to_du = (
        10 * 1.3807e-23 * 6.022e23 / 0.02894 * 273.15 / 9.80665 / 101325.
    )

    dp = -np.diff(pedges, axis=1) / 100
    ppm = inf.variables[key][:].copy()

    # Using O3 gradient as tropopause
    # commented out to use PV instead
    # dppm = np.concatenate([ppm[:,[0]] * 0, np.diff(ppm, axis=1)], axis=1)
    # mask = (dppm > 0.050).astype('i')
    # cmask = np.cumsum(mask, axis=1) > 0
    # ppm = np.ma.masked_where(cmask, ppm)
    if akpath is not None:
        # Use averaging kernel if available
        pp = (ppm * dp)
        dus = hPa_to_du * pp
        du = (dus * ak).sum(1, keepdims=True)
    else:
        pp = (ppm * dp).sum(1)  # 1e6 hPa
        du = hPa_to_du * pp

    # Create a template file from the input
    # for output results
    outf = inf.slice(LAY=0).subsetVariables([key])
    # copy in mean valid cloud fraction as a diagnostic
    outf.copyVariable(m2f.variables['CFRAC'][:], key='CFRAC')
    # Create a variable to store LAYER count in troposphere
    lays = outf.copyVariable(m2f.variables['CFRAC'][:], key='LAYERS')
    lays.units = 'count'
    lays.long_name = 'LAYERS'
    lays.var_desc = 'LAYERS'
    lays[:] = 0
    lays[:] = (~ppm.mask).sum(1, keepdims=True)

    # Store column in variable of same name as concentration
    # with dobson units
    ovar = outf.variables[key]
    ovar.units = 'dobson units'
    ovar[:, 0] = du

    # Save output
    outf.save(outpath, verbose=0)


ppm2du(cpath, m2path, m3path, outpath, key=key, akpath=akpath, akkey=akkey)
