import os
import numpy as np
from scipy.stats import binned_statistic_dd
import PseudoNetCDF as pnc
import gc
import argparse
import sys
import json
from warnings import warn


def getargs(inargs):
    parser = argparse.ArgumentParser()
    parser.add_argument('GRIDDESC', help='Path to GRIDDESC')
    parser.add_argument('GDNAM', help='Grid name in GRIDDESC')
    parser.add_argument('outpath', help='Output path')
    parser.add_argument('optpath', help='Options path for satellite product')
    parser.add_argument('inpaths', nargs='+', help='L2 input files')
    args = parser.parse_args(inargs)
    return args


def openpaths(inpaths, opts):
    """
    Return satellite retrieval file with geolocation and data in one file

    Arguments
    ---------
    inpaths : list
        paths to satellite files

    Returns
    -------
    omf : PseudoNetCDFFile
    """
    if inpaths[0].startswith('https://'):
        omfs = opendappaths(inpaths, opts)
    else:
        omfs = openhe5(inpaths, opts)
    omf = omfs[0].stack(omfs[1:], 'nTimes')
    return omf


def opendapquery(gf, short_name, daterange):
    """
    Returns a list of opendap paths in time/space domain

    Arguments
    ---------
    gf : PseudoNetCDFFile
        must have latitude/longitude variables
    short_name : str
        NASA satellite short name
    daterange : str
        temporal range as defined by CMR
    """
    from urllib.request import urlopen
    from urllib.parse import quote
    glon = gf.variables['longitude']
    glat = gf.variables['latitude']
    daterange = quote(daterange)
    bbox = f'{glon.min()},{glat.min()},{glon.max()},{glat.max()}'
    url = (
        'https://cmr.earthdata.nasa.gov/search/granules.json?' +
        f'short_name={short_name}' +
        f'&temporal[]={daterange}&bounding_box={bbox}&' +
        'page_size=1000&pretty=true'
    )
    r = urlopen(url)
    txt = r.read()
    results = json.loads(txt)
    entries = results['feed']['entry']
    opendappaths = [entry['links'][1]['href'] for entry in entries]
    return opendappaths


def opendappaths(inpaths, opts):
    omfs = []
    for inpath in inpaths:
        tmpf = pnc.pncopen(inpath, format='netcdf')
        omfi = pnc.PseudoNetCDFFile.from_ncvs(
            **{
                varkey: tmpf.variables[varkey]
                for varkey in opts['datakeys'] + opts['geokeys']
            }
        )
        omfs.append(omfi)

    return omfs


def openhe5(inpaths, opts):
    omfs = []
    for inpath in inpaths:
        tmpf = pnc.pncopen(inpath, format='netcdf')
        omfi = pnc.PseudoNetCDFFile.from_ncvs(
            **{
                varkey: tmpf[opts['datagrp']].variables[varkey]
                for varkey in opts['datakeys']
            }
        )
        omgfi = pnc.PseudoNetCDFFile.from_ncvs(
            **{
                varkey: tmpf[opts['geogrp']].variables[varkey]
                for varkey in opts['geokeys']
            }
        )
        datadims = opts.get('datadims', None)
        geodims = opts.get('geodims', None)
        if datadims is None:
            ddims = list(omfi.dimensions)
            datadims = dict(zip(ddims, ['nTimes', 'nXtrack', 'nLevels']))
            print('Dimension mapping heuristically')
            print({dk: len(dv) for dk, dv in omfi.dimensions.items()})
            print('Selected dimension mapping:', datadims)

        if geodims is None:
            gdims = list(omgfi.dimensions)
            geodims = dict(zip(gdims, ['nTimes', 'nXtrack', 'nLevels']))
            print('Dimension mapping heuristically')
            print({dk: len(dv) for dk, dv in omgfi.dimensions.items()})
            print('Selected dimension mapping:', geodims)

        omfi.renameDimensions(**datadims, inplace=True)
        omgfi.renameDimensions(**geodims, inplace=True)

        for geokey in opts['geokeys']:
            omfi.copyVariable(omgfi.variables[geokey], key=geokey)

        omfs.append(omfi)

    return omfs


def getunit(varv):
    attrs = varv.getncatts()
    for ukey in ('Units', 'units'):
        if ukey in attrs:
            return attrs[ukey].strip()
    else:
        return 'unknown'


def process(args):
    gf = pnc.pncopen(args.GRIDDESC, format='griddesc', GDNAM=args.GDNAM)
    outpath = args.outpath
    opts = eval(open(args.optpath, 'r').read())
    datakeys = opts['datakeys']
    if os.path.exists(outpath):
        print('Using cached', outpath, flush=True)
        return

    if args.inpaths[0].startswith('{'):
        dapopts = eval(' '.join(args.inpaths))
        args.inpaths = opendapquery(gf, **dapopts)

    omf = openpaths(args.inpaths, opts)

    for timekey in ['Time', 'time', 'TIME']:
        if timekey in omf.variables:
            tf = omf.subsetVariables([timekey]).renameVariable(timekey, 'time')
            tf.variables['time'].units = (
                "seconds since 1993-01-01 00:00:00+0000"
            )
            break
    else:
        tf = pnc.PseudoNetCDFFile()
        tf.createDimension('time', 1)
        t = tf.createVariable('time', 'd', ('time',))
        t.units = "seconds since 1993-01-01 00:00:00+0000"

    date = tf.getTimes()[0]
    gf.SDATE = int(date.strftime('%Y%j'))
    gf.STIME = 0
    gf.TSTEP = 240000
    LAT = omf.variables['Latitude'][:]
    LON = omf.variables['Longitude'][:]
    i, j = gf.ll2ij(LON, LAT, clean='mask')

    varos = [omf.variables[varkey] for varkey in datakeys]

    gbaddata = eval(opts['grndfilterexpr'], None, omf.variables)
    dbaddata = eval(opts['datafilterexpr'], None, omf.variables)
    baddata = gbaddata | dbaddata

    print('Remove {:.2f}'.format(baddata.mask.mean()))
    # isedge did not change main problem..
    # isedge = np.ones_like(i.mask)
    # isedge[..., 3:-3] = False
    mask2d = baddata.filled(True) | i.mask | j.mask
    if mask2d.all():
        print('No data; skipping', outpath, flush=True)
        return
    else:
        print('Making', outpath, flush=True)

    outf = gf.copy().subsetVariables(['DUMMY'])
    ndim, outshape = sorted([(varo.ndim, varo.shape) for varo in varos])[-1]
    if ndim > 2:
        nk = outshape[-1]
    else:
        nk = 1

    outf.createDimension('LAY', nk)
    twodkeys = []
    renamevars = opts.get('renamevars', {})
    for ki, varkey in enumerate(datakeys):
        outvarkey = renamevars.get(varkey, varkey)
        varv = omf.variables[varkey]
        if 'nTimes' not in varv.dimensions:
            continue
        varo = np.ma.masked_invalid(varv[:])
        print(varkey, outvarkey, flush=True)
        mask = (mask2d.T | varo.mask.T).T
        ol = np.ones(mask.shape)
        myi = np.ma.masked_where(mask, (i.T * ol.T).T).compressed() + 0.5
        myj = np.ma.masked_where(mask, (j.T * ol.T).T).compressed() + 0.5
        if varo.ndim == 2:
            myk = myj * 0 + .5
            twodkeys.append(outvarkey)
        else:
            myk = np.ma.masked_where(
                mask, np.indices(mask.shape)[-1]
            ).compressed() + 0.5

        if varo.ndim <= 3:
            loc = [myk, myj, myi]
            outdims = ('TSTEP', 'LAY', 'ROW', 'COL')
            bins = (
                np.arange(nk + 1), np.arange(gf.NROWS+1), np.arange(gf.NCOLS+1)
            )
        else:
            myk1, myk2 = np.indices(mask.shape)[-2:]
            myk1 = np.ma.masked_where(mask, myk1).compressed() + 0.5
            myk2 = np.ma.masked_where(mask, myk2).compressed() + 0.5
            loc = [myk1, myk2, myj, myi]
            outdims = ('TSTEP', 'LAY', 'LAY', 'ROW', 'COL')
            bins = (
                np.arange(nk + 1), np.arange(nk + 1),
                np.arange(gf.NROWS+1), np.arange(gf.NCOLS+1)
            )

        myvcd = np.ma.masked_where(mask, varo[:]).compressed()
        r = binned_statistic_dd(loc, myvcd, 'mean', bins=bins)
        c = binned_statistic_dd(loc, myvcd, 'count', bins=bins)
        var = outf.createVariable(
            outvarkey, 'f', outdims, missing_value=-9.000E36
        )
        var.var_desc = varkey.ljust(80)
        var.long_name = outvarkey.ljust(16)
        var.units = getunit(varv)
        var[:] = np.ma.masked_invalid(r[0])
        nvar = outf.createVariable(
            'N' + outvarkey, 'f', outdims, missing_value=-9.000E36
        )
        nvar.var_desc = ('Count ' + varkey).ljust(80)
        nvar.long_name = ('N' + outvarkey).ljust(16)
        nvar.units = 'none'
        nvar[:] = c[0]

    delattr(outf, 'VAR-LIST')
    # {dk: slice(None, None, -1) for dk in invertdims}
    pv = omf.variables[opts['pressurekey']]

    pu = getunit(pv).lower()

    if pu in ("hpa", "mb"):
        pfactor = 100.
    elif pu == ("pa", "pascal"):
        pfactor = 1.
    else:
        warn('Unknown unit {}; scale factor = 1'.format(pu))
        pfactor = 1.

    p = pv[:]
    dp = np.diff(p, axis=-1)
    if dp.ndim > 1:
        meanaxis = tuple(range(p.ndim - 1))
        meandp = dp.mean(meanaxis)
    else:
        meandp = dp

    if meandp.mean() > 0:
        p = p[..., ::-1]
        dp = np.diff(p, axis=-1)
        outf = outf.slice(LAY=slice(None, None, -1))
        # 2-D variables have data in layer 0
        # after inverting, it is in layerN
        # it must be inverted again
        for varkey in twodkeys:
            tmpv = outf.variables[varkey]
            tmpv[:] = tmpv[:, ::-1]

    hdp = dp / 2
    pedges = np.ma.concatenate(
        [
            p[..., :-1] - hdp,
            p[..., [-1]] - hdp[..., [-1]],
            p[..., [-1]] + hdp[..., [-1]],
        ],
        axis=-1
    )

    if pedges.ndim > 1:
        meanaxes = tuple(list(range(pedges.ndim-1)))
        pedges1d = pedges.mean(axis=meanaxes)
    else:
        pedges1d = pedges

    pedges1d = np.maximum(0, pedges1d) * pfactor
    ptop = outf.VGTOP = pedges1d[-1]
    psrf = pedges1d[0]
    sigma = (pedges1d[:] - ptop) / (psrf - ptop)

    outf.VGLVLS = sigma[:].astype('f')
    outf.VGTYP = 7
    del outf.variables['DUMMY']
    for k in list(outf.variables):
        klen = len(k)
        if klen > 15:
            print(k, 'too long', len(k))

    if hasattr(outf, 'VAR-LIST'):
        delattr(outf, 'VAR-LIST')

    outf.updatemeta()
    outf.FILEDESC = "cmaqsatproc output"
    outf.HISTORY = sys.argv[0] + ': ' + str(args)
    outf.save(outpath, verbose=1, complevel=1)
    gc.collect()


if __name__ == '__main__':
    args = getargs(None)
    process(args)
