import os
import numpy as np
from scipy.stats import binned_statistic_dd
import PseudoNetCDF as pnc
import gc
import argparse
import sys
import json
from warnings import warn


def reorderFileDims(infile, dims, inplace=True):
    if inplace:
        out = infile
    else:
        out = infile.copy(variables=False)

    for key, invar in infile.variables.items():
        out.variables[key] = reorderVarDims(invar, dims=dims, key=key)

    return out


def broadcastVar(var, oth):
    othdims = oth.dimensions
    vardimlens = dict(zip(var.dimensions, var.shape))
    othdimlens = dict(zip(oth.dimensions, oth.shape))
    for dk, dl in vardimlens.items():
        if dk not in othdims:
            raise KeyError(f'{dk} not in other variable')
        othdl = othdimlens[dk]
        if dl != othdl:
            raise ValueError(f'{dk}={dl} in var, but {othdl} in oth')

    invar = reorderVarDims(var, othdims)
    outvals = invar[:]
    for di, (dk, dl) in enumerate(zip(othdims, oth.shape)):
        if dk not in invar.dimensions:
            outvals = np.expand_dims(outvals, di).repeat(dl, di)
    return outvals


def reorderVarDims(var, dims, key=None):
    """
    Arguments
    ---------
    var: PseudoNetCDFVariable
    dims : iterable
        iterable of dimension names

    Returns
    -------
    outvar : PseudoNetCDFVariable
        var if dims matches or a new var with matching dims

    Notes
    -----
    """
    vardims = list(var.dimensions)
    newdims = (
        [dk for dk in dims if dk in vardims] +
        [dk for dk in vardims if dk not in dims]
    )
    if newdims == vardims:
        return var

    destax = [di for di, dk in enumerate(newdims)]
    sourceax = [vardims.index(dk) for dk in newdims]
    outdata = np.moveaxis(var[:], sourceax, destax)
    props = var.getncatts()
    if key is None:
        for k in ['Long_name', 'long_name', 'standard_name']:
            if k in props:
                key = getattr(var, k)
                break
        else:
            key = 'unknown'

    outvar = pnc.PseudoNetCDFVariable(
        None, key, outdata.dtype.char, newdims,
        values=outdata
    )
    outvar.setncatts(var.getncatts())
    return outvar


def getargs(inargs):
    parser = argparse.ArgumentParser()
    pa = parser.add_argument
    pa('-v', '--verbose', default=0, action='count')
    pa('-s', '--satpath', default=None, help='Store pregrid output')
    pa('GRIDDESC', help='Path to GRIDDESC')
    pa('GDNAM', help='Grid name in GRIDDESC')
    pa('outpath', help='IOAPI-output path')
    pa('optpath', help='Options path for satellite product')
    pa('inpaths', nargs='+', help='L2 input files')
    args = parser.parse_args(inargs)
    return args


def openpaths(inpaths, opts, verbose=0):
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
        if verbose > 1:
            print('Using opendappaths', flush=True)
        omfs = opendappaths(inpaths, opts, verbose)
    else:
        if verbose > 1:
            print('Using openhe5', flush=True)
        omfs = openhe5(inpaths, opts, verbose)
    omf = omfs[0].stack(omfs[1:], 'nTimes')
    return omf


def opendapquery(gf, verbose, daterange, **qopts):
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
    optstr = '&'.join([
        '{}={}'.format(k, v)
        for k, v in qopts.items()
    ])
    url = (
        'https://cmr.earthdata.nasa.gov/search/granules.json?' +
        optstr +
        f'&temporal[]={daterange}&bounding_box={bbox}&' +
        'page_size=1000&pretty=true'
    )
    if verbose > 0:
        print(url, flush=True)
    r = urlopen(url)
    txt = r.read()
    results = json.loads(txt)
    entries = results['feed']['entry']
    opendappaths = [entry['links'][1]['href'] for entry in entries]
    return opendappaths


def opendappaths(inpaths, opts, verbose):
    omfs = []
    dapdims = opts.get('opendapdims', None)
    for inpath in inpaths:
        if verbose > 1:
            print('Opening', inpath, flush=True)
        tmpf = pnc.pncopen(inpath, format='netcdf')
        omfi = pnc.PseudoNetCDFFile()
        for varkey in opts['datakeys'] + opts['geokeys']:
            if verbose > 2:
                print('Processing', varkey, flush=True)
            tmpv = tmpf.variables[varkey]
            for dim, dimlen in zip(tmpv.dimensions, tmpv.shape):
                if dim not in omfi.dimensions:
                    omfi.createDimension(dim, dimlen)
            dtype = tmpv.dtype
            # Aura OMI data is occasionaly stored as an int16
            # and scaled to a float32
            for propkey in ['scale_factor', 'add_offset']:
                if hasattr(tmpv, propkey):
                    stype = getattr(tmpv, propkey).dtype
                    if (
                        dtype.char in ('i', 'h') and
                        stype.char not in ('i', 'h')
                    ):
                        dtype = stype

            omfi.copyVariable(tmpv, key=varkey, dtype=dtype)

        if dapdims is not None:
            omfi.renameDimensions(**dapdims, inplace=True)

        omfs.append(omfi)

    return omfs


def openhe5(inpaths, opts, verbose):
    omfs = []
    for inpath in inpaths:
        if verbose > 1:
            print('Opening', inpath, flush=True)
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
            print('Dimension mapping heuristically', flush=True)
            print({dk: len(dv) for dk, dv in omfi.dimensions.items()})
            print('Selected dimension mapping:', datadims, flush=True)

        if geodims is None:
            gdims = list(omgfi.dimensions)
            geodims = dict(zip(gdims, ['nTimes', 'nXtrack', 'nLevels']))
            print('Dimension mapping heuristically', flush=True)
            print({dk: len(dv) for dk, dv in omgfi.dimensions.items()})
            print('Selected dimension mapping:', geodims, flush=True)

        for key in datadims:
            if key not in omfi.dimensions:
                print(f'Key {key} not found; in {omgi.dimensions}')

        omfi.renameDimensions(**datadims, inplace=True)

        for key in geodims:
            if key not in omgfi.dimensions:
                print(f'Key {key} not found; in {omgfi.dimensions}')

        omgfi.renameDimensions(**geodims, inplace=True)

        for geokey in opts['geokeys']:
            omfi.copyVariable(omgfi.variables[geokey], key=geokey)

        flipdimkeys = opts.get('flipdims', [])
        if len(flipdimkeys) > 0:
            flipslices = {
                k: slice(None, None, -1)
                for k in flipdimkeys if k in omfi.dimensions
            }
            omfi = omfi.sliceDimensions(**flipslices)
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
    if args.verbose > 1:
        print(f'Opening {args.GRIDDESC} GDNAM={args.GDNAM}', flush=True)
    gf = pnc.pncopen(args.GRIDDESC, format='griddesc', GDNAM=args.GDNAM)
    outpath = args.outpath
    if os.path.exists(outpath):
        print('Using cached', outpath, flush=True)
        return

    opts = eval(open(args.optpath, 'r').read())
    omf = subset(args, gf, opts)
    outf = grid(args, gf, opts, omf)
    outf.save(outpath, verbose=1, complevel=1)


def subset(args, gf, opts):
    """
    Arguments
    ---------
    args: namespace 
        must have inpaths, verbose, grndfilterexpr, datafilterexpr
        satpath and any requirements of openpaths
    gf : pnc.PseudoNetCDFFile
        griddesc file that implements IOAPI
    opts : mappable
        Product specific options

    Returns
    -------
    outf : PseudoNetCDFFile
        dimensions nTimes, nXtrack, and nLevels with datakeys dn geokeys
    """
    latkey = opts.get('Latitude', 'Latitude')
    lonkey = opts.get('Longitude', 'Longitude')
    if args.inpaths[0].startswith('{'):
        dapopts = eval(' '.join(args.inpaths))
        if args.verbose > 1:
            print('Dap options {}'.format(dapopts), flush=True)
        args.inpaths = opendapquery(gf, verbose=args.verbose, **dapopts)

    if args.verbose > 1:
        print(f'Opening inpaths {args.inpaths}', flush=True)
    omf = openpaths(args.inpaths, opts, verbose=args.verbose)

    if args.verbose > 1:
        print(f'Calculating i/j from {lonkey}/{latkey}')
    LAT = omf.variables[latkey][:]
    LON = omf.variables[lonkey][:]
    i, j = gf.ll2ij(LON, LAT, clean='mask')

    gbaddata = eval(opts['grndfilterexpr'], None, omf.variables)
    dbaddata = eval(opts['datafilterexpr'], None, omf.variables)
    baddata = gbaddata | dbaddata

    if args.verbose > 0:
        print('Bad data mask {:.2f}'.format(baddata.mask.mean()), flush=True)
        print('Bad data pct  {:.2f}'.format(baddata.mean()), flush=True)
        print(
            'Bad of valid  {:.2f}'.format(baddata.compressed().mean()),
            flush=True
        )

    omdims = ('nTimes', 'nXtrack')[:i.ndim]
    # isedge did not change main problem..
    # isedge = np.ones_like(i.mask)
    # isedge[..., 3:-3] = False
    mask2d = baddata.filled(True) | i.mask | j.mask
    omf.createVariable('BADDATA', 'i', omdims, values=mask2d.astype('i'))
    if mask2d.ndim == 1:
        badrow = mask2d
    elif mask2d.ndim == 2:
        badrow = mask2d.all(1)

    omf = omf.slice(nTimes=np.where(~badrow))
    mask2d = omf.variables['BADDATA'] == 1
    omf.HISTORY = sys.argv[0] + ': ' + str(args)
    omfm = omf.mask(where=mask2d, dims=('nTimes', 'nXtrack'))
    if 'nLevels' in omfm.dimensions:
        nz = len(omfm.dimensions['nLevels'])
        mask3d = mask2d[..., None].repeat(nz, 2)
        omfm = omfm.mask(where=mask3d, dims=('nTimes', 'nXtrack', 'nLevels'))
        
    if args.satpath is not None:
        omfm.save(args.satpath, complevel=1, verbose=1, format='NETCDF4')

    return omfm


def grid(args, gf, opts, omf):
    """
    Arguments
    ---------
    args: namespace 
        must have inpaths, verbose, grndfilterexpr, datafilterexpr
        satpath and any requirements of openpaths
    gf : pnc.PseudoNetCDFFile
        griddesc file that implements IOAPI
    opts : mappable
        Product specific options
    omf : pnc.PseudoNetCDFFile
        subset of data with masks applied

    Returns
    -------
    outf : PseudoNetCDFFile
        dimensions nTimes, nXtrack, and nLevels with datakeys dn geokeys
    """
    outpath = args.outpath
    datakeys = opts['datakeys']
    outkeys = opts.get('outkeys', datakeys)
    latkey = opts.get('Latitude', 'Latitude')
    lonkey = opts.get('Longitude', 'Longitude')
    timekey = opts.get('Time', 'Time')
    if args.verbose > 1:
        print(f'Calculating time', flush=True)
    for tkey in [timekey, 'Time', 'time', 'TIME']:
        if tkey in omf.variables:
            tf = omf.subsetVariables([tkey]).renameVariable(tkey, 'time')
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
    LAT = omf.variables[latkey][:]
    LON = omf.variables[lonkey][:]
    i, j = gf.ll2ij(LON, LAT, clean='mask')

    mask2d = omf.variables['BADDATA'][:] == 1
    if mask2d.all():
        print('No data; skipping', outpath, flush=True)
        return
    else:
        print('Making', outpath, flush=True)

    if args.verbose > 0:
        utchour = np.array([t.hour for t in tf.getTimes()])
        localhour = np.ma.masked_where(
            mask2d,
            utchour[:, None] + omf.variables[lonkey][:] / 15
        )
        ptiles = [0, 10, 25, 75, 90, 100]
        localhourpct = np.percentile(localhour.compressed(), ptiles)
        localhourpctstr = ' '.join([
            '{:5.2f}'.format(h) for h in localhourpct
        ])
        ptilestr = ' '.join([
            '{:5d}'.format(p) for p in ptiles
        ])
        print('Percentiles:', ptilestr)
        print('Local Time :', localhourpctstr)

    outf = gf.copy().subsetVariables(['DUMMY'])

    if 'nLevels' in omf.dimensions:
        nk = len(omf.dimensions['nLevels'])
    else:
        nk = 1

    outf.createDimension('LAY', nk)
    twodkeys = []
    renamevars = opts.get('renamevars', {})
    for ki, varkey in enumerate(outkeys):
        outvarkey = renamevars.get(varkey, varkey)
        if args.verbose > 1:
            print(f'Masking and gridding {varkey} as {outvarkey}', flush=True)
        varv = omf.variables[varkey]
        if 'nTimes' not in varv.dimensions:
            continue

        varo = np.ma.masked_invalid(
            reorderVarDims(varv, ('nTimes', 'nXtrack'), key=varkey)[:]
        )

        varmask = varo.mask
        mask = broadcastVar(mask2d, varo)
        if mask2d.shape == varmask.shape[:mask2d.ndim]:
            mask = (mask2d.T | varmask.T).T
        elif mask2d.shape == varmask.shape[-mask2d.ndim:]:
            mask = (mask2d | varmask)
        else:
            raise ValueError(
                f'Masks not aligned {mask2d.shape} and {varmask.shape}'
            )
        ol = np.ones(mask.shape)
        myi = np.ma.masked_where(mask, (i.T * ol.T).T).compressed() + 0.5
        myj = np.ma.masked_where(mask, (j.T * ol.T).T).compressed() + 0.5
        if varo.ndim <= 2:
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
    if args.verbose > 1:
        print('Calculating pressure for sigma approximation', flush=True)
    if opts['pressurekey'] is None:
        p = np.array([50000], dtype='f')
        pedges = np.array([101325, 0.], dtype='f')
        pedges1d = np.array([101325, 0], dtype='f')
    else:
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

        if 'nLevelEdges' in pv.dimensions:
            pedges = pv[:]
        else:
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

        # Ensure pedges is never negative
        # heuristic top identification could cause that problem.
        pedges1d = np.maximum(0, pedges1d) * pfactor

    ptop = outf.VGTOP = pedges1d[-1]
    psrf = pedges1d[0]

    if pedges.ndim == 1:
        # OMI ScatteringWtPressure is on a pressure
        # grid.
        outf.VGTYP = 4
        outf.VGLVLS = pedges1d.astype('f')
    else:
        # OMPROFOZ ProfileLevelPressure
        # is on a hybrid sigma/eta coordinate
        # sigma approximation is being used.
        sigma = (pedges1d[:] - ptop) / (psrf - ptop)
        outf.VGTYP = 7
        outf.VGLVLS = sigma[:].astype('f')

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

    gc.collect()
    return outf


if __name__ == '__main__':
    args = getargs(None)
    process(args)
