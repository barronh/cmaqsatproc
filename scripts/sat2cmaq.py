from glob import glob
import os
from datetime import datetime, timedelta
import numpy as np
from scipy.stats import binned_statistic_dd
import PseudoNetCDF as pnc
import gc
import sys


gdpath = sys.argv[1]
gdnam = sys.argv[2]
outpath = sys.argv[3]
optpath = sys.argv[4]
inpaths = sys.argv[5:]

gf = pnc.pncopen(gdpath, format='griddesc', GDNAM=gdnam)
exec(open(optpath, 'r').read())

def openpaths(inpaths):
    omfs = []
    omgfs = []
    for inpath in inpaths:
        omfi = pnc.pncopen(inpath, format='netcdf')[datagrp]
        omgfi = pnc.pncopen(inpath, format='netcdf')[geogrp]
        omgfi = pnc.PseudoNetCDFFile.from_ncf(omgfi)
        if geovars is not None:
            omgfi = pnc.PseudoNetCDFFile.from_ncvs(
                **omgfi.subsetVariables(geovars).variables
            )
        omfs.append(pnc.PseudoNetCDFFile.from_ncf(omfi))
        omgfs.append(omgfi)
    omf = omfs[0].stack(omfs[1:], datatgtdim)
    omgf = omgfs[0].stack(omgfs[1:], geotgtdim)
    return omf, omgf
    
def process():
    omf, omgf = openpaths(inpaths)
    if os.path.exists(outpath):
        print('Using cached', outpath, flush=True)
        return

    for timekey in ['Time', 'time', 'TIME']:
        if timekey in omgf.variables:
            tf = omgf.subsetVariables([timekey]).renameVariable(timekey, 'time')
            tf.variables['time'].units = "seconds since 1993-01-01 00:00:00+0000"
    else:
        tf = pnc.PseudoNetCDFFile()
        tf.createDimension('time', 1)
        t = tf.createVariable('time', 'd', ('time',))
        t.units = "seconds since 1993-01-01 00:00:00+0000"

    date = tf.getTimes()[0]
    gf.SDATE = int(date.strftime('%Y%j'))
    gf.STIME = int(date.strftime('%H%M%S'))
    gf.TSTEP = 240000
    LAT = omgf.variables['Latitude'][:]
    LON = omgf.variables['Longitude'][:]
    i, j = gf.ll2ij(LON, LAT, clean='mask')

    varos = [omf.variables[varkey] for varkey in varkeys]

    gbaddata = eval(grndfilterexpr, None, omgf.variables)
    dbaddata = eval(datafilterexpr, None, omf.variables)
    baddata = gbaddata | dbaddata

    print('Remove {:.2f}'.format(baddata.mask.mean()))
    # isedge did not change main problem..
    # isedge = np.ones_like(i.mask)
    # isedge[..., 3:-3] = False
    mask2d = baddata.filled(True) | i.mask | j.mask # | isedge
    if mask2d.all():
        print('No data; skipping', outpath, flush=True)
        return
    else:
        print('Making', outpath, flush=True)

    outf = gf.copy().subsetVariables(['DUMMY'])
    ndim, outshape =  sorted([(varo.ndim, varo.shape) for varo in varos])[-1]
    if ndim > 2:
        nk = outshape[-1]
    else:
        nk = 1
    outf.createDimension('LAY', nk)
    for ki, varkey in enumerate(varkeys):
        varv = omf.variables[varkey]
        varo = np.ma.masked_invalid(varv[:])
        print(varkey, flush=True)
        mask = (mask2d.T | varo.mask.T).T
        ol = np.ones(mask.shape)
        myi = np.ma.masked_where(mask, (i.T * ol.T).T).compressed() + 0.5
        myj = np.ma.masked_where(mask, (j.T * ol.T).T).compressed() + 0.5
        if varo.ndim == 2:
            myk = myj * 0 + .5
        else:
            myk = np.ma.masked_where(mask, np.indices(mask.shape)[-1]).compressed() + 0.5

        if varo.ndim <= 3:
            loc = [myk, myj, myi]
            outdims = ('TSTEP', 'LAY', 'ROW', 'COL')
            bins = (np.arange(nk + 1), np.arange(gf.NROWS+1), np.arange(gf.NCOLS+1))
        else:
            myk1, myk2 = np.indices(mask.shape)[-2:]
            myk1 = np.ma.masked_where(mask, myk1).compressed() + 0.5
            myk2 = np.ma.masked_where(mask, myk2).compressed() + 0.5
            loc = [myk1, myk2, myj, myi]
            outdims = ('TSTEP', 'LAY', 'LAY', 'ROW', 'COL')
            bins = (np.arange(nk + 1), np.arange(nk + 1), np.arange(gf.NROWS+1), np.arange(gf.NCOLS+1))
        myvcd = np.ma.masked_where(mask, varo[:]).compressed()
        r = binned_statistic_dd(loc, myvcd, 'mean', bins=bins)
        c = binned_statistic_dd(loc, myvcd, 'count', bins=bins)
        var = outf.createVariable(varkey, 'f', outdims, missing_value=-999)
        var.var_desc = varkey.ljust(80)
        var.long_name = varkey.ljust(16)
        var.units = varv.Units.ljust(16)
        var[:] = np.ma.masked_invalid(r[0])
        nvar = outf.createVariable('N' + varkey, 'f', outdims, missing_value=-999)
        nvar.var_desc = ('N' + varkey).ljust(80)
        nvar.long_name = ('N' + varkey).ljust(16)
        nvar.units = 'none'
        nvar[:] = c[0]
    delattr(outf, 'VAR-LIST')
    p = omf.variables[pressurekey][:]
    pedges = np.append(np.append(p[..., :-1] - np.diff(p)/2, p[..., -1] + np.diff(p)[..., -1]/2), 0)
    ptop = outf.VGTOP = pedges[-1]
    psrf = pedges[0]
    sigma = (pedges[:] - ptop) / (psrf - ptop)
    outf.VGLVLS = sigma[:]
    outf.VGTYP = 7
    outf.subsetVariables(varkeys + ['N' + varkey for varkey in varkeys]).save(outpath, verbose=1)
    gc.collect()

if __name__ == '__main__':
    process()
