__all__ = ['open_griddesc', 'open_ioapi']

import xarray as xr

known_overpasstimes = dict(
    aura=13.75,
    terra=10.5,
    aqua=13.5,
    metop_am=9.5,
    metop_pm=21.5,
)

default_griddesc_txt = b"""' '
'LATLON'
  1  0.0 0.0 0.0 0.0 0.0
'POLSTE_HEMI'
  6         1.000        45.000       -98.000       -98.000        90.000
'LamCon_40N_97W'
  2        33.000        45.000       -97.000       -97.000        40.000
'LamCon_25N_95W'
  2        25.000        25.000       -95.000       -95.000        25.000
' '
'US_1deg'
'LATLON'              -140.00        20.0      1.0      1.0   90   40 1
'US_0pt1deg'
'LATLON'              -140.00        20.0      0.1      0.1  900  400 1
'global_1deg'
'LATLON'              -180.00       -90.0      1.0      1.0  360  180 1
'global_0pt1deg'
'LATLON'              -180.00       -90.0      0.1      0.1 3600 1800 1
'global_2x2.5'
'LATLON'              -181.25       -89.0      2.5      2.0  144   89 1
'global_4x5'
'LATLON'              -182.50       -88.0      5.0      4.0   72   44 1
'108NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 108000.0 108000.0  187  187 1
'324NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 324000.0 324000.0   63   63 1
'1188NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 1188000. 1188000.   17   17 1
'972US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0 972000.0 972000.0   6    4 1
'324US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0 324000.0 324000.0   17   12 1
'108US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0 108000.0 108000.0   51   34 1
'36US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0  36000.0  36000.0  153  100 1
'12US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0  12000.0  12000.0  459  299 1
'4US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0   4000.0   4000.0 1377  897 1
'1US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0   1000.0   1000.0 5508 3588 1
'12US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0  12000.0  12000.0  396  246 1
'4US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0   4000.0   4000.0 1188  738 1
'1US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0   1000.0   1000.0 4752 2952 1
'36US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0  36000.0  36000.0  172  148 1
'12US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0  12000.0  12000.0  516  444 1
'4US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0   4000.0   4000.0 1548 1332 1
'1US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0   1000.0   1000.0 4644 3996 1
'108US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0 108000.0 108000.0   60   50 1
'NAQFC_CONUS'
'LamCon_25N_95W' -4226153.11044303 -834746.472325356 5079.0 5079.0  1473 1025 1
' '"""


def griddesc(griddesc_txt):
    from collections import OrderedDict
    gddefns = OrderedDict()
    prjdefns = OrderedDict()
    # lines with ' ' are separators
    reallines = [
        l for l in griddesc_txt.split('\n') if l.strip() not in ("''", "' '")
    ]
    # All definitions come in pairs of lines. The first is the name and the
    # second is either all numeric for projections or, for grids, a string
    #  projection name followed by numeric grid definitions
    prjattrkeys = ['GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT']
    gdattrkeys = [
        'PRJNAME', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS',
        'NROWS', 'NTHIK'
    ]
    for name, defn in zip(reallines[0:-1:2], reallines[1::2]):
        name = eval(name)
        values = defn.split()
        if "'" in defn:
            gddefn = dict(zip(gdattrkeys, values))
            gddefn['PRJNAME'] = eval(gddefn['PRJNAME'])
            gddefn['GDNAM'] = name
            gddefn['XORIG'] = float(gddefn['XORIG'])
            gddefn['YORIG'] = float(gddefn['YORIG'])
            gddefn['XCELL'] = float(gddefn['XCELL'])
            gddefn['YCELL'] = float(gddefn['YCELL'])
            gddefn['NCOLS'] = int(gddefn['NCOLS'])
            gddefn['NROWS'] = int(gddefn['NROWS'])
            gddefn['NTHIK'] = int(gddefn['NTHIK'])
            gddefns[name] = gddefn
        else:
            prjdefn = dict(zip(prjattrkeys, values))
            prjdefn['PRJNAME'] = name
            prjdefn['GDTYP'] = int(prjdefn['GDTYP'])
            prjdefn['P_ALP'] = float(prjdefn['P_ALP'])
            prjdefn['P_BET'] = float(prjdefn['P_BET'])
            prjdefn['P_GAM'] = float(prjdefn['P_GAM'])
            prjdefn['XCENT'] = float(prjdefn['XCENT'])
            prjdefn['YCENT'] = float(prjdefn['YCENT'])
            prjdefns[name] = prjdefn
    for name, defn in gddefns.items():
        defn.update(prjdefns[defn['PRJNAME']])
    return gddefns


def griddesc_from_attrs(attrs):
    """
    Create an xarray Dataset that defines a CMAQ grid using attributes
    supplied by the user. A minimum number of attributes is required:
    GDTYP, P_ALP, P_BET, P_GAM, XCENT, YCENT
    GDNAM, NCOLS, NROWS, XCELL, YCELL, XORIG, YORIG, NTHIK

    For attribute definitions, see https://www.cmascenter.org/ioapi/
    documentation/all_versions/html/GRIDS.html.

    Arguments
    ---------
    attrs : dict
        IOAPI attritbutes (e.g, attrs=dict(NCOLS=3, NROWS=2, ...))
    GDNAM: str
        Name of grid as defined in GRIDDESC file
    gdpath : str
        Path to GRIDDESC file. If None (default), use contents of
        default_griddesc_txt instead. default_griddesc_txt can be modified
        directly, although that is discouraged.

    Returns
    ---------
    gf : xarray.Dataset
        File with coordinates based on GDNAM in gdpath
    """
    import numpy as np
    attrs['crs'] = get_proj4string(attrs)
    outf = xr.Dataset(
        data_vars=dict(),
        coords=dict(
            ROW=np.arange(attrs['NROWS']) + 0.5,
            COL=np.arange(attrs['NCOLS']) + 0.5,
        ), attrs=attrs
    )
    return outf


def open_griddesc(GDNAM, gdpath=None):
    """
    Create an xarray Dataset that defines a CMAQ grid using GRIDDESC

    Arguments
    ---------
    GDNAM: str
        Name of grid as defined in GRIDDESC file
    gdpath : str
        Path to GRIDDESC file. If None (default), use contents of
        default_griddesc_txt instead. default_griddesc_txt can be modified
        directly, although that is discouraged.

    Returns
    ---------
    gf : xarray.Dataset
        File with coordinates based on GDNAM in gdpath
    """
    import os
    if gdpath is None:
        griddesc_txt = default_griddesc_txt.decode()
    else:
        if os.path.exists(gdpath):
            griddesc_txt = open(gdpath).read()
        else:
            griddesc_txt = gdpath

    gdattrs = griddesc(griddesc_txt)
    attrs = gdattrs[GDNAM]
    outf = griddesc_from_attrs(attrs)
    return outf


def get_proj4string(attrs):
    """
    Return a proj4 string from IOAPI attributes supplied by the user.

    Typical attributes required include:
    GDTYP, P_ALP, P_BET, P_GAM, XCENT, YCENT
    GDNAM, NCOLS, NROWS, XCELL, YCELL, XORIG, YORIG, NTHIK

    The file can have earth_radius to explicitly set the radius of the
    earth.

    For attribute definitions, see https://www.cmascenter.org/ioapi/
    documentation/all_versions/html/GRIDS.html.
    """
    import copy
    import os

    popts = copy.copy(attrs)
    ENV_IOAPI_ISPH = os.environ.get('IOAPI_ISPH', '6370000.')
    popts['R'] = popts.get('earth_radius', float(ENV_IOAPI_ISPH))
    popts['x_0'] = -attrs['XORIG']
    popts['y_0'] = -attrs['YORIG']

    if popts['GDTYP'] == 1:
        projtmpl = '+proj=lonlat +lat_0={YCENT} +lon_0={XCENT} +R={R} +no_defs'
    elif popts['GDTYP'] == 2:
        projtmpl = (
            '+proj=lcc +lat_0={YCENT} +lon_0={P_GAM} +lat_1={P_ALP}'
            + ' +lat_2={P_BET} +x_0={x_0} +y_0={y_0} +R={R} +to_meter={XCELL}'
            + ' +no_defs'
        )
    elif popts['GDTYP'] == 6:
        popts['lat_0'] = popts['P_ALP'] * 90
        projtmpl = (
            '+proj=stere +lat_0={lat_0} +lat_ts={P_BET} +lon_0={P_GAM}'
            + ' +x_0={x_0} +y_0={y_0} +R={R} +to_meter={XCELL} +no_defs'
        )
    elif popts['GDTYP'] == 7:
        projtmpl = (
            '+proj=merc +lat_ts=0 +lon_0={XCENT}'
            + ' +x_0={x_0} +y_0={y_0} +R={R} +to_meter={XCELL} +no_defs'
        )
    else:
        return None
    return projtmpl.format(**popts)


def open_ioapi(path, **kwargs):
    """
    Open an IOAPI file in NetCDF format using xarray and construct coordinate
    variables. (time based on TFLAG or properties, ROW/COL in projected space,
    and LAY based on VGTYP, VGLVLS, and VGTOP)

    Arguments
    ---------
    path : str
        Path to the IOAPI file in NetCDF format
    kwargs : mappable
        Passed to xr.open_dataset(path, **kwargs)

    Returns
    ---------
    qf : xarray.Dataset
        File with data and coordinates based on path
    """
    import xarray as xr

    qf = xr.open_dataset(path, **kwargs)
    qf.csp.add_coords(inplace=True)

    return qf


@xr.register_dataset_accessor("csp")
class CmaqSatProcAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def add_coords(self, inplace=True):
        """
        Add TSTEP, LAY, ROW, and COL coordinates based file.
        """
        import pandas as pd
        import warnings
        import numpy as np
        from datetime import datetime

        if inplace:
            qf = self._obj
        else:
            qf = self._obj.copy()

        if 'TFLAG' in qf.data_vars:
            tflag = qf.data_vars['TFLAG'][:, 0].values.copy()
            if (tflag[:, 0] < 1).all():
                tflag[:, 0] = 1970001
            times = pd.to_datetime([
                datetime.strptime(f'{JDATE}T{TIME:06d}', '%Y%jT%H%M%S')
                for JDATE, TIME in tflag
            ])
        elif 'TSTEP' in qf.coords:
            tstep = qf.coords['TSTEP']
            if tstep.size == 1:
                tstep = [tstep]
            times = [
                datetime(
                    t.dt.year, t.dt.month, t.dt.day, t.dt.hour, t.dt.minute,
                    t.dt.second
                ) for t in tstep
            ]
        else:
            SDATE = qf.SDATE
            if SDATE < 1:
                SDATE = 1970001
            date = datetime.strptime(f'{SDATE}T{qf.STIME:06d}', '%Y%jT%H%M%S')
            nt = qf.dims['TSTEP']
            dm = (qf.attrs['TSTEP'] % 10000) // 100
            ds = (qf.attrs['TSTEP'] % 100)
            dh = qf.attrs['TSTEP'] // 10000 + dm / 60 + ds / 3600.
            if nt == 1 and dh == 0:
                dh = 1
            times = pd.date_range(date, periods=nt, freq=f'{dh}H')

        qf.coords['TSTEP'] = times
        if 'VGLVLS' in qf.attrs:
            qf.coords['LAY'] = (qf.VGLVLS[1:] + qf.VGLVLS[:-1]) / 2

        crs = get_proj4string(qf.attrs)
        if crs is None:
            warnings.warn((
                'Unknown project ({GDTYP}); currently support lonlat (1),'
                + ' lcc (2), polar stereograpic (6), equatorial mercator (7)'
            ).format(GDTYP=qf.attrs['GDTYP']))
        else:
            qf.attrs['crs'] = crs
            row = np.arange(qf.dims['ROW']) + 0.5
            col = np.arange(qf.dims['COL']) + 0.5
            if qf.GDTYP == 1:
                row = row * qf.attrs['YCELL'] + qf.attrs['YORIG']
                col = col * qf.attrs['XCELL'] + qf.attrs['XORIG']
            qf.coords['ROW'] = row
            qf.coords['COL'] = col

        return qf

    def to_ioapi(self, reset_index=True, drop=True):
        """
        Infer standard IOAPI properties (including TFLAG) as possible.
        """
        import numpy as np
        import pandas as pd
        import xarray as xr

        now = pd.to_datetime('now')
        jnow = np.int32(now.strftime('%Y%j'))
        tnow = np.int32(now.strftime('%H%M%S'))
        outf = self._obj.copy()

        for k in outf.data_vars:
            if k == 'TFLAG':
                continue
            outvar = self._obj[k]
            dims = outvar.dims
            if 'LAY' not in dims:
                outvar = outvar.expand_dims(LAY=1)
            if 'TSTEP' not in dims:
                outvar = outvar.expand_dims(TSTEP=1)
            outf[k] = outvar.transpose('TSTEP', 'LAY', 'ROW', 'COL')
            defattrs = dict(
                units='unknown', long_name=k.ljust(16), var_desc=k.ljust(80)
            )
            defattrs.update(outf[k].attrs)
            outf[k].attrs.update(defattrs)

        outkeys = [k for k in outf.data_vars if k != 'TFLAG']
        nv = np.int32(len(outkeys))

        nl = np.int32(outf.dims.get('LAY', 1))
        vglvls = np.linspace(0, 1, nl + 1, dtype='f')[::-1]
        vgtyp = np.int32(-9999)
        vgtop = np.float32(5000)
        sdate = np.int32(-635)

        defattrs = {
            'EXEC_ID': 'NA'.ljust(80), 'IOAPI_VERSION': 'NA'.ljust(80),
            'UPNAM': 'cmaqsatproc'.ljust(16), 'FTYPE': np.int32(1), 'NLAYS': nl,
            'VGTOP': vgtop, 'VGLVLS': vglvls, 'VGTYP': vgtyp,
            'WDATE': jnow, 'CDATE': jnow, 'WTIME': tnow, 'CTIME': tnow,
            'SDATE': sdate, 'STIME': np.int32(0), 'TSTEP': np.int32(0),
            'NVARS': nv, 'VAR-LIST': ''.join([k.ljust(16) for k in outkeys]),
            'FILEDESC': 'Unknown'.ljust(80),
            'HISTORY': f'Created {now:%Y-%m-%dT%H:%M:%S}'.ljust(80),
        }

        # Overwrite any pre-existing values
        defattrs.update(self._obj.attrs)
        outf.attrs.update(defattrs)
        outf.csp.add_coords(inplace=True)

        if 'TFLAG' in self._obj.data_vars:
            tflag = self._obj['TFLAG']
            tstep = self._obj.attrs.get('TSTEP', 0)
        else:
            jday = (
                outf.coords['TSTEP'].dt.year * 1000
                + outf.coords['TSTEP'].dt.dayofyear
            )
            time = outf.coords['TSTEP'].dt.strftime('%H%M%S').astype('i')
            tflag = np.array([jday, time]).T[:, None, :].repeat(nv, 1)
            tflag = xr.DataArray(
                tflag, dims=('TSTEP', 'VAR', 'DATE-TIME',),
                attrs=dict(
                    units='<YYYYJJJ,HHMMSS>', long_name='TFLAG'.ljust(16),
                    var_desc='TFLAG'.ljust(80),
                )
            )
            dt = np.diff(outf.coords['TSTEP']).mean().astype('l') / 1e9
            dth = dt // 3600
            dtm = (dt % 3600) // 60
            dts = dt % 60
            tstep = f'{dth:.0f}{dtm:02.0f}{dts:02.0f}'

        outf['TFLAG'] = tflag
        outf.attrs['SDATE'] = np.int32(outf['TFLAG'][0, 0, 0])
        outf.attrs['STIME'] = np.int32(outf['TFLAG'][0, 0, 1])
        outf.attrs['TSTEP'] = np.int32(tstep)
        outf.attrs['WDATE'] = defattrs['WDATE']
        outf.attrs['WTIME'] = defattrs['WTIME']
        outf = outf[['TFLAG'] + outkeys]
        if reset_index:
            outf = outf.reset_index(
                ['TSTEP', 'LAY', 'ROW', 'COL'], drop=drop
            )
        return outf

    @property
    def proj4string(self):
        if not hasattr(self, '_proj4string'):
            self._proj4string = get_proj4string(self._obj.attrs)
        return self._proj4string

    @property
    def proj(self):
        if not hasattr(self, '_proj'):
            import pyproj
            self._proj = pyproj.Proj(self.proj4string)
        return self._proj

    @property
    def cno(self):
        if not hasattr(self, '_cno'):
            import pycno
            self._cno = pycno.cno(proj=self.proj)
        return self._cno

    @property
    def geodf(self):
        if not hasattr(self, '_geodf'):
            import pandas as pd
            import geopandas as gpd
            from shapely.geometry import box
            rows = self._obj['ROW'].values
            cols = self._obj['COL'].values
            midx = pd.MultiIndex.from_product([rows, cols])
            midx.names = 'ROW', 'COL'
            geoms = []
            for r in rows:
                for c in cols:
                    geoms.append(
                        box(c - 0.5, r - 0.5, c + 0.5, r + 0.5)
                    )
            self._geodf = gpd.GeoDataFrame(
                geometry=geoms, index=midx, crs=self.proj4string
            )
        return self._geodf

    @property
    def exterior(self):
        """
        Exterior polygon using row/col exterior points. This is archived
        from the dataframe once for efficiency.
        """
        if not hasattr(self, '_bbox'):
            import geopandas as gpd
            import numpy as np
            from shapely.geometry import Polygon

            nr = self._obj.attrs['NROWS']
            nc = self._obj.attrs['NCOLS']
            rowc = self._obj['ROW'].values
            colc = self._obj['COL'].values
            rows = np.append(rowc - 0.5, rowc[-1] + 0.5)
            cols = np.append(colc - 0.5, colc[-1] + 0.5)
            se = np.array([cols, cols * 0])
            ee = np.array([rows * 0 + nc + 1, rows])
            ne = np.array([cols[::-1], cols * 0 + nr + 1])
            we = np.array([rows * 0, rows[::-1]])
            points = np.concatenate([se, ee, ne, we], axis=1)
            self._bbox = gpd.GeoDataFrame(
                geometry=[Polygon(points.T)], crs=self.proj4string
            )

        return self._bbox

    def bbox(self, crs=4326):
        """
        Arguments
        ---------
        crs : scalar
            Coordinate reference system accepted by geopandas.

        Returns
        -------
        (swlon, swlat, nelon, nelat)
        """
        import numpy as np

        g = self.exterior.to_crs(crs)
        ge = g.geometry.iloc[0].envelope.exterior
        nelon, nelat = np.max(ge.xy, axis=1)
        swlon, swlat = np.min(ge.xy, axis=1)
        if self._obj.attrs['GDTYP'] == 6:
            nelat = 90
            swlon = -180
            nelon = 180
        return (swlon, swlat, nelon, nelat)

    def get_tz(self, method='longitude'):
        """
        Return the timezone.

        Arguments
        ---------
        method : str
            Currently only supports 'longitude', but more to come.

        Returns
        -------
        tz : xr.DataArray
            time zone as decimal hours since 0Z
        """
        import numpy as np
        import xarray as xr
        if method != 'longitude':
            raise ValueError('only longitude is supported at this time')
        x = self._obj['COL'].values
        y = self._obj['ROW'].values
        X, Y = np.meshgrid(x, y)
        lon, lat = self.proj(X, Y, inverse=True)
        tz = xr.DataArray((lon / 15), dims=('ROW', 'COL'))
        return tz

    def get_lst(self, times=None, method='longitude'):
        """
        Calculate LST from times in UTC using method.

        Arguments
        ---------
        times : array-like
            Times in UTC, otherwise, times will be read from TSTEP coordinate
            variable.
        method : str
            Method supported by get_tz

        Returns
        -------
        lst_hour : xr.DataArray
            Times in hours using LST offset from UTC.
        """
        import numpy as np
        import xarray as xr
        if times is None:
            times = self._obj['TSTEP']

        tz = self.get_tz()
        if isinstance(times, xr.DataArray):
            utc_hour = (
                times.dt.hour + times.dt.minute / 60. + times.dt.second / 3600.
            )
        else:
            utc_hour = np.array([
                time.hour + time.minute / 60 + time.second / 3600
                for time in times
            ])
        utc_hour = xr.DataArray(
            utc_hour,
            dims=('TSTEP',)
        )
        lst_hour = (tz + utc_hour).transpose('TSTEP', 'ROW', 'COL') % 24
        return lst_hour

    def is_overpass(self, times=None, method='longitude', satellite=None):
        """
        Identify if times UTC are consistent with an overpass of satellite.

        Arguments
        ---------
        times : array-like
            Times in UTC, otherwise, times will be read from TSTEP coordinate
            variable.
        method : str
            Method supported by get_tz
        satellite : None or str or dict
            If None, then all satellites in known_overpasstimes are processed.
            If str, the one satellite in known_overpasstimes is processed.
            If dict, then keys are satellites and values are overpass times.

        Returns
        -------
        lst_hour : xr.DataArray
            Times in hours using LST offset from UTC.
        """
        import xarray as xr
        LST = self.get_lst(times=times, method=method)
        if satellite is None:
            # Overpass in UTC space
            overpasstimes = known_overpasstimes
        else:
            if isinstance(satellite, dict):
                overpasstimes = satellite
            else:
                overpasstimes = {satellite: known_overpasstimes[satellite]}
        isoverpass = {
            key: ((LST >= (midh - 1)) & (LST <= (midh + 1)))
            for key, midh in overpasstimes.items()
        }
        return xr.Dataset(isoverpass)

    def mean_overpass(self, satellite, method='longitude', times=None):
        """
        Arguments
        ---------
        satellite: str
            Satellite to use as overpass
        method: str
            Method for converting UTC times to LST
        times: None or array-like
            Times in UTC

        Returns
        -------
        ovpf : xr.Dataset
            Average of all data at times within 1 hour plus or minus of the
            overpass time.
        """
        import xarray as xr
        inputf = self._obj
        isoverpass = self.is_overpass(
            times, method=method, satellite=satellite
        )[satellite]
        keys = [
            key for key, var in inputf.variables.items()
            if var.dims == ('TSTEP', 'LAY', 'ROW', 'COL')
        ]
        if isinstance(inputf, xr.Dataset):
            output = inputf[keys].where(isoverpass).mean('TSTEP')
            for key in keys:
                output[key].attrs.update(inputf[key].attrs)
            output.attrs.update(inputf.attrs)
        else:
            output = xr.Dataset({
                key: xr.DataArray(
                    var, dims=var.dimensions, attrs=var.getncatts()
                ).where(isoverpass).mean('TSTEP')
                for key, var in inputf.variables.items()
                if key in keys
            })
            output.attrs.update(inputf.getncatts())
        output.coords['TSTEP'] = inputf['TSTEP'].isel(TSTEP=0)
        return output

    def mole_per_m2(self, metf=None, add=True):
        """
        Arguments
        ---------
        metf : xr.Dataset
            File with ZF and DENS, or ZF, PRES, and TEMP variables.
        add : bool
            If True, add MOL_PER_M2 as a variable to self.

        Returns
        -------
        MOL_PER_M2 : xr.DataArray
        """
        import warnings
        # copied from
        # https://github.com/USEPA/CMAQ/blob/main/CCTM/src/ICL/fixed/const/
        # CONST.EXT
        R = 8.314459848
        MWAIR = 0.0289628
        if metf is None:
            metf = self._obj
        if 'ZF' in metf.variables:
            ZF = metf['ZF']
            # Layer 1 and the difference above it.
            DZ = xr.concat([
                ZF.isel(LAY=slice(None, 1)), ZF.diff('LAY', n=1)
            ], dim='LAY')
            DZ.attrs.update(ZF.attrs)
            DZ.attrs['long_name'] = 'DZ'.ljust(16)
            DZ.attrs['var_desc'] = 'diff(ZF)'.ljust(80)
            if 'DENS' in metf.variables:
                MOL_PER_M2 = metf['DENS'] / MWAIR * DZ
            elif (
                'PRES' in metf.variables and 'TEMP' in metf.variables
            ):
                P = metf['PRES']
                T = metf['TEMP']
                MOL_PER_M2 = P / R / T * DZ
                if 'Q' in metf.variables:
                    MWH2O = 0.01801528
                    # Q = gH2O/kgAir
                    Q = metf['Q']
                    MOLH2O_PER_MOLAIR = Q * 1000 / MWH2O * MWAIR
                    MOL_PER_M2 = MOL_PER_M2 * (1 - MOLH2O_PER_MOLAIR)
                else:
                    warnings.warn('Using wet mole density')
            else:
                raise KeyError(
                    'You file has ZF, but must also have DENS or PRES/TEMP'
                    + '(optionally Q)'
                )
        else:
            raise KeyError('Must have ZF and DENS or PRES/TEMP')
        MOL_PER_M2.attrs.update(dict(
            long_name='MOL_PER_M2'.ljust(16),
            var_desc='air areal density'.ljust(80),
            units='mole/m**2'.ljust(16)
        ))
        if add:
            self._obj['MOL_PER_M2'] = MOL_PER_M2

        return MOL_PER_M2

    def apply_ak(self, key, ak):
        """
        Apply averaging kernel to variable identified by key.

        Arguments
        ---------
        key : str
            Variable key to apply averaging kernel to.
        ak : xr.DataArray
            Averaging kernel.

        Returns
        -------
        out : xr.DataArray
            Output of KEY * AK).sum(LAY)
        """
        pcd = self._obj[key]
        validak = ~ak.isnull()
        validcell = validak.any('LAY')
        out = (pcd * ak).sum('LAY').where(validcell)
        out.attrs.update(pcd.attrs)
        return out
