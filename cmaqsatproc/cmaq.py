__all__ = ['CMAQGrid']

from .utils import centertobox

known_overpasstimes = dict(
    aura=13.75,
    terra=10.5,
    aqua=13.5,
    metop_am=9.5,
    metop_pm=21.5,
)


def _default_griddesc(GDNAM):
    import PseudoNetCDF as pnc
    import tempfile

    with tempfile.NamedTemporaryFile() as gdfile:
        gdfile.write(b"""' '
'LATLON'
  1  0.0 0.0 0.0 0.0 0.0
'POLSTE_HEMI'
  6         1.000        45.000       -98.000       -98.000        90.000
'LamCon_40N_97W'
  2        33.000        45.000       -97.000       -97.000        40.000
' '
'US_1deg'
'LATLON'               -140.0        20.0      1.0      1.0   90   40 1
'US_0pt1deg'
'LATLON'               -140.0        20.0      0.1      0.1  900  400 1
'global_1deg'
'LATLON'               -180.0       -90.0      1.0      1.0  360  180 1
'global_0pt1deg'
'LATLON'               -180.0       -90.0      0.1      0.1 3600 1800 1
'global_2x2.5'
'LATLON'               -181.5       -89.0      2.5      2.0  144   89 1
'global_4x5'
'LATLON'               -182.5       -88.0      5.0      4.0   72   44 1
'108NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 108000.0 108000.0  187  187 1
'324NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 324000.0 324000.0   63   63 1
'1188NHEMI2'
'POLSTE_HEMI'     -10098000.0 -10098000.0 1188000. 1188000.   17   17 1
'1188US1'
'LamCon_40N_97W'   -2556000.0  -1728000.0 1188000. 1188000.    6    4 1
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
'LamCon_40N_97W'   -2556000.0  -1728000.0   4000.0   4000.0 5508 3588 1
'12US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0  12000.0  12000.0  396  246 1
'4US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0  12000.0  12000.0 1188  738 1
'1US2'
'LamCon_40N_97W'   -2412000.0  -1620000.0  12000.0  12000.0 4752 2952 1
'36US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0  36000.0  36000.0  172  148 1
'108US3'
'LamCon_40N_97W'   -2952000.0  -2772000.0 108000.0 108000.0   60   50 1
' '""")
        gdfile.flush()
        outf = pnc.pncopen(gdfile.name, format='griddesc', GDNAM=GDNAM)
    return outf


class CMAQGrid:
    @classmethod
    def from_gf(cls, gf):
        out = cls.__new__(cls)
        out.GDNAM = gf.GDNAM.strip()
        out.gf = gf
        out.proj4string = out.gf.getproj(withgrid=True, projformat='proj4')
        return out

    def __init__(self, GDNAM, gdpath=None):
        """
        Arguments
        ---------
        GDNAM : str
            Name of domain
        gdpath : str or None
            Path to GRIDDESC file. If None, then a default is provided with
            access to typical EPA domains (12US1, 12US2, 36US3, 108NHEMI2) and
            a few test domains (1188NHEMI2, 108US1).
        """
        import PseudoNetCDF as pnc
        import warnings
        with warnings.catch_warnings(record=True):
            if gdpath is None:
                gf = _default_griddesc(GDNAM)
                gf.HISTORY = 'From GRIDDESC'
            else:
                gf = pnc.pncopen(
                    gdpath, format='griddesc', GDNAM=GDNAM
                )
        self.GDNAM = GDNAM
        self.gf = gf.subset(['DUMMY'])
        self.gf.SDATE = 1970001
        self.gf.updatetflag(overwrite=True)
        self.proj4string = self.gf.getproj(withgrid=True, projformat='proj4')

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
            self._cno = pycno.cno(self.proj)
        return self._cno
            
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

            nr = self.gf.NROWS
            nc = self.gf.NCOLS
            rows = np.arange(nr + 1)
            cols = np.arange(nc + 1)
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
        if self.gf.GDTYP == 6:
            nelat = 90
            swlon = -180
            nelon = 180
        return (swlon, swlat, nelon, nelat)

    def to_ioapi(
        self, df, TSTEP=0, LAY=0, ROW='ROW', COL='COL', rename=None, **kwds
    ):
        """
        Create an IOAPI file from a dataframe.

        Arguments
        ---------
        df : pandas.Dataframe
            Columns represent output data. attrs of the dataframe and columns
            will be used to construct metadata for the ioapi file.
        TSTEP : int
        LAY : int
        ROW : str
            Index name to use for IOAPI rows
        COL : str
            Index name to use for IOAPI cols
        rename : mappable
            Optional dictionary to rename df columns to make IOAPI variable
            names. Useful for making sure columns are less than 16 characters.
        kwds : mappable
            Keywords passed through to the template function.

        Returns
        -------
        outf : PseudoNetCDF.cmaqfiles.ioapi_base
            Output file
        """
        import numpy as np
        if rename is None:
            rename = {key: key for key in df.columns}
        elif callable(rename):
            rename = {key: rename(key) for key in df.columns}

        if isinstance(TSTEP, str):
            TSTEP = df.index.get_level_values(TSTEP)
        if isinstance(LAY, str):
            LAY = df.index.get_level_values(LAY)
            vglvls_mid, layer_idx = np.unique(LAY, return_inverse=True)
            vglvls_edge = np.interp(
                np.arange(vglvls_mid.size + 1) - 0.5,
                np.arange(vglvls_mid.size),
                vglvls_mid
            )
            LAY = layer_idx
        else:
            vglvls_edge = np.asarray([1, 0])
        if isinstance(ROW, str):
            ROW = df.index.get_level_values(ROW)
        if isinstance(COL, str):
            COL = df.index.get_level_values(COL)

        if 'vglvls' not in kwds:
            kwds['vglvls'] = vglvls_edge

        varkeys = sorted(rename)
        outkeys = [rename[key] for key in varkeys]
        outf = self.template(outkeys, **kwds)
        outf.setncatts(df.attrs)
        for varkey in varkeys:
            outkey = rename.get(varkey, varkey)
            varo = outf.variables[outkey]
            ds = df[varkey]
            varo[TSTEP, LAY, ROW, COL] = ds.values
            varo.setncatts(ds.attrs)

        return outf

    def template(self, varkeys, vglvls=(1, 0), vgtop=5000, ntsteps=1, **kwds):
        """
        Create an empty IOAPI file with variables that have arbitrary values

        Arguments
        ---------
        varkeys : list
            List of strings to make as variables.
        vglvls : iterable
            Iterable of IOAPI level edges
        vgtop : scalar
            Top of the IOAPI model
        ntsteps : int
            Number of output time steps

        Returns
        -------
        outf : PseudoNetCDF.cmaqfiles.ioapi_base
            NetCDF-like file with arbitrary values for all varkeys
        """
        import numpy as np

        outf = self.gf.renameVariable('DUMMY', varkeys[0]).mask(values=0)
        nz = len(vglvls) - 1
        if nz > 1:
            outf = outf.slice(LAY=[0] * nz)
        if ntsteps > 1:
            outf = outf.slice(TSTEP=[0] * nz)
        outf.VGLVLS = np.asarray(vglvls, dtype='f')
        outf.VGTOP = np.float32(vgtop)
        outf.setncatts(kwds)
        del outf.dimensions['nv']
        del outf.Conventions
        tmpvar = outf.variables[varkeys[0]]
        del tmpvar.coordinates
        tmpvar.long_name = varkeys[0].ljust(16)
        tmpvar.var_desc = varkeys[0].ljust(80)
        tmpvar.units = 'unknown'
        for varkey in varkeys[1:]:
            newvar = outf.copyVariable(outf.variables[varkeys[0]], key=varkey)
            newvar.long_name = varkey.ljust(16)
            newvar.var_desc = varkey.ljust(80)

        outf.updatemeta()
        outf.updatetflag(overwrite=True)
        outf.variables.move_to_end('TFLAG', last=False)
        return outf

    @property
    def geodf(self):
        if not hasattr(self, '_geodf'):
            import pandas as pd
            import geopandas as gpd
            import numpy as np
            nr = self.gf.NROWS
            nc = self.gf.NCOLS
            rows = np.arange(nr)
            cols = np.arange(nc)
            midx = pd.MultiIndex.from_product([rows, cols])
            midx.names = 'ROW', 'COL'
            geoms = []
            for r in rows:
                for c in cols:
                    geoms.append(
                        centertobox(xc=c + 0.5, yc=r + 0.5, width=1, height=1)
                    )
            self._geodf = gpd.GeoDataFrame(
                geometry=geoms, index=midx, crs=self.proj4string
            )
        return self._geodf

    def get_tz(self, method='longitude'):
        import numpy as np
        import xarray as xr
        if method != 'longitude':
            raise ValueError('only longitude is supported at this time')
        i = np.arange(self.gf.NCOLS)
        j = np.arange(self.gf.NROWS)
        I, J = np.meshgrid(i, j)
        lon, lat = self.gf.ij2ll(I, J)
        tz = xr.DataArray((lon / 15), dims=('ROW', 'COL'))
        return tz

    def get_lst(self, times, method='longitude'):
        import numpy as np
        import xarray as xr
        tz = self.get_tz()
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

    def is_overpass(self, times, method='longitude', satellite=None):
        import xarray as xr
        LST = self.get_lst(times, method=method)
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

    def mean_overpass(self, inputf, satellite, method='longitude', times=None):
        import xarray as xr
        if times is None:
            try:
                times = inputf.getTimes()
            except Exception:
                times = inputf['TSTEP']
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
        return output

    @classmethod
    def mole_per_m2(self, metf, add=True):
        import warnings
        # copied from
        # https://github.com/USEPA/CMAQ/blob/main/CCTM/src/ICL/fixed/const/
        # CONST.EXT
        R = 8.314459848
        MWAIR = 0.0289628
        if 'ZF' in metf.variables:
            DZ = metf['ZF'].copy()
            DZ[1:] = DZ[1:] - metf['ZF'][:-1].values
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
        if add:
            metf['MOL_PER_M2'] = MOL_PER_M2

        return MOL_PER_M2
