__all__ = ['CMAQGrid']

from .utils import centertobox


class CMAQGrid:
    def __init__(self, gdpath, GDNAM):
        import PseudoNetCDF as pnc
        self._gf = pnc.pncopen(
            gdpath, format='griddesc', GDNAM=GDNAM
        ).subset(['DUMMY'])
        self._gf.SDATE = 1970001
        self._gf.updatetflag(overwrite=True)
        self.proj = self._gf.getproj(withgrid=True, projformat='proj4')
        self._geodf = None
        self._bbox = None

    @property
    def exterior(self):
        """
        Exterior polygon using row/col exterior points. This is archived
        from the dataframe once for efficiency.
        """
        if self._bbox is None:
            import geopandas as gpd
            import numpy as np
            from shapely.geometry import Polygon

            nr = self._gf.NROWS
            nc = self._gf.NCOLS
            rows = np.arange(nr + 1)
            cols = np.arange(nc + 1)
            se = np.array([cols, cols * 0])
            ee = np.array([rows * 0 + nc + 1, rows])
            ne = np.array([cols[::-1], cols * 0 + nr + 1])
            we = np.array([rows * 0, rows[::-1]])
            points = np.concatenate([se, ee, ne, we], axis=1)
            self._bbox = gpd.GeoDataFrame(
                geometry=[Polygon(points.T)], crs=self.proj
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
        if self._gf.GDTYP == 6:
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
        if rename is None:
            rename = {key: key for key in df.columns}
        if isinstance(TSTEP, str):
            TSTEP = df.index.get_level_values(TSTEP)
        if isinstance(LAY, str):
            LAY = df.index.get_level_values(LAY)
        if isinstance(ROW, str):
            ROW = df.index.get_level_values(ROW)
        if isinstance(COL, str):
            COL = df.index.get_level_values(COL)

        if 'vglvls' not in kwds:
            import numpy as np
            vglvls_mid = np.unique(np.asarray(LAY))
            xp = np.arange(vglvls_mid.size)
            x = np.arange(xp.size + 1) - 0.5
            kwds['vglvls'] = np.interp(x, xp, vglvls_mid)

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
        Create an empty IOAPI file with masked variables

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
            NetCDF-like file with masked valeus for all varkeys
        """
        import numpy as np

        outf = self._gf.renameVariable('DUMMY', varkeys[0]).mask(values=0)
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
        if self._geodf is None:
            import pandas as pd
            import geopandas as gpd
            import numpy as np
            nr = self._gf.NROWS
            nc = self._gf.NCOLS
            rows = np.arange(nr)
            cols = np.arange(nc)
            midx = pd.MultiIndex.from_product([rows, cols])
            midx.names = 'ROW', 'COL'
            geoms = []
            for r in rows:
                for c in cols:
                    geoms.append(
                        centertobox(xc=c+0.5, yc=r+0.5, width=1, height=1)
                    )
            self._geodf = gpd.GeoDataFrame(
                geometry=geoms, index=midx, crs=self.proj
            )
        return self._geodf
