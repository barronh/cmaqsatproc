__all__ = ['MODISL3', 'MOD04', 'MOD04_3K', 'MOD04_L2']

from ..core import satellite
from ...utils import centertobox


class MOD04(satellite):
    @property
    def short_description(self):
        desc = """OMNO2 filters valid pixels and provides weights for
destination polygon:
* valid = Land_Ocean_Quality_Flag > 1 and value == value
* pixel_area = corners from centroids +- delta
* weight = intx_fracarea
"""
        return desc

    def __init__(self, path, targetkey='Optical_Depth_Land_And_Ocean', **kwds):
        satellite.__init__(self, path, targetkey=targetkey, **kwds)

    @property
    def required_args(self):
        return (
            'Longitude', 'Latitude', 'dx', 'dy', 'Land_Ocean_Quality_Flag',
            'Optical_Depth_Land_And_Ocean'
        )

    @property
    def valid_index(self):
        if self._valididx is None:
            key = self.attrs['targetkey']
            df = self.to_dataframe(
                'Land_Ocean_Quality_Flag', key, valid_only=False
            )
            self._valididx = df.query(
                f'Land_Ocean_Quality_Flag > 1 and {key} == {key}'
            )
        return self._valididx

    @property
    def ds(self):
        """
        The ds property contains a dataset, which must contain sufficient
        data for geodf to define a polygon.
        """
        if self._ds is None:
            import xarray as xr
            ds = xr.open_dataset(self.path)
            self.ds = ds
        return self._ds

    @ds.setter
    def ds(self, ds):
        import xarray as xr
        import numpy as np

        dx = ds.Longitude.copy()
        dx[:] = np.nan
        dx.name = 'dx'
        diff_lon = np.abs(ds.Longitude.diff('Cell_Across_Swath'))
        diff_lon_mean = diff_lon.rolling(Cell_Across_Swath=2).mean().isel(
            Cell_Across_Swath=slice(1, None)
        )
        dx[:] = xr.concat([
            diff_lon.isel(
                Cell_Across_Swath=[0]
            ),
            diff_lon_mean,
            diff_lon.isel(
                Cell_Across_Swath=[-1]
            )
        ], dim='Cell_Across_Swath').values
        dy = ds.Latitude.copy()
        dy.name = 'dy'
        dy = ds.Latitude.copy()
        dy.name = 'dy'
        dy[:] = np.nan
        diff_lat = np.abs(ds.Latitude.diff('Cell_Along_Swath'))
        diff_lat_mean = diff_lat.rolling(Cell_Along_Swath=2).mean().isel(
            Cell_Along_Swath=slice(1, None)
        )
        dy[:] = xr.concat([
            diff_lat.isel(
                Cell_Along_Swath=[0]
            ),
            diff_lat_mean,
            diff_lat.isel(
                Cell_Along_Swath=[-1]
            )
        ], dim='Cell_Along_Swath').values
        ds['dx'] = dx
        ds['dy'] = dy
        self._ds = ds

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd

        if self._geodf is None:
            df = self.to_dataframe(
                'Longitude', 'Latitude', 'dx', 'dy'
            ).join(self.valid_index[[]], how='inner')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                # geometry=gpd.points_from_xy(df['Longitude'], df['Latitude']),
                geometry=df.apply(
                    lambda row: centertobox(
                        row['Longitude'], row['Latitude'], row['dx'], row['dy']
                    ),
                    axis=1
                ), crs=4326
            )
        return self._geodf

    def weights(self, othdf, option='fracarea', **kwds):
        return satellite.weights(self, othdf, option=option, **kwds)


MOD04_L2 = MOD04
MOD04_3K = MOD04


class MODISL3(satellite):
    @property
    def short_description(self):
        desc = """OMNO2 filters valid pixels and provides weights for
destination polygon:
* valid = Land_Ocean_Quality_Flag > 1 and value == value
* pixel_area = corners from centroids +- delta
* weight = intx_fracarea
"""
        return desc

    @classmethod
    def _to_xarray(cls, path):
        import pyproj
        import numpy as np
        from osgeo import gdal
        import xarray as xr
        ds = gdal.Open(path)
        sds = ds.GetSubDatasets()
        outds = xr.merge([
            cls._getda(sdn, defn)
            for sdn, defn in sds
            if (
                sdn.endswith('AOD_QA')
                or sdn.endswith('AOD_Uncertainty')
                or sdn.endswith('Optical_Depth_055')
            )
        ])
        attrs = outds.Optical_Depth_055.attrs
        nc = int(attrs['DATACOLUMNS'])
        nr = int(attrs['DATAROWS'])
        proj = pyproj.Proj(attrs['crs_wkt'])
        wx, sy = proj(
            float(attrs['WESTBOUNDINGCOORDINATE']),
            float(attrs['SOUTHBOUNDINGCOORDINATE'])
        )
        ex, ny = proj(
            float(attrs['EASTBOUNDINGCOORDINATE']),
            float(attrs['NORTHBOUNDINGCOORDINATE'])
        )
        GT = attrs['geo_transform']
        xbnds = np.linspace(wx, ex, nc + 1)
        ybnds = np.linspace(sy, ny, nr + 1)
        x = (xbnds[:-1] + xbnds[1:]) / 2
        y = (ybnds[:-1] + ybnds[1:]) / 2
        X_pixel = np.arange(nc)
        Y_line = np.arange(nr)
        x = GT[0] + X_pixel * GT[1] + Y_line * GT[2]
        y = GT[3] + X_pixel * GT[4] + Y_line * GT[5]
        outds.coords['x'] = x
        outds.coords['y'] = y
        outds['crs'] = xr.DataArray(
            0,
            dims=(),
            attrs=dict(crs_wkt=attrs['crs_wkt'])
        )
        outds['QualityLevel'] = xr.DataArray(
            (
                (outds['AOD_QA'].values.astype('<H') >> 8)
                & np.array(int('1111', 2), dtype='<H')
            ),
            dims=('nCandidate', 'y', 'x'),
            name='QualityLevel',
            attrs=dict(
                notes=(
                    "see https://amt.copernicus.org/articles/11/5741/2018/"
                    + "amt-11-5741-2018.pdf"
                ),
                long_name='QualityLevel', units='none',
                description=(
                    '0: "best", 1: "Water sediment"; 2: "n/a", 3: "one cloud",'
                    + ' 4: ">1 clouds", 5: "no retrieval",'
                    + ' 6: "no retrieval near", 7: "Climatology AOD (0.02)",'
                    + ' 8: "No retrieval glint",'
                    + ' 9: "Retrieval low due to glint", 10: "AOD near coast",'
                    + ' 11: "– Land, research quality: AOD retrieved but CM is'
                    + 'possibly cloudy"'
                )
            )
        )
        # outds.coords['x_bnds'] = xbnds
        # outds.coords['y_bnds'] = ybnds
        return outds

    @classmethod
    def _getda(self, path, defn):
        import gdal
        import xarray as xr
        import numpy as np
        ds0 = gdal.Open(path)
        data = ds0.ReadAsArray()
        attrs = ds0.GetMetadata_Dict()
        attrs['crs_wkt'] = ds0.GetProjection()
        attrs['geo_transform'] = ds0.GetGeoTransform()
        if 'valid_range' in attrs:
            vmin, vmax = eval(attrs['valid_range'])
            data = np.ma.masked_greater(np.ma.masked_less(data, vmin), vmax)
        finaldata = (
            float(attrs.get('add_offset', 0))
            + np.ma.masked_values(data, -28672)
            * float(attrs.get('scale_factor', 1))
        )
        return xr.DataArray(
            finaldata,
            dims=('nCandidate', 'y', 'x'),
            attrs=attrs,
            name=path.split(':')[-1]
        )

    @property
    def required_args(self):
        return (
            'x', 'y', 'QualityLevel', 'Optical_Depth_055', 'crs'
        )

    @classmethod
    def open_path(cls, path, **kwds):
        """
        Initialize a satellite object where the datastore will be opened
        from the path
        """
        ds = cls._to_xarray(path)
        obj = cls.from_dataset(ds, **kwds)
        return obj

    @property
    def valid_index(self):
        if self._valididx is None:
            df = self.to_dataframe('QualityLevel', 'Optical_Depth_055')
            self._valididx = df.query(
                'QualityLevel == 0 and Optical_Depth_055 == Optical_Depth_055'
            )
        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd

        try:
            crs = self.ds.crs.attrs['crs_wkt']
        except Exception:
            crs = (
                'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom'
                + ' spheroid",DATUM["Not specified (based on custom spheroid)"'
                + ',SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM['
                + '"Greenwich",0],UNIT["degree",0.0174532925199433]],'
                + 'PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center"'
                + ',0],PARAMETER["false_easting",0],PARAMETER["false_northing"'
                + ',0],UNIT["Meter",1]]'
            )
        if self._geodf is None:
            df = self.to_dataframe(
                'x', 'y',
            ).join(self.valid_index[[]], how='inner').groupby(
                ['x', 'y']
            ).first()
            df['x'] = df.index.get_level_values('x')
            df['y'] = df.index.get_level_values('y')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=df.apply(lambda row: centertobox(
                    df['x'], df['y'], width=1000, height=1000
                )),
                crs=crs
            )
        return self._geodf

    # def weights(self, *args, option='fracarea', **kwds):
    #     return satellite.weights(self, *args, option=option, **kwds)
    @classmethod
    def process(
        cls, links, grid, varkeys2d=('Optical_Depth_Land_And_Ocean',),
        varkeys3d=None, verbose=0
    ):
        return satellite.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )
