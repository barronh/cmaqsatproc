__all__ = ['MODISL3', 'MOD04', 'MOD04_3K', 'MOD04_L2', 'modis_readers']

from ..core import satellite
from ...utils import centertobox
from . import modis_readers


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

    def __init__(self, path, targetkey='Optical_Depth_Land_And_Ocean', **attrs):
        satellite.__init__(self, path, targetkey=targetkey, **attrs)

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
            inkeys = ['Land_Ocean_Quality_Flag', key]
            if self._bbox is not None:
                inkeys.extend(['Longitude', 'Latitude'])

            df = self.to_dataframe(
                *inkeys, valid_only=False
            )
            df = df.query(
                f'Land_Ocean_Quality_Flag > 1 and {key} == {key}'
            )
            if self._bbox is not None:
                swlon, swlat, nelon, nelat = self._bbox
                df = df.query(
                    f'Longitude >= {swlon} and Longitude <= {nelon}'
                    + f' and Latitude >= {swlat} and Latitude <= {nelat}'
                )
            self._valididx = df
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

    @property
    def required_args(self):
        return (
            'x', 'y', 'QualityLevel', 'Optical_Depth_055', 'crs'
        )

    @property
    def ds(self):
        """
        The ds property contains a dataset, which must contain sufficient
        data for geodf to define a polygon. The ds may be superseded by
        a df if set.
        """
        if self._ds is None:
            self._ds = modis_readers.gdalhdf4_to_xrdataset(self.path)
        return self._ds

    @property
    def valid_index(self):
        if self._valididx is None:
            inkeys = ['QualityLevel', 'Optical_Depth_055']
            if self._bbox is not None:
                inkeys.extend(['x', 'y'])
            df = self.to_dataframe(*inkeys, valid_only=False)
            df = df.query(
                'QualityLevel == 0 and Optical_Depth_055 == Optical_Depth_055'
            )
            if self._bbox is not None:
                import pyproj
                bbox = self._bbox
                proj = pyproj.Proj(self.crs)
                minx, miny = proj(bbox[0], bbox[1])
                maxx, maxy = proj(bbox[2], bbox[3])
                # xs = [minx, maxx]
                # ys = [miny, maxy]
                # minx = max(xs)
                # maxx = max(xs)
                # miny = max(ys)
                # maxy = max(ys)
                df = df.query(
                    f'x >= {minx} and x <= {maxx}'
                    + f' and y >= {miny} and y <= {maxy}'
                )
            self._valididx = df
        return self._valididx

    @property
    def crs(self):
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
        return crs

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd

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
                    row['x'], row['y'].values.T, width=1000, height=1000
                )),
                crs=self.crs
            )
        return self._geodf

    # def weights(self, *args, option='fracarea', **kwds):
    #     return satellite.weights(self, *args, option=option, **kwds)
    _varkeys2d = ('Optical_Depth_Land_And_Ocean',)
