__all__ = ['TropOMI']
from ..core import satellite
from ...utils import EasyDataFramePolygon


class TropOMI(satellite):
    @property
    def short_description(self):
        desc = """TropOMI filters valid pixels and provides weights for
destination polygon:
* valid = qa_value > 0.5
* pixel_area = corners from latitude_bounds, longitude_bounds or interpolation
* weight = fractional area
"""
        return desc

    @property
    def required_args(self):
        return (
            'LL_Longitude', 'LU_Longitude', 'UL_Longitude', 'UU_Longitude',
            'LL_Latitude', 'LU_Latitude', 'UL_Latitude', 'UU_Latitude',
            'qa_value'
        )

    @property
    def ds(self):
        import xarray as xr
        if self._ds is None:
            ds = xr.open_dataset(self.path, decode_times=False)
            rename = {
                k: k.replace('PRODUCT_', '')
                for k in ds.variables
            }
            dimrename = {
                k: k.replace('PRODUCT_', '')
                for k in tmp.ds.dims
            }
            rename.update(dimrename)
            self._ds = ds.rename(rename)
        return self._ds

    @ds.setter
    def ds(self, ds):
        self._ds = ds
        inlatkey = 'SUPPORT_DATA_GEOLOCATIONS_latitude_bounds'
        inlonkey = 'SUPPORT_DATA_GEOLOCATIONS_longitude_bounds'
        if not (
            inlatkey in self._ds.variables
            and inlonkey in self._ds.variables
        ) and self._attrs.get('usepoly', True):
            import warnings
            import xarray as xr
            import numpy as np
            warnings.warn(f'Missing {inlatkey} or {inlonkey} - approximating')
            ns = self._ds.dims['scanline']
            ng = self._ds.dims['ground_pixel']
            lonc = self._ds['longitude'].interp(
                scanline=xr.DataArray(
                    np.arange(ns + 1) - 0.5, dims=('scanline_edge',)
                ),
                ground_pixel=xr.DataArray(
                    np.arange(ng + 1) - 0.5, dims=('ground_pixel_edge',)
                )
            ).values
            latc = self._ds['latitude'].interp(
                scanline=xr.DataArray(
                    np.arange(ns + 1) - 0.5, dims=('scanline_edge',)
                ),
                ground_pixel=xr.DataArray(
                    np.arange(ng + 1) - 0.5, dims=('ground_pixel_edge',)
                )
            ).values
            self._ds[inlatkey] = xr.DataArray(
                np.stack([
                    latc[:, :-1, :-1], latc[:, 1:, :-1],
                    latc[:, 1:, 1:], latc[:, :-1, 1:]
                ], axis=3),
                dims=(
                    'time', 'scanline', 'ground_pixel',
                    'corner'
                )
            )
            self._ds[inlonkey] = xr.DataArray(
                np.stack([
                    lonc[:, :-1, :-1], lonc[:, 1:, :-1],
                    lonc[:, 1:, 1:], lonc[:, :-1, 1:]
                ], axis=3),
                dims=(
                    'time', 'scanline', 'ground_pixel',
                    'corner'
                )
            )
        if self._attrs.get('usepoly', True):
            for idx, key in [(0, 'LL'), (1, 'UL'), (2, 'UU'), (3, 'LU')]:
                latkey = f'{key}_Latitude'
                if latkey not in ds.data_vars:
                    self._ds[latkey] = self._ds[inlatkey][..., idx]
                lonkey = f'{key}_Longitude'
                if lonkey not in ds.data_vars:
                    self._ds[lonkey] = self._ds[inlonkey][..., idx]

    @property
    def valid_index(self):
        if self._valididx is None:
            if self._bbox is None:
                self._valididx = self.to_dataframe(
                    'qa_value', valid_only=False
                ).fillna(1).query(
                    'qa_value > 0.5'
                )
            else:
                swlon, swlat, nelon, nelat = self._bbox
                df = self.to_dataframe(
                    'qa_value', 'longitude',
                    'latitude', valid_only=False
                ).fillna(1).query('qa_value > 0.5')
                self._valididx = df.query(
                    f'longitude >= {swlon}'
                    + f' and longitude <= {nelon}'
                    + f' and latitude >= {swlat}'
                    + f' and latitude <= {nelat}'
                )

        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd

        if self._geodf is None:
            # Use polygons soon. Very slow at present.
            if self._attrs.get('usepoly', True):
                df = self.to_dataframe(
                    'LL_Latitude', 'LL_Longitude',
                    'LU_Latitude', 'LU_Longitude',
                    'UL_Latitude', 'UL_Longitude',
                    'UU_Latitude', 'UU_Longitude',
                ).join(self.valid_index[[]], how='inner')
                self._geodf = gpd.GeoDataFrame(
                    df[[]],
                    geometry=EasyDataFramePolygon(df),
                    crs=4326
                )
            else:
                df = self.to_dataframe(
                    'longitude', 'latitude',
                ).join(self.valid_index[[]], how='inner')
                self._geodf = gpd.GeoDataFrame(
                    df[[]],
                    geometry=gpd.points_from_xy(
                        df['longitude'], df['latitude']
                    ), crs=4326
                )

        return self._geodf

    _zkey = 'layer'

# If usepoly is False, then use equal weighting...
#    def weights(self, othdf, option='equal', clip=None):
#        return satellite.weights(self, othdf, option=option, clip=clip)


class TropOMICO(TropOMI):
    _varkeys2d = ('carbonmonoxide_total_column',)
    _varkeys3d = (
        'SUPPORT_DATA_DETAILED_RESULTS_column_averaging_kernel',
    )


class TropOMINO2(TropOMI):
    _varkeys2d = (
        'nitrogendioxide_tropospheric_column',
        'air_mass_factor_troposphere'
    )
    _varkeys3d = ('averaging_kernel',)


class TropOMIHCHO(TropOMI):
    _varkeys2d = (
        'formaldehyde_tropospheric_vertical_column',
        (
            'SUPPORT_DATA_DETAILED_RESULTS_'
            + 'formaldehyde_tropospheric_air_mass_factor'
        )
    )
    _varkeys3d = ('averaging_kernel',)


class TropOMICH4(TropOMI):
    _varkeys2d = ('methane_mixing_ratio',)
    _varkeys3d = ('averaging_kernel',)
