__all__ = ['TropOMI']
from ..core import satellite
from ...utils import EasyDataFramePolygon


class TropOMI(satellite):
    @property
    def short_description(self):
        desc = """TropOMI filters valid pixels and provides weights for
destination polygon:
* valid = qa_value > 0.5
* pixel_area = not calculatedcorners from latitude_bounds, longitude_bounds
* weight = if centroid is in grid cell, 1
"""
        return desc

    @property
    def required_args(self):
        return (
            'LL_Longitude', 'LU_Longitude', 'UL_Longitude', 'UU_Longitude',
            'LL_Latitude', 'LU_Latitude', 'UL_Latitude', 'UU_Latitude',
            'PRODUCT_qa_value'
        )

    @property
    def ds(self):
        import xarray as xr
        if self._ds is None:
            self.ds = xr.open_dataset(self.path, decode_times=False)
        return self._ds

    @ds.setter
    def ds(self, ds):
        self._ds = ds
        inlatkey = 'PRODUCT_SUPPORT_DATA_GEOLOCATIONS_latitude_bounds'
        inlonkey = 'PRODUCT_SUPPORT_DATA_GEOLOCATIONS_longitude_bounds'
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
            self._valididx = self.to_dataframe(
                'PRODUCT_qa_value', valid_only=False
            ).fillna(1).query(
                'PRODUCT_qa_value > 0.5'
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
                df = self.valid_index
                self._geodf = gpd.GeoDataFrame(
                    df[[]],
                    geometry=gpd.points_from_xy(
                        df['PRODUCT_longitude'], df['PRODUCT_latitude']
                    ), crs=4326
                )

        return self._geodf

#    def weights(self, othdf, option='equal', clip=None):
#        return satellite.weights(self, othdf, option=option, clip=clip)
class TropOMICO(TropOMI):
    @classmethod
    def process(
        cls, links, grid,
        varkeys2d=('PRODUCT_carbonmonoxide_total_column',), varkeys3d=None,
        verbose=0
    ):
        return TropOMI.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )

class TropOMINO2(TropOMI):
    @classmethod
    def process(
        cls, links, grid,
        varkeys2d=('PRODUCT_nitrogendioxide_tropospheric_column',),
        varkeys3d=None, verbose=0
    ):
        return TropOMI.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )

class TropOMIHCHO(TropOMI):
    @classmethod
    def process(
        cls, links, grid,
        varkeys2d=('PRODUCT_formaldehyde_tropospheric_vertical_column',),
        varkeys3d=None, verbose=0
    ):
        return TropOMI.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )

class TropOMICH4(TropOMI):
    @classmethod
    def process(
        cls, links, grid,
        varkeys2d=('PRODUCT_methane_mixing_ratio',),
        varkeys3d=None, verbose=0
    ):
        return TropOMI.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )
