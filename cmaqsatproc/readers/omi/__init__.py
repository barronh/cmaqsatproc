__all__ = ['OMNO2', 'OMHCHO', 'OMNO2d']

from ..core import satellite
from ...utils import centertobox


class OMNO2(satellite):
    @property
    def required_args(self):
        return (
            'LL_Longitude', 'LU_Longitude', 'UL_Longitude', 'UU_Longitude',
            'LL_Latitude', 'LU_Latitude', 'UL_Latitude', 'UU_Latitude',
            'CloudFraction', 'VcdQualityFlags', 'XTrackQualityFlags',
            'ColumnAmountNO2Std'
        )

    @property
    def ds(self):
        if self._ds is None:
            import xarray as xr
            self.ds = xr.open_dataset(self.path)

        return self._ds

    @ds.setter
    def ds(self, ds):
        import warnings
        self._ds = ds
        for corneri, ckey in zip(range(4), ('LL', 'UL', 'UU', 'LU')):
            lonkey = f'{ckey}_Longitude'
            latkey = f'{ckey}_Latitude'
            if lonkey not in self._ds.data_vars:
                if 'FoV75CornerLongitude' not in self._ds.data_vars:
                    warnings.warn(
                        'Corners must be calculated and FoV75CornerLongitude'
                        + ' is missing'
                    )
                else:
                    lonx = self._ds['FoV75CornerLongitude'].isel(
                        nCorners=corneri
                    )
                    self._ds[lonkey] = lonx

            if latkey not in self._ds.data_vars:
                if 'FoV75CornerLatitude' not in self._ds.data_vars:
                    warnings.warn(
                        'Corners must be calculated and FoV75CornerLatitude'
                        + ' is missing'
                    )
                else:
                    latx = self._ds['FoV75CornerLatitude'].isel(
                        nCorners=corneri
                    )
                    self._ds[latkey] = latx

    @property
    def valid_index(self):
        if self._valididx is None:
            df = self.to_dataframe(
                'CloudFraction', 'VcdQualityFlags', 'XTrackQualityFlags',
                valid_only=False
            ).fillna(1)
            df['FinalFlag'] = (
                df['VcdQualityFlags'].astype('i') & 1
            ).astype('i')
            self._valididx = df.query(
                'FinalFlag == 0 and XTrackQualityFlags == 0'
                + ' and CloudFraction <= 0.3'
            )
        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd
        from shapely.geometry import Polygon
        if self._geodf is None:
            df = self.to_dataframe(
                'LL_Latitude', 'LL_Longitude',
                'LU_Latitude', 'LU_Longitude',
                'UL_Latitude', 'UL_Longitude',
                'UU_Latitude', 'UU_Longitude',
            ).join(self.valid_index[[]], how='inner')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=df.apply(
                    lambda row: Polygon([
                        [row['LL_Longitude'], row['LL_Latitude']],
                        [row['LU_Longitude'], row['LU_Latitude']],
                        [row['UU_Longitude'], row['UU_Latitude']],
                        [row['UL_Longitude'], row['UL_Latitude']],
                        [row['LL_Longitude'], row['LL_Latitude']],
                    ]), axis=1
                ),
                crs=4326
            )
        return self._geodf

    def weights(self, *args, **kwds):
        """
        Combines intx_fracarea with
        * area_factor = (1- (area - area_min) / area_max)
        * uncertainty_factor = 1 / ColumnAmountNO2Std

        See https://acdisc.gesdisc.eosdis.nasa.gov/data/
        Aura_OMI_Level3/OMNO2d.003/doc/README.OMNO2.pdf
        section 6.3 for more detail
        """
        wgtdf = satellite.weights(self, *args, **kwds)
        area = wgtdf['geometry'].area
        area_min = area.min()
        area_max = area.max()
        area_factor = 1 - (area - area_min) / area_max
        uncertainty_factor = 1 / self.to_dataframe('ColumnAmountNO2Std').loc[
            wgtdf.index, 'ColumnAmountNO2Std'
        ]
        wgtdf['weights'] = (
            wgtdf['intx_fracarea'] * area_factor * uncertainty_factor
        )
        return wgtdf


class OMHCHO(satellite):
    @property
    def required_args(self):
        return (
            'LL_Longitude', 'LU_Longitude', 'UL_Longitude', 'UU_Longitude',
            'LL_Latitude', 'LU_Latitude', 'UL_Latitude', 'UU_Latitude',
            'MainDataQualityFlag', 'ColumnUncertainty'
        )

    @property
    def ds(self):
        if self._ds is None:
            import xarray as xr
            self._ds = xr.open_dataset(self.path)
            slices = [('L', slice(None, -1)), ('U', slice(1, None))]
            for tkey, tslice in slices:
                for xkey, xslice in slices:
                    lonx = xr.DataArray(
                        self._ds['PixelCornerLongitudes'].isel(
                            nTimes_1=tslice, nXtrack_1=xslice
                        ).values,
                        dims=('nTimes', 'nXtrack')
                    )
                    self._ds[f'{tkey}{xkey}_Longitude'] = lonx
                    latx = xr.DataArray(
                        self._ds['PixelCornerLatitudes'].isel(
                            nTimes_1=tslice, nXtrack_1=xslice
                        ).values,
                        dims=('nTimes', 'nXtrack')
                    )
                    self._ds[f'{tkey}{xkey}_Latitude'] = latx

        return self._ds

    @property
    def valid_index(self):
        if self._valididx is None:
            df = self.to_dataframe('MainDataQualityFlag', valid_only=False)
            self._valididx = df.query('MainDataQualityFlag == 0')
        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd
        from shapely.geometry import Polygon
        if self._geodf is None:
            df = self.to_dataframe(
                'LL_Latitude', 'LL_Longitude',
                'LU_Latitude', 'LU_Longitude',
                'UL_Latitude', 'UL_Longitude',
                'UU_Latitude', 'UU_Longitude',
            ).join(self.valid_index[[]], how='inner')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=df.apply(
                    lambda row: Polygon([
                        [row['LL_Longitude'], row['LL_Latitude']],
                        [row['LU_Longitude'], row['LU_Latitude']],
                        [row['UU_Longitude'], row['UU_Latitude']],
                        [row['UL_Longitude'], row['UL_Latitude']],
                        [row['LL_Longitude'], row['LL_Latitude']],
                    ]), axis=1
                ),
                crs=4326
            )
        return self._geodf

    def weights(self, *args, **kwds):
        """
        Combines intx_fracarea with
        * area_factor = (1- (area - area_min) / area_max)
        * uncertainty_factor = 1 / ColumnUncertainty

        The area_factor is currently being used instead of the Spatial
        Response (S)

        To do, replace area_factor with spatial response
        For weights used by OMHCHOd see
        https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/
          OMHCHOd.003/doc/README_OMHCHOd_v003.pdf
        For explanation of the Spatial Response function see section 4.1 of
         https://amt.copernicus.org/articles/11/6679/2018/
        For area factor, see section 6.3 of
          https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMNO2d.003/
        doc/README.OMNO2.pdf
        """
        wgtdf = satellite.weights(self, *args, **kwds)
        if wgtdf.shape[0] == 0:
            return wgtdf
        area = wgtdf['geometry'].area
        area_min = area.min()
        area_max = area.max()
        area_factor = 1 - (area - area_min) / area_max
        uncertainty_factor = 1 / self.to_dataframe('ColumnUncertainty').loc[
            wgtdf.index, 'ColumnUncertainty'
        ]
        wgtdf['weights'] = (
            wgtdf['intx_fracarea'] * area_factor * uncertainty_factor
        )
        return wgtdf


class OMNO2d(satellite):
    def __init__(
        self, path, targetkey='ColumnAmountNO2TropCloudScreened', **kwds
    ):
        """
        Same as satellite, but takes two optional arguments.

        Arguments
        ---------
        path : str
            path to file
        targetkey : str
            key in file to use for filtering cells
        bbox : list
            Optional, bounding box to prefilter the dataset
            [swlon, swlat, nelon, nelat]
        **kwds : mappable
            Same see satellite.__init__
        """
        satellite.__init__(self, path, targetkey=targetkey, **kwds)

    @property
    def ds(self):
        if self._ds is None:
            import xarray as xr

            ds = xr.open_dataset(self.path)
            attrs = self._attrs
            if 'bbox' in attrs:
                lonslice = slice(attrs['bbox'][0], attrs['bbox'][2])
                latslice = slice(attrs['bbox'][1], attrs['bbox'][3])
                self._ds = ds.sel(lon=lonslice, lat=latslice)
            else:
                self._ds = ds

        return self._ds

    @property
    def valid_index(self):
        if self._valididx is None:
            key = self.attrs['targetkey']
            df = self.to_dataframe(key, valid_only=False)
            self._valididx = df.query(f'{key} == {key}')
        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd
        import numpy as np
        
        if self._geodf is None:
            # Updated to dataframe method. Needs checking
            lldf = self.to_dataframe('lon', 'lat', valid_only=False)
            dlon = np.diff(lldf['lon'].values).mean()
            dlat = np.diff(lldf['lat']).mean()
            df = self.to_dataframe(
                'lon', 'lat',
            ).join(self.valid_index[[]], how='inner')
            df['lon'] = df.index.get_level_values('lon')
            df['lat'] = df.index.get_level_values('lat')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=df.apply(
                    lambda row: centertobox(
                        row['lon'], row['lat'], width=dlon, height=dlat
                    ), axis=1
                ),
                crs=4326
            )
        return self._geodf

    def weights(self, *args, **kwds):
        """
        Combines intx_fracarea with Weight from dataset.
        See https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/
        OMNO2d.003/doc/README.OMNO2.pdf
        section 6.3 for more detail
        """
        wgtdf = satellite.weights(self, *args, **kwds)
        wgtdf['weights'] = (
            wgtdf['intx_fracarea']
            * self.to_dataframe('Weight').loc[wgtdf.index, 'Weight']
        )
        return wgtdf
