__all__ = ['IASI_NH3']

from ..core import satellite
from ...utils import centertobox, EasyDataFramePolygon


class IASI_NH3(satellite):
    @property
    def short_description(self):
        desc = """IASI_NH3 filters valid pixels and provides weights for
destination polygon:
* valid = Not missing, AMPM == 0 (AM)
* pixel_area = diameter assumed to be 0.1 degree, which is approximately
  the 12 km discussed in Clerbaux (doi: 10.5194/acp-9-6041-2009).
* weights = area / uncertainty
"""
        return desc

    @property
    def required_args(self):
        return (
            'longitude', 'latitude', 'nh3_total_column',
            'nh3_total_column_uncertainty'
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

    def to_dataframe(self, *keys, valid_only=True):
        """
        Export keys to a dataframe for only valid values.

        Arguments
        ---------
        keys : list
            List of keys to pull from the self.ds or self.df for all valid
            pixels.
        valid_only : bool
            Default True. Only export a dataframe where valid_idx.

        Returns
        -------
        df : pandas.Dataframe
            Dataframe with columns keys
        """
        # Mostly copied from core, but the dataset input is not CF compliant
        # There is no valid coordinate. Theoretically, orbit/scanline/pixel
        # should be uniquely identifying. The file, however, is setup with
        # a time coordiante...
        #
        # * time is duplicated many times
        # * orbit_number has invalid values, which prevents us from using
        #   it as an index
        props = {}
        if self.df is None:
            df = self.ds[list(keys)].to_dataframe().reset_index()[list(keys)]
            for key in keys:
                if key in self.ds.data_vars:
                    props[key] = self.ds[key].attrs
        else:
            df = self.df[list(keys)]
            for key in keys:
                props[key] = self.df[key].attrs

        if valid_only:
            outdf = df.join(
                self.valid_index[[]], how='inner'
            )
        else:
            outdf = df

        for key in keys:
            if key in props:
                outdf[key].attrs.update(props[key])

        return outdf

    @property
    def valid_index(self):
        if self._valididx is None:
            df = self.to_dataframe(
                'nh3_total_column', 'AMPM',
                valid_only=False
            )
            self._valididx = df.query(
                'nh3_total_column == nh3_total_column'
                + ' and AMPM == 0'
            )
        return self._valididx

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry
        """
        import geopandas as gpd
        if self._geodf is None:
            df = self.to_dataframe(
                'latitude', 'longitude'
            ).join(self.valid_index[[]], how='inner')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=gpd.points_from_xy(
                    df['longitude'], df['latitude']
                ).buffer(0.1),
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
        # Currently thowing an warning due to calculation in spherical space.
        # Using area in a relative form, so for now that is okay.
        area = wgtdf['geometry'].area
        area_factor = area
        uncertainty_factor = 1 / self.to_dataframe(
            'nh3_total_column_uncertainty'
        ).loc[
            wgtdf.index, 'nh3_total_column_uncertainty'
        ]
        wgtdf['weights'] = (
            wgtdf['intx_fracarea'] * area_factor * uncertainty_factor
        )
        return wgtdf

    @classmethod
    def process(
        cls, links, grid,
        varkeys2d=('nh3_total_column', 'nh3_total_column_uncertainty',),
        varkeys3d=None, verbose=0
    ):
        return satellite.process(
            links, grid, varkeys2d=varkeys2d, varkeys3d=varkeys3d,
            verbose=verbose
        )
