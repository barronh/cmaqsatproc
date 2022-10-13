__all__ = ['IASI_NH3']

from ..core import satellite
from ...utils import EasyDataFramePoint


class IASI_NH3(satellite):
    __doc__ = """IASI Ammonia Processor
    * bbox subsets each pixel
    * valid = (nh3_total_column > -999) & (AMPM == 0)
    * geometry is a 0.1 degree buffer
    * pixel_area = diameter assumed to be 0.1 degree, which is approximately
      the 12 km discussed in Clerbaux (doi: 10.5194/acp-9-6041-2009).
    * weights = area
    """
    _defaultkeys = ('nh3_total_column',)
    _geokeys = ('latitude', 'longitude')

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a IASI_NH3 file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: IASI_NH3
            Satellite processing instance
        """
        import xarray as xr
        import numpy as np

        ds = xr.open_dataset(path, decode_cf=False, **kwargs)
        ds['valid'] = (
            (ds['nh3_total_column'] > -999)
            & (ds['AMPM'] == 0)
        )
        ds.coords['time'] = np.arange(ds.dims['time'])
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            df = ds[['latitude', 'longitude']].to_dataframe()
            df = df.query(
                f'longitude >= {swlon} and longitude <= {nelon}'
                + f' and latitude >= {swlat} and latitude <= {nelat}'
            )
            times = df.index.get_level_values('time').values
            ds = ds.sel(time=times)
            ds['valid'] = ds['valid'] & (
                (ds['latitude'] >= swlat) & (ds['latitude'] <= nelat)
                & (ds['longitude'] >= swlon) & (ds['longitude'] <= nelon)
            )

        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        sat = cls()
        sat.path = path
        sat.bbox = bbox
        sat.ds = ds
        return sat

    def to_dataframe(
        self, *varkeys, valid=True, geo=False, default_keys=False
    ):
        """
        Arguments
        ---------
        varkeys : iterable
            Keys to use in the dataframe. Defaults to class._defaultkeys that
            is class specific
        valid : bool
            If true, only return valid pixels.
        geo : bool
            Add geometry to output a geopandas.GeoDataFrame

        Returns
        -------
        df : pandas.DataFrame or geopandas.GeoDataFrame
        """
        df = satellite.to_dataframe(
            self, *varkeys, valid=valid, geo=False, default_keys=default_keys
        )
        if not (geo is False):
            import geopandas as gpd
            if geo is True:
                geo = EasyDataFramePoint

            gdf = self.to_dataframe(*self._geokeys, valid=valid)
            gdf['x'] = gdf['longitude']
            gdf['y'] = gdf['latitude']
            gdf = gpd.GeoDataFrame(
                gdf, geometry=geo(gdf).buffer(0.1), crs=4326
            )
            df = gdf[['geometry']].join(df)
        return df

    @classmethod
    def cmaq_process(
        cls, qf, l3, satellite='metop_am'
    ):
        """
        Process CMAQ as though it were observed by IASI, which is simply based
        on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., NH3)
        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        satellite : str
            Key for satellite overpass metop_am or metop_pm.

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        saodkey = 'nh3_total_column'
        overf = qf.csp.mean_overpass(satellite=satellite).where(
            ~l3[saodkey].isnull()
        )
        return overf
