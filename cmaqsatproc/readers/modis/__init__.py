__all__ = ['MODISL3', 'MOD04', 'MOD04_3K', 'MOD04_L2', 'modis_readers']

from ..core import satellite
from ...utils import EasyDataFramePoint, grouped_weighted_avg
from . import modis_readers


class MOD04(satellite):
    __doc__ = """
    MOD04 (MOD04_3K or MOD04_L2) filters valid pixels and provides weights for
    destination polygon:
    * valid = Land_Ocean_Quality_Flag > 1 and value == value
    * pixel_area = corners from centroids +- delta
    * weight = intx_fracarea
    """
    _defaultkeys = ('Optical_Depth_Land_And_Ocean',)

    @classmethod
    def shorten_name(cls, key):
        """
        Provide a short name for long keys. This is useful for renaming
        variables to fit IOAPI 16 character restrictions.

        Arguments
        ---------
        key : str
            Original variable name.

        Returns
        -------
        shortkey : str
            Shortened key
        """
        key = key.replace('Optical_Depth_Land_And_Ocean', 'LAND_OCEAN_AOD')
        key = key.replace('Land_Ocean_Quality_Flag', 'LAND_OCEAN_QA')
        return key

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Opens dataset and defines valid where the pixel is within bbox,
        flags are appropriate, and inteprolates pixel centers to corners.

        Arguments
        ---------
        path : str
            Path to a MOD04 OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: MOD04
            Satellite processing instance
        """
        import xarray as xr
        import numpy as np
        ds = xr.open_dataset(path, **kwargs)
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            df = ds[['Longitude', 'Latitude']].to_dataframe()
            df = df.query(
                f'Longitude >= {swlon} and Longitude <= {nelon}'
                + f' and Latitude >= {swlat} and Latitude <= {nelat}'
            )
            cas = df.index.get_level_values('Cell_Along_Swath').values
            ds = ds.sel(Cell_Along_Swath=slice(cas.min(), cas.max() + 1))

        Cell_Along_Swath_Edges = np.concatenate(
            [
                ds.Cell_Along_Swath[:1],
                (
                    ds.Cell_Along_Swath[1:].values
                    + ds.Cell_Along_Swath[:-1].values
                ) / 2,
                ds.Cell_Along_Swath[-1:]
            ], axis=0
        )
        Cell_Across_Swath_Edges = np.concatenate(
            [
                ds.Cell_Across_Swath[:1],
                (
                    ds.Cell_Across_Swath[1:].values
                    + ds.Cell_Across_Swath[:-1].values
                ) / 2,
                ds.Cell_Across_Swath[-1:]
            ], axis=0
        )
        lat_edges = ds.Latitude.interp(
            Cell_Along_Swath=Cell_Along_Swath_Edges,
            Cell_Across_Swath=Cell_Across_Swath_Edges
        )
        lon_edges = ds.Longitude.interp(
            Cell_Along_Swath=Cell_Along_Swath_Edges,
            Cell_Across_Swath=Cell_Across_Swath_Edges
        )
        corner_slices = {
            'll': (slice(None, -1), slice(None, -1)),
            'lu': (slice(None, -1), slice(1, None)),
            'ul': (slice(1, None), slice(None, -1)),
            'uu': (slice(1, None), slice(1, None)),
        }
        for corner, (xslice, aslice) in corner_slices.items():
            ds[f'{corner}_x'] = xr.DataArray(
                lon_edges.isel(
                    Cell_Along_Swath=aslice, Cell_Across_Swath=xslice
                ),
                dims=ds.Longitude.dims,
                coords=ds.Longitude.coords,
            )
            ds[f'{corner}_y'] = xr.DataArray(
                lat_edges.isel(
                    Cell_Along_Swath=aslice, Cell_Across_Swath=xslice
                ),
                dims=ds.Longitude.dims,
                coords=ds.Longitude.coords,
            )
        ds['valid'] = (
            (ds['Land_Ocean_Quality_Flag'] > 1)
            & (~np.isnan(ds['Optical_Depth_Land_And_Ocean']))
        )
        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat


class MOD04_L2(MOD04):
    __doc__ = """
    MOD04_L3 filters valid pixels and provides weights for destination polygon:
    * valid = Land_Ocean_Quality_Flag > 1 and value == value
    * pixel_area = corners from centroids +- delta
    * weight = intx_fracarea
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "MOD04_L2".

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download or OpenDAP
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'MOD04_L2')
        return MOD04.cmr_links(method=method, **kwargs)


class MOD04_3K(MOD04):
    __doc__ = """
    MOD04_3K filters valid pixels and provides weights for destination polygon:
    * valid = Land_Ocean_Quality_Flag > 1 and value == value
    * pixel_area = corners from centroids +- delta
    * weight = intx_fracarea
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "MOD04_3K".

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download or OpenDAP
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'MOD04_3K')
        return MOD04.cmr_links(method=method, **kwargs)


class MODISL3(satellite):
    __doc__ = """
    MODISL3 filters valid pixels and provides weights for
    destination polygon:
    * valid = Land_Ocean_Quality_Flag > 1 and value == value
    * pixel_area = corners from centroids +- delta
    * weight = intx_fracarea
    """
    _defaultkeys = ('Optical_Depth_055',)
    _geokeys = ('x', 'y')

    def to_dataframe(self, *varkeys, valid=True, geo=False, default_keys=False):
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

            gdf = self.to_dataframe(*self._geokeys)
            gdf['x'] = gdf.index.get_level_values('x')
            gdf['y'] = gdf.index.get_level_values('y')
            gdf = gpd.GeoDataFrame(
                gdf, geometry=geo(gdf), crs=self.ds.crs.crs_wkt
            )
            df = gdf[['geometry']].join(df)
        return df

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Opens dataset (remote or local). Defines valid as QualityLevel == 0.
        Interpolates pixel centers to corners. And adjusts valid to only
        include pixels in bbox.

        Arguments
        ---------
        path : str
            Path to a MOD04 OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: MODISL3
            Satellite processing instance
        """
        import pyproj
        import xarray as xr
        ds = modis_readers.gdalhdf4_to_xrdataset(path, **kwargs)
        ds['valid'] = ds['QualityLevel'] == 0
        sat = cls()
        sat.proj = pyproj.Proj(ds.crs.crs_wkt)
        corner_slices = {
            'll': (slice(None, -1), slice(None, -1)),
            'ul': (slice(1, None), slice(None, -1)),
            'uu': (slice(1, None), slice(1, None)),
            'lu': (slice(None, -1), slice(1, None)),
        }
        xbounds, ybounds = xr.broadcast(
            ds['x_bounds'], ds['y_bounds']
        )
        for corner, (xslice, yslice) in corner_slices.items():
            ds[f'{corner}_x'] = xr.DataArray(
                xbounds.isel(x_bounds=xslice, y_bounds=yslice),
                dims=('x', 'y'),
                coords=dict(
                    x=ds.coords['x'],
                    y=ds.coords['y'],
                )
            )
            ds[f'{corner}_y'] = xr.DataArray(
                ybounds.isel(x_bounds=xslice, y_bounds=yslice),
                dims=('x', 'y'),
                coords=dict(
                    x=ds.coords['x'],
                    y=ds.coords['y'],
                )
            )
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            swx, swy = sat.proj(swlon, swlat)
            nex, ney = sat.proj(nelon, nelat)
            ds['valid'] = (
                ds['valid']
                & (ds['x'] > swx) & (ds['x'] < nex)
                & (ds['y'] > swy) & (ds['y'] < ney)
            )

        sat.ds = ds
        sat.bbox = bbox
        return sat

    def to_level3(
        self, *varkeys, grid, griddims=None, weighting='equal',
        as_dataset=True, verbose=0
    ):
        """
        Arguments
        ---------
        varkeys : iterable
            See to_dataframe
        grid : geopandas.GeoDataFrame or xr.Dataset
            Defines the grid used as the L3 destination
        griddims : iterable
            Defaults to grid.index.names
        weighting : str
            Passed as option to self.add_weights

        Returns
        -------
        outputs : dict
            Dictionary of outputs by output dimensions
        """
        import xarray as xr

        if isinstance(grid, xr.Dataset) and weighting == 'equal':
            df = self.to_dataframe('Optical_Depth_055')
            lon, lat = self.proj(
                df.index.get_level_values('x').values,
                df.index.get_level_values('y').values,
                inverse=True
            )
            i, j = grid.csp.proj(lon, lat)
            # find nearest row/col
            df['ROW'] = grid['ROW'].sel(ROW=j, method='nearest')
            df['COL'] = grid['COL'].sel(COL=i, method='nearest')
            df['weight'] = 1
            df.set_index(['ROW', 'COL'], append=True, inplace=True)
            outdf = grouped_weighted_avg(df, df['weight'], ['ROW', 'COL'])
            return outdf

        return satellite.to_level3(
            self, *varkeys, grid=grid, griddims=griddims, weighting=weighting,
            as_dataset=as_dataset, verbose=verbose
        )
