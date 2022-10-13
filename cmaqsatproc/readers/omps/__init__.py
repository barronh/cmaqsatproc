__all__ = ['OMPS_NPP_L2', 'OMPS_NPP_NMNO2_L2', 'OMPS_NPP_NMTO3_L2']
from ..core import satellite


class OMPS_NPP_L2(satellite):
    __doc__ = """
    Default OMPS satellite processor.
    """
    _defaultkeys = ()

    @classmethod
    def _open_hierarchical_dataset(cls, path, bbox=None, **kwargs):
        """
        Convenience function to open h5 datasets with SienceData and
        GeolocationData variables promoted to the main dataset.

        Arguments
        ---------
        path : str
            Path to a OMPS L2 OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: OMPS_NPP_L2
            Satellite processing instance
        """
        import xarray as xr
        datakey = 'ScienceData'
        geokey = 'GeolocationData'
        # inpkey = 'InputPointers'

        dss = [
            xr.open_dataset(path, group=datakey, **kwargs),
            xr.open_dataset(path, group=geokey, **kwargs),
            # xr.open_dataset(path, group=inpkey, **kwargs),
        ]
        ds = xr.merge(dss)
        ds = cls.prep_dataset(ds, bbox=bbox, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a OMPS L2 OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: OMPS_NPP_L2
            Satellite processing instance
        """
        import xarray as xr

        if path.startswith('http'):
            print(
                'Currently, h5netcdf with h5pyd fails to read from NASA due to'
                + ' redirection. netCDF4 fails because all dimensions are'
                + ' UNLIMITED.'
            )
        kwargs.setdefault('engine', 'h5netcdf')
        ds = xr.open_dataset(path, **kwargs)
        if len(ds.data_vars) == 0:
            return cls._open_hierarchical_dataset(
                path, bbox=bbox, **kwargs
            )
        rename = {
            k: k.replace('ScienceData_', '').replace(
                'GeolocationData_', ''
            ).replace(
                'InputPointers_', ''
            )
            for k in ds.variables
        }
        rename.update({k: k.replace('ScienceData_', '') for k in ds.dims})
        ds = ds.rename(rename)
        ds = cls.prep_dataset(ds, bbox=bbox, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        sat.bbox = bbox
        return sat

    @classmethod
    def prep_dataset(cls, ds, bbox=None, path=None):
        """
        Defines valid based on Latitude/Longitude and bbox, and inteprolates
        pixel centers to corners.

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        path : str
            Unused.

        Returns
        -------
        ds : xarray.Dataset
        """
        import xarray as xr
        import numpy as np
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            lldf = ds[['Latitude', 'Longitude']].to_dataframe().query(
                f'Latitude >= {swlat} and Latitude <= {nelat}'
                + f' and Longitude >= {swlon} and Longitude <= {nelon}'
            )
            sval = lldf.index.get_level_values('DimAlongTrack').unique()
            # Not subsetting pixel dimension, because I want to use
            # the cross ways dimension in the interpolation to corners
            if len(sval) < 0:
                raise ValueError(f'{path} has no values in {bbox}')

            ds = ds.sel(
                DimAlongTrack=slice(max(sval.min() - 1, 0), sval.max() + 1)
            )

        DimAlongTrack = ds.DimAlongTrack.values
        DimAlongTrack_edges = xr.DataArray(
            np.concatenate([
                DimAlongTrack[:1], DimAlongTrack[1:] - 0.5, DimAlongTrack[-1:]
            ]),
            dims=('DimAlongTrack',)
        )
        pixel = ds.DimCrossTrack.values
        pixel_edges = xr.DataArray(
            np.concatenate([pixel[:1], pixel[1:] - 0.5, pixel[-1:]]),
            dims=('DimCrossTrack',)
        )
        if (
            'LatitudeCorner' in ds.variables
            and 'LongitudeCorner' in ds.variables
        ):
            lat_bnds = ds.LatitudeCorner
            lon_bnds = ds.LongitudeCorner
            dims = ('DimAlongTrack', 'DimCrossTrack')
            corner_slices = {
                'll': 0, 'ul': 1, 'uu': 2, 'lu': 3,
            }
            for cornerkey, corner_slice in corner_slices.items():
                ds[f'{cornerkey}_y'] = lat_bnds.isel(DimCorners=corner_slice)
                ds[f'{cornerkey}_x'] = lon_bnds.isel(DimCorners=corner_slice)
        else:
            lat_edges = ds.Latitude.interp(
                DimAlongTrack=DimAlongTrack_edges, DimCrossTrack=pixel_edges,
            )
            lon_edges = ds.Longitude.interp(
                DimAlongTrack=DimAlongTrack_edges, DimCrossTrack=pixel_edges,
            )
            dims = ('DimAlongTrack', 'DimCrossTrack')
            corner_slices = {
                'll': (slice(None, -1), slice(None, -1)),
                'lu': (slice(None, -1), slice(1, None)),
                'uu': (slice(1, None), slice(1, None)),
                'ul': (slice(1, None), slice(None, -1)),
            }
            coords = {
                'DimAlongTrack': ds.coords['DimAlongTrack'],
                'DimCrossTrack': ds.coords['DimCrossTrack'],
            }
            for cornerkey, corner_slice in corner_slices.items():
                ds[f'{cornerkey}_y'] = xr.DataArray(
                    lat_edges[corner_slice], dims=dims, coords=coords
                )
                ds[f'{cornerkey}_x'] = xr.DataArray(
                    lon_edges[corner_slice], dims=dims, coords=coords
                )
        ds['cn_x'] = ds['Longitude']
        ds['cn_y'] = ds['Latitude']

        if bbox is not None:
            ds['valid'] = (
                (ds['Latitude'] >= swlat) & (ds['Latitude'] <= nelat)
                & (ds['Longitude'] >= swlon) & (ds['Longitude'] <= nelon)
            )
        else:
            ds['valid'] = ds['Latitude'] > -999

        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        return ds


class OMPS_NPP_NMNO2_L2(OMPS_NPP_L2):
    _defaultkeys = (
        'ColumnAmountNO2', 'ColumnAmountNO2tropo', 'ColumnAmountNO2strat',
        'CloudFraction'
    )
    __doc__ = """
    Default OMPS NO2 satellite processor.
    * valid when
      * GroundPixelQualityFlags == 0
      * PixelQualityFlags == 0
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "OMPS_NPP_NMNO2_L2".

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
        kwargs.setdefault('short_name', 'OMPS_NPP_NMNO2_L2')
        return OMPS_NPP_L2.cmr_links(method=method, **kwargs)

    @classmethod
    def prep_dataset(cls, ds, bbox=None, path=None):
        """
        Defines valid based on Latitude/Longitude and bbox, and inteprolates
        pixel centers to corners. Then further subsets valid pixels where
        GroundPixel and Pixel QualityFlags are 0.

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        path : str
            Unused.

        Returns
        -------
        ds : xarray.Dataset
        """
        OMPS_NPP_L2.prep_dataset(ds, bbox=bbox, path=path)
        ds['valid'] = (
            ds['valid']
            & (ds['GroundPixelQualityFlags'] == 0)
            & (ds['PixelQualityFlags'] == 0)
        )
        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        return ds


class OMPS_NPP_NMTO3_L2(OMPS_NPP_L2):
    _defaultkeys = ('ColumnAmountO3')
    __doc__ = """
    Default OMPS O3 satellite processor.
    * valid when
      * GroundPixelQualityFlags == 0
      * QualityFlags == 0
    """

    @classmethod
    def prep_dataset(cls, ds, bbox=None, path=None):
        """
        Defines valid based on Latitude/Longitude and bbox, and inteprolates
        pixel centers to corners. Then further subsets valid pixels where
        GroundPixelQualityFlags and QualityFlags are 0.

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        path : str
            Unused.

        Returns
        -------
        ds : xarray.Dataset
        """
        ds = OMPS_NPP_L2.prep_dataset(ds, bbox=bbox, path=path)
        ds['valid'] = (
            ds['valid']
            & (ds['GroundPixelQualityFlags'] == 0)
            & (ds['QualityFlags'] == 0)
        )
        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        return ds

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "OMPS_NPP_NMTO3_L2".

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
        kwargs.setdefault('short_name', 'OMPS_NPP_NMTO3_L2')
        return OMPS_NPP_L2.cmr_links(method=method, **kwargs)
