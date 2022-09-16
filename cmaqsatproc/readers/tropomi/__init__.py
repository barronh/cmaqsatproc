__all__ = ['TropOMI']
from ..core import satellite


class TropOMI(satellite):
    __doc__ = """
    Default TropOMI satellite processor.
    * bbox subsets the scanline dimensions
    """
    _defaultkeys = ()

    @classmethod
    def _open_hierarchical_dataset(cls, path, bbox=None, isvalid=0.5, **kwargs):
        import xarray as xr
        datakey = 'PRODUCT'
        geokey = 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS'
        detkey = 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS'
        inpkey = 'PRODUCT/SUPPORT_DATA/INPUT_DATA'

        dss = [
            xr.open_dataset(path, group=datakey, **kwargs),
            xr.open_dataset(path, group=geokey, **kwargs),
            xr.open_dataset(path, group=detkey, **kwargs),
            xr.open_dataset(path, group=inpkey, **kwargs),
        ]
        ds = xr.merge(dss)
        ds = cls.prep_dataset(ds, bbox=bbox, isvalid=isvalid, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def open_dataset(cls, path, bbox=None, isvalid=0.5, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a TropOMI OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        isvalid : float
            Value of qa_value to use for valid pixels
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: TropOMI
            Satellite processing instance
        """
        import xarray as xr
        ds = xr.open_dataset(path, **kwargs)
        if len(ds.dims) == 0:
            return cls._open_hierarchical_dataset(
                path, bbox=bbox, isvalid=isvalid, **kwargs
            )
        rename = {
            k: k.replace('PRODUCT_', '').replace(
                'SUPPORT_DATA_DETAILED_RESULTS_', ''
            ).replace(
                'SUPPORT_DATA_GEOLOCATIONS_', ''
            ).replace(
                'SUPPORT_DATA_INPUT_DATA_', ''
            ).replace('METADATA_QA_STATISTICS_', '')
            for k in ds.variables
        }
        rename.update({k: k.replace('PRODUCT_', '') for k in ds.dims})
        ds = ds.rename(rename)
        ds = cls.prep_dataset(ds, bbox=bbox, isvalid=isvalid, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        sat.bbox = bbox
        return sat

    @classmethod
    def prep_dataset(cls, ds, bbox=None, isvalid=0.5, path=None):
        import xarray as xr
        import numpy as np
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            lldf = ds[['latitude', 'longitude']].to_dataframe().query(
                f'latitude >= {swlat} and latitude <= {nelat}'
                + f' and longitude >= {swlon} and longitude <= {nelon}'
            )
            sval = lldf.index.get_level_values('scanline').unique()
            # Not subsetting pixel dimension, because I want to use
            # the cross ways dimension in the interpolation to corners
            if len(sval) < 0:
                raise ValueError(f'{path} has no values in {bbox}')

            ds = ds.sel(scanline=slice(sval.min() - 1, sval.max() + 1))

        scanline = ds.scanline.values
        scanline_edges = xr.DataArray(
            np.concatenate([scanline[:1], scanline[1:] - 0.5, scanline[-1:]]),
            dims=('scanline',)
        )
        pixel = ds.ground_pixel.values
        pixel_edges = xr.DataArray(
            np.concatenate([pixel[:1], pixel[1:] - 0.5, pixel[-1:]]),
            dims=('ground_pixel',)
        )
        if (
            'latitude_bounds' in ds.variables
            and 'longitude_bounds' in ds.variables
        ):
            lat_bnds = ds.latitude_bounds
            lon_bnds = ds.longitude_bounds
            dims = ('time', 'scanline', 'ground_pixel')
            corner_slices = {
                'll': 0, 'ul': 1, 'uu': 2, 'lu': 3,
            }
            for cornerkey, corner_slice in corner_slices.items():
                ds[f'{cornerkey}_y'] = lat_bnds.sel(corner=corner_slice)
                ds[f'{cornerkey}_x'] = lon_bnds.sel(corner=corner_slice)
        else:
            lat_edges = ds.latitude.interp(
                scanline=scanline_edges, ground_pixel=pixel_edges,
            )
            lon_edges = ds.longitude.interp(
                scanline=scanline_edges, ground_pixel=pixel_edges,
            )
            dims = ('time', 'scanline', 'ground_pixel')
            corner_slices = {
                'll': (slice(None,), slice(None, -1), slice(None, -1)),
                'lu': (slice(None,), slice(None, -1), slice(1, None)),
                'uu': (slice(None,), slice(1, None), slice(1, None)),
                'ul': (slice(None,), slice(1, None), slice(None, -1)),
            }
            coords = {
                'time': ds.coords['time'],
                'scanline': ds.coords['scanline'],
                'ground_pixel': ds.coords['ground_pixel'],
            }
            for cornerkey, corner_slice in corner_slices.items():
                ds[f'{cornerkey}_y'] = xr.DataArray(
                    lat_edges[corner_slice], dims=dims, coords=coords
                )
                ds[f'{cornerkey}_x'] = xr.DataArray(
                    lon_edges[corner_slice], dims=dims, coords=coords
                )
        ds['valid'] = ds['qa_value'] > isvalid
        return ds

    @classmethod
    def shorten_name(cls, key):
        key = key.replace('SUPPORT_DATA_DETAILED_RESULTS', 'SUP')
        key = key.replace('SUPPORT_DATA_INPUT_DATA', 'IN')
        key = key.replace('SUPPORT_DATA_GEOLOCATIONS', 'GEO')
        key = key.replace('total', 'TOT')
        key = key.replace('tropospheric', 'TROP')
        key = key.replace('stratospheric', 'STRAT')
        key = key.replace('column', 'VCD')
        key = key.replace('averaging_kernel', 'AK')
        key = key.replace('carbonmonoxide', 'CO')
        key = key.replace('nitrogendioxide', 'NO2')
        key = key.replace('surface', 'sfc')
        key = key.replace('pressure', 'pres')
        return key


class TropOMICO(TropOMI):
    _defaultkeys = ('carbonmonoxide_total_column',)
    __doc__ = """
    TropOMICO satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """


class TropOMINO2(TropOMI):
    _defaultkeys = (
        'air_mass_factor_troposphere', 'nitrogendioxide_tropospheric_column',
        'nitrogendioxide_total_column'
    )
    __doc__ = """
    TropOMINO2 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """


class TropOMICH4(TropOMI):
    _defaultkeys = ('methane_mixing_ratio', 'averaging_kernel')
    __doc__ = """
    TropOMICH4 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """
