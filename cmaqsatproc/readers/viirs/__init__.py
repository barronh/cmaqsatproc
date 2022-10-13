__all__ = ['L2_VIIRS_SNPP', 'AERDB_L2_VIIRS_SNPP', 'AERDT_L2_VIIRS_SNPP']
# https://ladsweb.modaps.eosdis.nasa.gov/opendap/allData/5110/AERDB_L2_VIIRS_SNPP/contents.html
from .. import satellite


class L2_VIIRS_SNPP(satellite):
    __doc__ = """
    VIIRS SNPP
    * valid = Land_Ocean_Quality_Flag > isvalid (default 2)
    * pixel corners are based on interpolated lat/lon
    """
    _defaultkeys = ('Optical_Depth_Land_And_Ocean', 'Land_Ocean_Quality_Flag')

    @classmethod
    def _open_hierarchical_dataset(cls, path, bbox=None, isvalid=2, **kwargs):
        """
        Convenience function to promote groups geophysical_data and
        geolocation_data from groups into the main xarray.Dataset object.

        Arguments
        ---------
        path : str
            Path to a L2_VIIRS_SNPP OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            Minimum value of flag for valid date (flag>=isvalid)
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: L2_VIIRS_SNPP
            Satellite processing instance
        """
        import xarray as xr
        datakey = 'geophysical_data'
        geokey = 'geolocation_data'

        dss = [
            xr.open_dataset(path, group=datakey, **kwargs),
            xr.open_dataset(path, group=geokey, **kwargs),
        ]
        ds = xr.merge(dss)
        ds = cls.prep_dataset(ds, bbox=bbox, isvalid=isvalid, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def open_dataset(cls, path, bbox=None, isvalid=2, **kwargs):
        """
        Open a local or remote path as a VIIRS satellite processor.

        Arguments
        ---------
        path : str
            Path to a L2_VIIRS_SNPP OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            Minimum value of flag for valid date (flag>=isvalid)
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: L2_VIIRS_SNPP
            Satellite processing instance
        """
        import xarray as xr

        ds = xr.open_dataset(path, **kwargs).reset_coords()
        if len(ds.dims) == 0:
            return cls._open_hierarchical_dataset(
                path, bbox=bbox, isvalid=isvalid, **kwargs
            )
        ds = cls.prep_dataset(ds, bbox=bbox, path=path, isvalid=isvalid)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def cmr_links(cls, method='opendap', **kwds):
        """
        Reimplemenation of satellite.cmr_links to account for html links in
        VIIRS OpenDAP links in the CMR.

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download or OpenDAP or s3
        """
        from copy import copy
        from ...utils import getcmrlinks
        kwds = copy(kwds)
        down_f = (
            lambda x: (
                'opendap' not in x['href']
                and (
                    x['href'].endswith('he5')
                    or x['href'].endswith('nc')
                )
                and x['href'].startswith('http')
            )
        )
        s3_f = (
            lambda x: (
                'opendap' not in x['href']
                and (
                    x['href'].endswith('he5')
                    or x['href'].endswith('nc')
                )
                and x['href'].startswith('s3')
            )
        )
        open_f = (
            lambda x: (
                'opendap' in x['href'] and x['href'].endswith('.nc.html')
            )
        )
        kwds.setdefault(
            'filterfunc', {
                'opendap': open_f, 's3': s3_f, 'download': down_f
            }[method]
        )

        links = getcmrlinks(**kwds)
        if method == 'opendap':
            links = [link.replace('.html', '') for link in links]
        return links


class AERDT_L2_VIIRS_SNPP(L2_VIIRS_SNPP):
    __doc__ = """
    AERDT_L2_VIIRS_SNPP
    * valid = Land_Ocean_Quality_Flag > isvalid (default 2)
    * pixel corners are based on interpolated lat/lon
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around L2_VIIRS_SNPP.cmr_links with short_name set to
        AERDT_L2_VIIRS_SNPP.

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download or OpenDAP or s3
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'AERDT_L2_VIIRS_SNPP')
        return L2_VIIRS_SNPP.cmr_links(method=method, **kwargs)

    @classmethod
    def prep_dataset(cls, ds, bbox=None, isvalid=2, path=None):
        """
        Applies spatial subset based on Latitude and Longitude, also
        interpolates pixel centers to corners and applies valid flags:
            * Land_Ocean_Quality_Flag >= isvalid
            * ~Optical_Depth_Land_And_Ocean.isnull()

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            Value greater than or equal to 2 are valid.
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
            df = ds[['latitude', 'longitude']].to_dataframe()
            df = df.query(
                f'latitude >= {swlat} and latitude <= {nelat}'
                + f' and longitude >= {swlon} and longitude <= {nelon}'
            )
            lines = df.index.get_level_values('number_of_lines_8x8').unique()
            if len(lines) < 0:
                raise ValueError(f'{path} has no values in {bbox}')
            ds = ds.isel(
                number_of_lines_8x8=slice(lines.min(), lines.max() + 1)
            )
        scanline = ds.number_of_lines_8x8.values
        scanline_edges = xr.DataArray(
            np.concatenate([scanline[:1], scanline[1:] - 0.5, scanline[-1:]]),
            dims=('number_of_lines_8x8',)
        )
        pixel = ds.number_of_pixels_8x8.values
        pixel_edges = xr.DataArray(
            np.concatenate([pixel[:1], pixel[1:] - 0.5, pixel[-1:]]),
            dims=('number_of_pixels_8x8',)
        )

        lat_edges = ds.latitude.interp(
            number_of_lines_8x8=scanline_edges,
            number_of_pixels_8x8=pixel_edges,
        )
        lon_edges = ds.longitude.interp(
            number_of_lines_8x8=scanline_edges,
            number_of_pixels_8x8=pixel_edges,
        )
        dims = ('number_of_lines_8x8', 'number_of_pixels_8x8')
        corner_slices = {
            'll': (slice(None, -1), slice(None, -1)),
            'lu': (slice(None, -1), slice(1, None)),
            'uu': (slice(1, None), slice(1, None)),
            'ul': (slice(1, None), slice(None, -1)),
        }
        coords = {
            'number_of_lines_8x8': ds.coords['number_of_lines_8x8'],
            'number_of_pixels_8x8': ds.coords['number_of_pixels_8x8'],
        }
        for cornerkey, corner_slice in corner_slices.items():
            ds[f'{cornerkey}_y'] = xr.DataArray(
                lat_edges[corner_slice], dims=dims, coords=coords
            )
            ds[f'{cornerkey}_x'] = xr.DataArray(
                lon_edges[corner_slice], dims=dims, coords=coords
            )
        ds['valid'] = (
            (ds['Land_Ocean_Quality_Flag'] >= isvalid)
            & (~ds['Optical_Depth_Land_And_Ocean'].isnull())
        )
        return ds

    @classmethod
    def cmaq_process(cls, qf, l3):
        """
        Process CMAQ as though it were observed by VIIRS, which is simply based
        on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., NH3)
        l3 : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        saodkey = 'Optical_Depth_Land_And_Ocean'
        overf = qf.csp.mean_overpass(satellite='aura').where(
            ~l3[saodkey].isnull()
        )
        return overf


class AERDB_L2_VIIRS_SNPP(L2_VIIRS_SNPP):
    __doc__ = """
    AERDB_L2_VIIRS_SNPP
    * valid = Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate not None
    * pixel corners are based on interpolated lat/lon
    """

    _defaultkeys = (
        'Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate',
        'Aerosol_Optical_Thickness_550_Ocean_Best_Estimate',
        'Aerosol_Optical_Thickness_550_Land_Best_Estimate',
        'Aerosol_Optical_Thickness_QA_Flag_Ocean',
        'Aerosol_Optical_Thickness_QA_Flag_Land'
    )

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around L2_VIIRS_SNPP.cmr_links with short_name set to
        AERDT_L2_VIIRS_SNPP.

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download or OpenDAP or s3
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'AERDB_L2_VIIRS_SNPP')
        return L2_VIIRS_SNPP.cmr_links(method=method, **kwargs)

    @classmethod
    def prep_dataset(cls, ds, bbox=None, isvalid=2, path=None):
        """
        Applies spatial subset based on Latitude and Longitude, also
        interpolates pixel centers to corners and applies valid flags:
            * AOT_QA_Flag_Land >= isvalid, or
            * AOT_QA_Flag_Ocean >= isvalid
            * and ~AOT_550_Land_Ocean_Best_Estimate.isnull()
            * where AOT = Aerosol_Optical_Thickness

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            Value greater than or equal to 2 are valid.
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
            df = ds[['Latitude', 'Longitude']].to_dataframe()
            df = df.query(
                f'Latitude >= {swlat} and Latitude <= {nelat}'
                + f' and Longitude >= {swlon} and Longitude <= {nelon}'
            )
            lines = df.index.get_level_values('Idx_Atrack').unique()
            if len(lines) < 0:
                raise ValueError(f'{path} has no values in {bbox}')
            ds = ds.sel(
                Idx_Atrack=slice(lines.min(), lines.max() + 1)
            )
        scanline = ds.Idx_Atrack.values
        scanline_edges = xr.DataArray(
            np.concatenate([scanline[:1], scanline[1:] - 0.5, scanline[-1:]]),
            dims=('Idx_Atrack',)
        )
        pixel = ds.Idx_Xtrack.values
        pixel_edges = xr.DataArray(
            np.concatenate([pixel[:1], pixel[1:] - 0.5, pixel[-1:]]),
            dims=('Idx_Xtrack',)
        )

        lat_edges = ds.Latitude.interp(
            Idx_Atrack=scanline_edges,
            Idx_Xtrack=pixel_edges,
        )
        lon_edges = ds.Longitude.interp(
            Idx_Atrack=scanline_edges,
            Idx_Xtrack=pixel_edges,
        )
        dims = ('Idx_Atrack', 'Idx_Xtrack')
        corner_slices = {
            'll': (slice(None, -1), slice(None, -1)),
            'lu': (slice(None, -1), slice(1, None)),
            'uu': (slice(1, None), slice(1, None)),
            'ul': (slice(1, None), slice(None, -1)),
        }
        coords = {
            'Idx_Atrack': ds.coords['Idx_Atrack'],
            'Idx_Xtrack': ds.coords['Idx_Xtrack'],
        }
        ds['cn_x'] = ds['Longitude']
        ds['cn_y'] = ds['Latitude']
        for cornerkey, corner_slice in corner_slices.items():
            ds[f'{cornerkey}_y'] = xr.DataArray(
                lat_edges[corner_slice], dims=dims, coords=coords
            )
            ds[f'{cornerkey}_x'] = xr.DataArray(
                lon_edges[corner_slice], dims=dims, coords=coords
            )
        ds['valid'] = (
            (
                (ds['Aerosol_Optical_Thickness_QA_Flag_Land'] >= isvalid)
                | (ds['Aerosol_Optical_Thickness_QA_Flag_Ocean'] >= isvalid)
            )
            & (~ds[
                'Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'
            ].isnull())
        )
        return ds

    @classmethod
    def cmaq_process(cls, qf, satl3f):
        """
        Process CMAQ as though it were observed by VIIRS, which is simply based
        on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., NH3)
        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        saodkey = 'Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'
        overf = qf.csp.mean_overpass(satellite='aura').where(
            ~satl3f[saodkey].isnull()
        )
        return overf
