__all__ = [
    'TropOMI', 'S5P_L2__NO2___', 'S5P_L2__CO____', 'S5P_L2__HCHO__',
    'S5P_L2__CH4___'
]
from ..core import satellite
import numpy as np


# Hard-coded to match current TropOMI implementation. This allows users
# to process less data.
tm5_constant_a = np.array([
    [0.0000000e+00, 2.2835938e+01, 4.2441406e+02, 1.3875469e+03,
     3.0572656e+03, 5.5643828e+03, 8.1633750e+03, 1.1901340e+04,
     1.4898453e+04, 1.7471840e+04, 1.9290227e+04, 2.0361816e+04,
     2.0337863e+04, 1.9859391e+04, 1.9031289e+04, 1.8308434e+04,
     1.7008789e+04, 1.5508257e+04, 1.3881331e+04, 1.2766873e+04,
     1.1116662e+04, 9.5626826e+03, 8.6085254e+03, 7.3118691e+03,
     6.1560742e+03, 4.4908174e+03, 3.3817437e+03, 2.2654316e+03,
     1.4237701e+03, 8.2396783e+02, 4.2759250e+02, 1.9133856e+02,
     6.9520576e+01, 1.8608931e+01],
    [2.2835938e+01, 4.2441406e+02, 1.3875469e+03, 3.0572656e+03,
     5.5643828e+03, 8.1633750e+03, 1.1901340e+04, 1.4898453e+04,
     1.7471840e+04, 1.9290227e+04, 2.0361816e+04, 2.0337863e+04,
     1.9859391e+04, 1.9031289e+04, 1.8308434e+04, 1.7008789e+04,
     1.5508257e+04, 1.3881331e+04, 1.2766873e+04, 1.1116662e+04,
     9.5626826e+03, 8.6085254e+03, 7.3118691e+03, 6.1560742e+03,
     4.4908174e+03, 3.3817437e+03, 2.2654316e+03, 1.4237701e+03,
     8.2396783e+02, 4.2759250e+02, 1.9133856e+02, 6.9520576e+01,
     1.8608931e+01, 0.0000000e+00]
], dtype='f').T

tm5_constant_b = np.array([
    [1.00000e+00, 9.91984e-01, 9.69513e-01, 9.31881e-01, 8.73929e-01,
     7.90717e-01, 7.04669e-01, 5.76692e-01, 4.66003e-01, 3.58254e-01,
     2.63242e-01, 1.68910e-01, 1.11505e-01, 7.79580e-02, 5.17730e-02,
     3.80260e-02, 2.23550e-02, 1.18060e-02, 5.37800e-03, 2.85700e-03,
     8.90000e-04, 1.99000e-04, 5.90000e-05, 0.00000e+00, 0.00000e+00,
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00],
    [9.91984e-01, 9.69513e-01, 9.31881e-01, 8.73929e-01, 7.90717e-01,
     7.04669e-01, 5.76692e-01, 4.66003e-01, 3.58254e-01, 2.63242e-01,
     1.68910e-01, 1.11505e-01, 7.79580e-02, 5.17730e-02, 3.80260e-02,
     2.23550e-02, 1.18060e-02, 5.37800e-03, 2.85700e-03, 8.90000e-04,
     1.99000e-04, 5.90000e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00]
], dtype='f').T


class TropOMI(satellite):
    __doc__ = """
    Default TropOMI satellite processor.
    * bbox subsets the scanline dimensions
    """
    _defaultkeys = ()

    @classmethod
    def _open_hierarchical_dataset(
        cls, path, bbox=None, isvalid=0.75, **kwargs
    ):
        """
        Convenience function to promote PRODUCT GEOLOCATIONS, DETAILED_RESULTS,
        and INPUT_DATA groups' variables to the main xarray.Dataset object.

        Arguments
        ---------
        path : str
            Path to a TropOMI he5-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            Pixels are valid with the qa_value is greater than isvalid
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: TropOMI
            Satellite processing instance
        """
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
    def open_dataset(cls, path, bbox=None, isvalid=0.75, **kwargs):
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
        kwargs.setdefault('decode_coords', False)
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
    def prep_dataset(cls, ds, bbox=None, isvalid=0.75, path=None):
        """
        Defines pixels as valid when they are within bbox and interpolates
        latitude and longitude from pixel centers to corners.

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        isvalid : float
            When qa_value is greater than isvalid, the pixel is valid.
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
        ds['cn_x'] = ds['longitude']
        ds['cn_y'] = ds['latitude']

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
                ds[f'{cornerkey}_y'] = lat_bnds.isel(corner=corner_slice)
                ds[f'{cornerkey}_x'] = lon_bnds.isel(corner=corner_slice)
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
        if bbox is not None:
            ds['valid'] = ds['valid'] & (
                (ds['latitude'] >= swlat) & (ds['latitude'] <= nelat)
                & (ds['longitude'] >= swlon) & (ds['longitude'] <= nelon)
            )

        # The HCHO OpenDAP file constants have dimensions ('time', 'layer')
        for vkey in ['tm5_constant_a', 'tm5_constant_b']:
            if vkey in ds.data_vars:
                vvar = ds[vkey]
                if 'time' in vvar.dims:
                    ds[vkey] = vvar.isel(time=0)

        if not ds['valid'].any():
            import warnings
            warnings.warn('No valid pixels')

        return ds

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

    @classmethod
    def cmaq_sw(cls, overf, satl3f, amfkey, tropopausekey=None):
        """
        Calculate scattering weights as averaging kernel multiplied by the
        total air mass factor. Then, interpolate to CMAQ vertical grid, based
        on PRES (pressure in Pa). The scattering weights are set to 0 above
        the tropopause.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key for the total air mass factor.
        tropopauskey : str
            Key for the tropopause when available. Otherwise, use None.

        Returns
        -------
        q_sw_trop : xr.DataArray
            Tropospheric scattering Weights on the CMAQ grid
        """
        from ...utils import coord_interp
        import xarray as xr

        if 'tm5_constant_b' in satl3f:
            tm5_b = satl3f['tm5_constant_b']
            if 'vertices' in tm5_b.dims:
                tm5_b = tm5_b.mean('vertices')
        else:
            tm5_b = xr.DataArray(
                tm5_constant_b.mean(1), dims=('layer',),
                coords=dict(layer=satl3f['layer'])
            )
        if 'tm5_constant_a' in satl3f:
            tm5_a = satl3f['tm5_constant_a']
            if 'vertices' in tm5_a.dims:
                tm5_a = tm5_a.mean('vertices')
        else:
            tm5_a = xr.DataArray(
                tm5_constant_a.mean(1), dims=('layer',),
                coords=dict(layer=satl3f['layer'])
            )

        pres = tm5_b * satl3f['surface_pressure'] + tm5_a
        ak = satl3f['averaging_kernel']
        if tropopausekey is not None:
            ak = ak.where(satl3f[tropopausekey] >= satl3f['layer'])

        sw_trop = ak * satl3f[amfkey]

        q_sw_trop = coord_interp(
            overf['PRES'], pres, sw_trop,
            indim='layer', outdim='LAY', ascending=False
        )
        q_sw_trop.attrs.update(satl3f['averaging_kernel'].attrs)
        q_sw_trop.attrs['var_desc'] = 'AK * AMF'.ljust(80)
        return q_sw_trop

    @classmethod
    def cmaq_amf(cls, overf, satl3f, amfkey, key):
        """
        Calculate an alternative Air Mass Factor (AMF) using satellite
        scattering weights and the CMAQ vertical profile as a partial column
        density.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key of the total air mass factor
        key : str
            Key of the partial column density variable from CMAQ, which must
            have a LAY dimension that describes teh vertical coordinate.

        Returns
        -------
        cmaqamf : xr.DataArray
            Air Mass Factor on the CMAQ grid
        """
        q_sw_trop = cls.cmaq_sw(overf, satl3f, amfkey=amfkey)
        denom = overf[key].sum('LAY')
        q_amf = (q_sw_trop * overf[key]).sum('LAY') / denom
        q_amf = q_amf.where((denom != 0) & ~q_sw_trop.isnull().all('LAY'))
        q_amf.attrs.update(q_sw_trop.attrs)
        return q_amf


class S5P_L2__CO____(TropOMI):
    _defaultkeys = (
        'carbonmonoxide_total_column', 'column_averaging_kernel', 'layer'
    )
    __doc__ = """
    TropOMICO satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.75)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "S5P_L2__CO____" or "S5P_L2__CO_____HiR".

        The HiR product started 2019-08-06T02:41:41.000Z

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
        hirstartdate = '2019-08-06'
        querydate = kwargs.get('temporal', hirstartdate)[:10]
        if querydate >= hirstartdate:
            defshortname = 'S5P_L2__CO_____HiR'
        else:
            defshortname = 'S5P_L2__CO____'
        kwargs.setdefault('short_name', defshortname)
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_process(cls, qf, l3, key='CO'):
        """
        Process CMAQ as though it were observed by TropOMI, which is simply
        based on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., CO)
        l3 : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        import xarray as xr

        skey = 'carbonmonoxide_total_column'
        overf = qf.csp.mean_overpass(satellite='aura').where(
            ~l3[skey].isnull()
        )
        n_per_m2 = overf.csp.mole_per_m2(add=True)
        tgtvar = overf[key]
        if tgtvar.units.strip().startswith('ppm'):
            vmr = tgtvar / 1e6
        elif tgtvar.units.strip().startswith('ppb'):
            vmr = tgtvar / 1e9
        elif tgtvar.units.strip().startswith('ppt'):
            vmr = tgtvar / 1e12
        overf['CO_MOL_PER_M2'] = n_per_m2 * vmr
        overf['CO_MOL_PER_M2'].attrs.update(overf[key].attrs)
        overf['CO_MOL_PER_M2'].attrs['units'] = 'mole/m**2'.ljust(16)
        overf['VCDCO_CMAQ'] = overf['CO_MOL_PER_M2'].sum('LAY')

        # For CO, the model should be integrated to meter-based vertical levels
        # https://sentinels.copernicus.eu/documents/247904/3541451/Sentinel-5P-
        # Carbon-Monoxide-Level-2-Product-Readme-File.pdf/f8942626-ffb6-4951-90
        # fc-a16b6589e39e?t=1658386616101
        dz = np.diff(l3['layer']).mean() / 2.
        tops = np.append(0, l3['layer'] + dz)
        akvar = l3['column_averaging_kernel']
        if akvar.attrs['units'].strip() == 'm':
            # Older than v2.0.4 uses needs to be divided by 1000m
            # which is the depth of each layer.
            ak = akvar[:] / 1000
        else:
            ak = akvar[:]

        nl = l3.dims['layer']
        CO_MOL_PER_M2_Z = xr.DataArray(
            np.zeros((nl, overf.dims['ROW'], overf.dims['COL']), dtype='f'),
            dims=('layer', 'ROW', 'COL'), name='CO_MOL_PER_M2_Z'
        )
        overf['VCDCO_CMAQ_TOMI'] = xr.DataArray(
            np.zeros((overf.dims['ROW'], overf.dims['COL']), dtype='f'),
            dims=('ROW', 'COL'), name='VCDCO_CMAQ_TOMI'
        )
        skip = ak.isnull().all('layer')
        # Currently looping over rows, which is inefficient.
        for ri, r in enumerate(overf['ROW']):
            # print('.', end='', flush=True)
            if skip[ri].all():
                continue
            for ci, c in enumerate(overf['COL']):
                if skip[ri, ci]:
                    continue
                cdfp = np.append(
                    0, np.cumsum(overf['CO_MOL_PER_M2'][:, ri, ci])
                )
                zfp = np.append(0, overf['ZF'][:, ri, ci])
                cdf = np.interp(tops, zfp, cdfp)
                pmf = np.diff(cdf)
                CO_MOL_PER_M2_Z[:, ri, ci] = pmf
                overf['VCDCO_CMAQ_TOMI'][ri, ci] = (
                    CO_MOL_PER_M2_Z[:, ri, ci] * ak[ri, ci]
                ).sum()

        # overf['CO_MOL_PER_M2_Z'] = CO_MOL_PER_M2_Z.where(~skip)
        overf['VCDCO_CMAQ_TOMI'] = overf['VCDCO_CMAQ_TOMI'].where(~skip)
        return overf


class S5P_L2__NO2___(TropOMI):
    _defaultkeys = (
        'air_mass_factor_total', 'air_mass_factor_troposphere',
        'nitrogendioxide_tropospheric_column', 'nitrogendioxide_total_column',
        'averaging_kernel', 'tm5_tropopause_layer_index', 'surface_pressure'
    )
    __doc__ = """
    TropOMINO2 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.75)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "S5P_L2__NO2___" or "S5P_L2__NO2____HiR".

        "S5P_L2__NO2____HiR" starts 2019-08-06T02:41:41.000Z

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
        hirstartdate = '2019-08-06'
        querydate = kwargs.get('temporal', hirstartdate)[:10]
        if querydate >= hirstartdate:
            defshortname = 'S5P_L2__NO2____HiR'
        else:
            defshortname = 'S5P_L2__NO2___'
        kwargs.setdefault('short_name', defshortname)
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_sw(
        cls, overf, satl3f, amfkey='air_mass_factor_total',
        tropopausekey='tm5_tropopause_layer_index'
    ):
        """
        Calculate scattering weights as averaging kernel multiplied by the
        total air mass factor. Then, interpolate to CMAQ vertical grid, based
        on PRES (pressure in Pa). The scattering weights are set to 0 above
        the tropopause.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key for the total air mass factor.
        tropopausekey : str
            Key for the tropopause index.

        Returns
        -------
        q_sw_trop : xr.DataArray
            Tropospheric scattering Weights on the CMAQ grid
        """
        return TropOMI.cmaq_sw(
            overf, satl3f, amfkey=amfkey, tropopausekey=tropopausekey
        )

    @classmethod
    def cmaq_amf(
        cls, overf, satl3f, amfkey='air_mass_factor_total', key='NO2_PER_M2'
    ):
        """
        Calculate an alternative Air Mass Factor (AMF) using satellite
        scattering weights and the CMAQ vertical profile as a partial column
        density.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key of the total air mass factor
        key : str
            Key of the partial column density variable from CMAQ, which must
            have a LAY dimension that describes teh vertical coordinate.

        Returns
        -------
        cmaqamf : xr.DataArray
            Air Mass Factor on the CMAQ grid
        """
        return TropOMI.cmaq_amf(overf, satl3f, amfkey=amfkey, key=key)

    @classmethod
    def cmaq_ak(cls, overf, satl3f, amfkey='air_mass_factor_total'):
        """
        Calculate an averaging kernel (AK) that would process CMAQ as though
        it were observed by the satellite. In this case, the averaging kernel
        is the scattering weights divided by the tropospheric air mass factor.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        q_ak : xr.DataArray
            Averaging kernel on the CMAQ grid
        """
        q_sw_trop = cls.cmaq_sw(overf, satl3f, amfkey=amfkey)
        q_ak = q_sw_trop / satl3f['air_mass_factor_troposphere']
        q_ak.attrs.update(q_sw_trop.attrs)
        return q_ak

    @classmethod
    def cmaq_process(cls, qf, satl3f, key='NO2'):
        """
        Process CMAQ as though it were observed by TropOMI and recalculate
        TropOMI tropospheric columns with the CMAQ AMF. This process relies
        on cmaq_ak and cmaq_amf.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., NO2), PRES, DENS, and ZF
            variables with a LAY dimension describing the vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        key : str
            Composition key

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ and CMAQ-like satellite.
        """
        # OMI is on the aura satellite, so we create an average overpass
        overf = qf.csp.mean_overpass(satellite='aura')
        n_per_m2 = overf.csp.mole_per_m2(add=True)
        tgtvar = overf[key]
        if tgtvar.units.strip().startswith('ppm'):
            vmr = tgtvar / 1e6
        elif tgtvar.units.strip().startswith('ppb'):
            vmr = tgtvar / 1e9
        elif tgtvar.units.strip().startswith('ppt'):
            vmr = tgtvar / 1e12
        overf['NO2_PER_M2'] = n_per_m2 * vmr
        overf['NO2_PER_M2'].attrs.update(overf[key].attrs)
        overf['NO2_PER_M2'].attrs['units'] = 'mole/m**2'.ljust(16)

        ak = overf['NO2_AK_CMAQ'] = cls.cmaq_ak(overf, satl3f)
        overf['NO2_SW_CMAQ'] = cls.cmaq_sw(overf, satl3f)
        # uses AK for tropopause
        overf['VCDNO2_CMAQ'] = overf.csp.apply_ak('NO2_PER_M2', ak / ak)
        # uses AK for vertical weighting
        overf['VCDNO2_CMAQ_TOMI'] = overf.csp.apply_ak('NO2_PER_M2', ak)
        # Recalculate satellite
        amf = overf['AMF_CMAQ'] = cls.cmaq_amf(overf, satl3f, key='NO2_PER_M2')

        overf['VCDNO2_TOMI_CMAQ'] = (
            satl3f['nitrogendioxide_tropospheric_column']
            * satl3f['air_mass_factor_troposphere'] / amf
        )
        overf['VCDNO2_TOMI_CMAQ'].attrs.update(
            satl3f['nitrogendioxide_tropospheric_column'].attrs
        )
        overf['VCDNO2_TOMI_CMAQ'].attrs['var_desc'] = (
            'TropOMI_NO2 x TropOMI_AMF / CMAQ_AMF'
        ).ljust(80)

        return overf


class S5P_L2__CH4___(TropOMI):
    _defaultkeys = ('methane_mixing_ratio', 'averaging_kernel')
    __doc__ = """
    TropOMICH4 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.75)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "S5P_L2__CH4___" or "S5P_L2__CH4____HiR".

        The HiR product started 2019-08-06T02:41:41.000Z

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
        hirstartdate = '2019-08-06'
        querydate = kwargs.get('temporal', hirstartdate)[:10]
        if querydate >= hirstartdate:
            defshortname = 'S5P_L2__CH4____HiR'
        else:
            defshortname = 'S5P_L2__CH4___'
        kwargs.setdefault('short_name', defshortname)
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_process(cls, qf, l3, key='CH4'):
        """
        Process CMAQ as though it were observed by TropOMI, which is simply
        based on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., ECH4)
        l3 : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        skey = 'methane_mixing_ratio'
        overf = qf.csp.mean_overpass(satellite='aura').where(
            ~l3[skey].isnull()
        )
        return overf


class S5P_L2__HCHO__(TropOMI):
    __doc__ = """
    TropOMIHCHO satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.75)
    """
    _defaultkeys = (
        'formaldehyde_profile_apriori', 'averaging_kernel',
        'formaldehyde_tropospheric_air_mass_factor',
        'formaldehyde_tropospheric_vertical_column',
        'surface_pressure'
    )

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where short_name is set to
        "S5P_L2__HCHO__" or "S5P_L2__HCHO___HiR".

        "S5P_L2_HCHO___HiR" starts 2019-08-06T02:41:41.000Z

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
        hirstartdate = '2019-08-06'
        querydate = kwargs.get('temporal', hirstartdate)[:10]
        if querydate >= hirstartdate:
            defshortname = 'S5P_L2__HCHO___HiR'
        else:
            defshortname = 'S5P_L2__HCHO__'
        kwargs.setdefault('short_name', defshortname)
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_sw(
        cls, overf, satl3f, amfkey='formaldehyde_tropospheric_air_mass_factor'
    ):
        """
        Calculate scattering weights as averaging kernel multiplied by the
        total air mass factor. Then, interpolate to CMAQ vertical grid, based
        on PRES (pressure in Pa). The scattering weights are set to 0 above
        the tropopause.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key for the total air mass factor.

        Returns
        -------
        q_sw_trop : xr.DataArray
            Tropospheric scattering Weights on the CMAQ grid
        """
        return TropOMI.cmaq_sw(overf, satl3f, amfkey=amfkey)

    @classmethod
    def cmaq_amf(
        cls, overf, satl3f,
        amfkey='formaldehyde_tropospheric_air_mass_factor', key='FORM_PER_M2'
    ):
        """
        Calculate an alternative Air Mass Factor (AMF) using satellite
        scattering weights and the CMAQ vertical profile as a partial column
        density.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        amfkey : str
            Key of the total air mass factor
        key : str
            Key of the partial column density variable from CMAQ, which must
            have a LAY dimension that describes teh vertical coordinate.

        Returns
        -------
        cmaqamf : xr.DataArray
            Air Mass Factor on the CMAQ grid
        """
        return TropOMI.cmaq_amf(overf, satl3f, amfkey=amfkey, key=key)

    @classmethod
    def cmaq_ak(cls, overf, satl3f):
        """
        Calculate an averaging kernel (AK) that would process CMAQ as though
        it were observed by the satellite. In this case, the averaging kernel
        is the scattering weights divided by the tropospheric air mass factor.

        Arguments
        ---------
        overf : xarray.Dataset
            Must have PRES variable with LAY dimension that describes the
            vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).

        Returns
        -------
        q_ak : xr.DataArray
            Averaging kernel on the CMAQ grid
        """
        q_sw_trop = cls.cmaq_sw(overf, satl3f)
        q_ak = q_sw_trop / satl3f['formaldehyde_tropospheric_air_mass_factor']
        return q_ak

    @classmethod
    def cmaq_process(cls, qf, satl3f, key='FORM'):
        """
        Process CMAQ as though it were observed by OMI and recalculate TropOMI
        tropospheric columns with the CMAQ AMF. This process relies on cmaq_ak
        and cmaq_amf.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., FORM), PRES, DENS, and ZF
            variables with a LAY dimension describing the vertical coordinate.

        satl3f : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        key : str
            Composition key

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ and CMAQ-like satellite.
        """
        # OMI is on the aura satellite, so we create an average overpass
        overf = qf.csp.mean_overpass(satellite='aura')
        n_per_m2 = overf.csp.mole_per_m2(add=True)
        tgtvar = overf[key]
        if tgtvar.units.strip().startswith('ppm'):
            vmr = tgtvar / 1e6
        elif tgtvar.units.strip().startswith('ppb'):
            vmr = tgtvar / 1e9
        elif tgtvar.units.strip().startswith('ppt'):
            vmr = tgtvar / 1e12
        overf['FORM_PER_M2'] = n_per_m2 * vmr

        ak = overf['FORM_AK_CMAQ'] = cls.cmaq_ak(overf, satl3f)
        overf['FORM_SW_CMAQ'] = cls.cmaq_sw(overf, satl3f)
        # uses AK for tropopause
        overf['VCDFORM_CMAQ'] = overf.csp.apply_ak('FORM_PER_M2', ak / ak)
        # uses AK for vertical weighting
        overf['VCDFORM_CMAQ_TOMI'] = overf.csp.apply_ak('FORM_PER_M2', ak)
        # Recalculate satellite
        amf = overf['AMF_CMAQ'] = cls.cmaq_amf(
            overf, satl3f, key='FORM_PER_M2'
        )

        overf['VCDHCHO_TOMI_CMAQ'] = (
            satl3f['formaldehyde_tropospheric_vertical_column']
            * satl3f['formaldehyde_tropospheric_air_mass_factor'] / amf
        )
        overf['VCDHCHO_TOMI_CMAQ'].attrs.update(
            satl3f['formaldehyde_tropospheric_vertical_column'].attrs
        )
        overf['VCDHCHO_TOMI_CMAQ'].attrs['var_desc'] = (
            'TropOMI_NO2 x TropOMI_AMF / CMAQ_AMF'
        ).ljust(80)

        return overf
