__all__ = ['TropOMI', 'TropOMINO2', 'TropOMICO', 'TropOMIHCHO', 'TropOMICH4']
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
        cls, path, bbox=None, isvalid=0.5, **kwargs
    ):
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
    def cmaq_sw(cls, overf, outputs, amfkey, tropopausekey=None):
        from ...utils import coord_interp
        import xarray as xr

        if 'tm5_constant_b' in outputs:
            tm5_b = outputs['tm5_constant_b']
            if 'vertices' in tm5_b.dims:
                tm5_b = tm5_b.mean('vertices')
        else:
            tm5_b = xr.DataArray(
                tm5_constant_b.mean(1), dims=('layer',),
                coords=[('layer', outputs['layer'])]
            )
        if 'tm5_constant_a' in outputs:
            tm5_a = outputs['tm5_constant_a']
            if 'vertices' in tm5_a.dims:
                tm5_a = tm5_a.mean('vertices')
        else:
            tm5_a = xr.DataArray(
                tm5_constant_a.mean(1), dims=('layer',),
                coords=[('layer', outputs['layer'])]
            )

        pres = tm5_b * outputs['surface_pressure'] + tm5_a
        ak = outputs['averaging_kernel']
        if tropopausekey is not None:
            ak = ak.where(outputs[tropopausekey] >= outputs['layer'])

        sw_trop = ak * outputs[amfkey]

        q_sw_trop = coord_interp(
            overf['PRES'], pres, sw_trop,
            indim='layer', outdim='LAY', ascending=False
        )
        q_sw_trop.attrs.update(outputs['averaging_kernel'].attrs)
        q_sw_trop.attrs['var_desc'] = 'AK * AMF'.ljust(80)
        return q_sw_trop

    @classmethod
    def cmaq_amf(cls, overf, outputs, amfkey, key):
        """
        Calculates the Tropospheric Averaging Kernel
        """
        q_sw_trop = cls.cmaq_sw(overf, outputs, amfkey=amfkey)
        denom = overf[key].sum('LAY')
        q_amf = (q_sw_trop * overf[key]).sum('LAY') / denom
        q_amf = q_amf.where((denom != 0) & ~q_sw_trop.isnull().all('LAY'))
        q_amf.attrs.update(q_sw_trop.attrs)
        return q_amf


class TropOMICO(TropOMI):
    _defaultkeys = (
        'carbonmonoxide_total_column', 'column_averaging_kernel', 'layer'
    )
    __doc__ = """
    TropOMICO satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'S5P_L2__CO____')
        return TropOMI.cmr_links(method=method, **kwargs)


class TropOMINO2(TropOMI):
    _defaultkeys = (
        'air_mass_factor_total', 'air_mass_factor_troposphere',
        'nitrogendioxide_tropospheric_column', 'nitrogendioxide_total_column',
        'averaging_kernel', 'tm5_tropopause_layer_index', 'surface_pressure'
    )
    __doc__ = """
    TropOMINO2 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'S5P_L2__NO2___')
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_sw(
        cls, overf, outputs, amfkey='air_mass_factor_total',
        tropopausekey='tm5_tropopause_layer_index'
    ):
        return TropOMI.cmaq_sw(
            overf, outputs, amfkey=amfkey, tropopausekey=tropopausekey
        )

    @classmethod
    def cmaq_amf(
        cls, overf, outputs, amfkey='air_mass_factor_total', key='NO2_PER_M2'
    ):
        return TropOMI.cmaq_amf(overf, outputs, amfkey=amfkey, key=key)

    @classmethod
    def cmaq_ak(cls, overf, outputs, amfkey='air_mass_factor_total'):
        """
        Calculates the Tropospheric Averaging Kernel
        """
        q_sw_trop = cls.cmaq_sw(overf, outputs, amfkey=amfkey)
        q_ak = q_sw_trop / outputs['air_mass_factor_troposphere']
        q_ak.attrs.update(q_sw_trop.attrs)
        return q_ak

    @classmethod
    def cmaq_process(cls, qf, satl3f, key='NO2'):
        """
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
        overf['NO2_PER_M2'].attrs.update(overf['NO2'].attrs)
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


class TropOMICH4(TropOMI):
    _defaultkeys = ('methane_mixing_ratio', 'averaging_kernel')
    __doc__ = """
    TropOMICH4 satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'S5P_L2__CH4___')
        return TropOMI.cmr_links(method=method, **kwargs)


class TropOMIHCHO(TropOMI):
    __doc__ = """
    TropOMIHCHO satellite processor.
    * bbox subsets the scanline dimensions
    * valid = qa_value >= threshold (default 0.5)
    """
    _defaultkeys = (
        'formaldehyde_profile_apriori', 'averaging_kernel',
        'formaldehyde_tropospheric_air_mass_factor',
        'formaldehyde_tropospheric_vertical_column',
        'surface_pressure'
    )

    @classmethod
    def cmr_links(cls, method='opendap', **kwargs):
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('short_name', 'S5P_L2__HCHO__')
        return TropOMI.cmr_links(method=method, **kwargs)

    @classmethod
    def cmaq_sw(
        cls, overf, outputs, amfkey='formaldehyde_tropospheric_air_mass_factor'
    ):
        return TropOMI.cmaq_sw(overf, outputs, amfkey=amfkey)

    @classmethod
    def cmaq_amf(
        cls, overf, outputs,
        amfkey='formaldehyde_tropospheric_air_mass_factor', key='FORM_PER_M2'
    ):
        return TropOMI.cmaq_amf(overf, outputs, amfkey=amfkey, key=key)

    @classmethod
    def cmaq_ak(cls, overf, outputs):
        """
        Calculates the Tropospheric Averaging Kernel
        """
        q_sw_trop = cls.cmaq_sw(overf, outputs)
        q_ak = q_sw_trop / outputs['formaldehyde_tropospheric_air_mass_factor']
        return q_ak

    @classmethod
    def cmaq_process(cls, qf, satl3f, key='FORM'):
        """
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
