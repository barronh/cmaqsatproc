__all__ = ['TEMPO_L2', 'TEMPO_L2_NO2', 'TEMPO_L2_HCHO']
from ..core import satellite
import xarray as xr

# GEOS-5 Ap [hPa] for 72 levels (73 edges)
geos5_hyb_a = xr.DataArray(
    [0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01,
     2.609201e+01, 3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01,
     5.805321e+01, 6.436264e+01, 7.062198e+01, 7.883422e+01, 8.909992e+01,
     9.936521e+01, 1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
     1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02, 2.032590e+02,
     2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
     2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02,
     9.236572e+01, 7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01,
     4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01, 1.979160e+01,
     1.645710e+01, 1.364340e+01, 1.127690e+01, 9.292942e+00, 7.619842e+00,
     6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00, 2.620211e+00,
     2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
     6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01,
     1.594950e-01, 1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02,
     3.270000e-02, 2.000000e-02, 1.000000e-02],
    dims=('swt_level_edge',), attrs={'units': 'hPa'}
)
# Bp [unitless] for 72 levels (73 edges)
geos5_hyb_b = xr.DataArray(
    [1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01,
     8.989080e-01, 8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01,
     7.919469e-01, 7.706375e-01, 7.493782e-01, 7.211660e-01, 6.858999e-01,
     6.506349e-01, 6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
     4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01, 2.467411e-01,
     2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
     6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
     0.000000e+00, 0.000000e+00, 0.000000e+00],
    dims=('swt_level_edge',), attrs={'units': '1'}
)

geos5_hybm_a = xr.DataArray(
    (geos5_hyb_a[:-1] + geos5_hyb_a[1:]) / 2,
    dims=('swt_level',), attrs={'units': 'hPa'}
)
geos5_hybm_b = xr.DataArray(
    (geos5_hyb_b[:-1] + geos5_hyb_b[1:]) / 2,
    dims=('swt_level',), attrs={'units': '1'}
)


class TEMPO_L2(satellite):
    __doc__ = """
    Default TEMPO_L2 satellite processor.
    * bbox subsets the mirror_step dimensions
    """
    _defaultkeys = ()

    @classmethod
    def open_dataset(
        cls, path, bbox=None, min_data_quality_flag=2,
        max_eff_cloud_fraction=0.2, max_solar_zenith_angle=70,
        **kwargs
    ):
        """
        Arguments
        ---------
        path : str
            Path to a TEMPO file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        min_data_quality_flag : int
            Minimum acceptable value of main_data_quality_flag:
            - 0=normal
            - 1=suspicious
            - 2=bad
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: TropOMI
            Satellite processing instance
        """
        import xarray as xr
        verbose = kwargs.get('verbose', 0)
        rf = xr.open_dataset(path, group='/')
        sf = xr.open_dataset(path, group='support_data')
        gf = xr.open_dataset(path, group='geolocation')
        pf = xr.open_dataset(path, group='product')
        if verbose > 0:
            print('root:')
            print('  -' + '\n  -'.join(list(rf.data_vars)))
            print('support_data:')
            print('  -' + '\n  -'.join(list(sf.data_vars)))
            print('geolocation:')
            print('  -' + '\n  -'.join(list(gf.data_vars)))
            print('product:')
            print('  -' + '\n  -'.join(list(pf.data_vars)))
        # default suspicious or better
        keep = pf['main_data_quality_flag'] <= min_data_quality_flag
        # within lon/lat box
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            inlat = (gf['latitude'] >= swlat) & (gf['latitude'] <= nelat)
            inlon = (gf['longitude'] >= swlon) & (gf['longitude'] <= nelon)
            keep = keep & inlon & inlat
        # not too cloudy
        keep = keep & (sf['eff_cloud_fraction'] <= max_eff_cloud_fraction)
        # not too dark
        keep = keep & (gf['solar_zenith_angle'] <= max_solar_zenith_angle)
        if not keep.any():
            raise ValueError('No valid pixels')
        # Gather files
        fs = [rf, pf, gf, sf]
        # Only keep mirror_steps that have at least one good value.
        xkeep = keep.any('mirror_step')
        fs = [f.sel(xtrack=xkeep) for f in fs]
        valid = keep.sel(xtrack=xkeep)
        # Mask all bad values
        ds = xr.merge(fs)  # .where(valid)
        # Keep track of valid measure.
        ds['valid'] = valid
        # add corners
        lat_bnds = ds.latitude_bounds
        lon_bnds = ds.longitude_bounds
        corner_slices = {
            'll': 0, 'ul': 1, 'uu': 2, 'lu': 3,
        }
        for cornerkey, corner_slice in corner_slices.items():
            ds[f'{cornerkey}_y'] = lat_bnds.isel(corner=corner_slice)
            ds[f'{cornerkey}_x'] = lon_bnds.isel(corner=corner_slice)
        ds['cn_x'] = ds.longitude
        ds['cn_y'] = ds.latitude
        pressure_mid = ds['surface_pressure'] * geos5_hybm_b + geos5_hybm_a
        pressure_mid.attrs.update(units='hPa', long_name='midlev_pressure')
        ds['midlev_pressure'] = pressure_mid
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        sat.bbox = bbox
        return sat

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
        key = key.replace('amf_cloud_fraction', 'amf_cloud_frac')
        key = key.replace('amf_cloud_pressure', 'amf_cloud_pres')
        key = key.replace('amf_diagnostic_flag', 'amf_diag_flag')
        key = key.replace('amf_stratosphere', 'amf_strat')
        key = key.replace('eff_cloud_fraction', 'eff_cloud_frac')
        key = key.replace('fitted_slant_column', 'fit_slant_col')
        key = key.replace('fitted_slant_column_uncertainty', 'fit_slant_col')
        key = key.replace('ground_pixel_quality_flag', 'pixel_qa_flag')
        key = key.replace('longitude_bounds', 'lon_bnds')
        key = key.replace('latitude_bounds', 'lat_bnds')
        key = key.replace('main_data_quality_flag', 'main_qa_flag')
        key = key.replace('relative_azimuth_angle', 'rel_azimuth')
        key = key.replace('scattering_weights', 'scattering_wts')
        key = key.replace('snow_ice_fraction', 'snow_ice_frac')
        key = key.replace('solar_azimuth_angle', 'solar_azimuth')
        key = key.replace('solar_zenith_angle', 'solar_zenith')
        key = key.replace('surface_pressure', 'surface_pres')
        key = key.replace('temperature_profile', 'temperature')
        key = key.replace('tropopause_pressure', 'tropopause_pres')
        key = key.replace('vertical_column_stratosphere', 'vert_col_strat')
        key = key.replace('vertical_column_total', 'vert_col_tot')
        key = key.replace('vertical_column_total_uncertainty', 'vert_col_tot_u')
        key = key.replace('vertical_column_troposphere', 'vert_col_trop')
        key = key.replace(
            'vertical_column_troposphere_uncertainty', 'vert_col_trop_u'
        )
        key = key.replace('viewing_azimuth_angle', 'viewing_azimuth')
        key = key.replace('viewing_zenith_angle', 'viewing_zenith')
        return key


class TEMPO_NO2_L2(TEMPO_L2):
    __doc__ = """
    Default TEMPO satellite processor.
    * bbox subsets the mirror_step dimensions
    """
    _defaultkeys = (
        'vertical_column_troposphere', 'tropopause_pressure',
        'amf_troposphere', 'amf_stratosphere',
        'scattering_weights', 'gas_profile', 'midlev_pressure'
    )

    @classmethod
    def cmr_links(cls, method='s3', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where concept_id is set to
        "C3685896872-LARC_CLOUD", which is the V04 product. For V03, use
        "C2930725014-LARC_CLOUD".

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download, s3, or OpenDAP
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('concept_id', 'C3685896872-LARC_CLOUD')
        return TEMPO_L2.cmr_links(method=method, **kwargs)


class TEMPO_HCHO_L2(TEMPO_L2):
    __doc__ = """
    Default TEMPO satellite processor.
    * bbox subsets the mirror_step dimensions
    """
    _defaultkeys = (
        'vertical_column_troposphere', 'tropopause_pressure',
        'amf_troposphere', 'amf_stratosphere',
        'scattering_weights', 'gas_profile', 'midlev_pressure'
    )

    @classmethod
    def cmr_links(cls, method='s3', **kwargs):
        """
        Thin wrapper around satellite.cmr_links where concept_id is set to
        "C3685912035-LARC_CLOUD", which is the V04 product. For V03, use
        "C2930730944-LARC_CLOUD"

        Arguments
        ---------
        method : str
            'opendap', 'download', or 's3'.

        Returns
        -------
        links : list
            List of links for download, s3, or OpenDAP
        """
        from copy import copy
        kwargs = copy(kwargs)
        kwargs.setdefault('concept_id', 'C3685912035-LARC_CLOUD')
        return TEMPO_L2.cmr_links(method=method, **kwargs)
