__all__ = ['TEMPO_L2', 'TEMPO_NO2_L2', 'TEMPO_HCHO_L2']
from ..core import satellite
import xarray as xr

# GEOS-5 Ap [hPa] for 72 levels (73 edges)
geos5_hybi_a = xr.DataArray(
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
geos5_hybi_b = xr.DataArray(
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
    (geos5_hybi_a[:-1] + geos5_hybi_a[1:]) / 2,
    dims=('swt_level',), attrs={'units': 'hPa'}
)
geos5_hybm_b = xr.DataArray(
    (geos5_hybi_b[:-1] + geos5_hybi_b[1:]) / 2,
    dims=('swt_level',), attrs={'units': '1'}
)


def _cmaq2tempo(
    tf, qf, method='nearest', pkey=None, crs=None, droppvar=False
):
    """
    Interpolate CMAQ to TEMPO (horizontally and/or vertically).

    Arguments
    ---------
    tf : xr.Dataset
        TEMPO dataset which must have:
        - longitude/latitude for horizontal interpolation
        - either midlev_pressure or surface_pressure for vertical interpolation
    qf : xr.Dataset or xr.DataArray
        CMAQ dataset or DataArray. Must have crs as an attribute if crs is not
        provided as a keyword
    method : str
        xarray sel method keyword.
    pkey : str
        If provided, pkey is pressure in Pa and will be used to vertically
        interpolate CMAQ to the TEMPO swt_level.
    crs : str
        Optional, overrides qf.attrs['crs']
    droppvar : bool
        If True, do not interpolate or output coordinate variable.
        If False, keep the PRES variable in the output.

    Returns
    -------
    tqf : xr.Dataset or xr.DataArray
        Same type as input with dimensions converted from (TSTEP, LAY, ROW,
        COL) to (mirror_step, xtrack, LAY)
    """
    qhq = ('COL' in qf.dims and 'ROW' in qf.dims)
    thq = ('COL' in tf.dims and 'ROW' in tf.dims)
    qvq = 'LAY' in qf.dims
    tvq = 'LAY' in tf.dims
    if qhq and not thq:
        tqf = _cmaq2tempoh(tf, qf, method=method, crs=crs)
    else:
        tqf = qf
    if pkey is not None and (qvq and not tvq):
        PRESPa = tqf[pkey]
        if droppvar:
            tqf = tqf.drop_vars(pkey)
        ptqf = _cmaq2tempov(tf, tqf, PRESPa)
    else:
        ptqf = tqf

    return ptqf


def _cmaq2tempoh(tf, qf, method='nearest', crs=None):
    """
    Interpolate CMAQ to TEMPO (horizontally and/or vertically).

    Arguments
    ---------
    tf : xr.Dataset
        TEMPO dataset which must have:
        - longitude/latitude for horizontal interpolation
    qf : xr.Dataset or xr.DataArray
        CMAQ dataset or DataArray. Must have crs as an attribute if crs is not
        provided as a keyword
    method : str
        xarray sel method keyword.
    crs : str
        Optional, overrides qf.attrs['crs']

    Returns
    -------
    tqf : xr.Dataset or xr.DataArray
        Same type as input with dimensions converted from (TSTEP, LAY, ROW,
        COL) to (mirror_step, xtrack, LAY)
    """
    import xarray as xr
    import pyproj
    if crs is None:
        crs = qf.attrs['crs']

    # Regrid CMAQ to mirror_step and xtrack (from TSTEP, LAY, ROW, COL)
    lon = tf['longitude']
    lat = tf['latitude']
    x, y = pyproj.Proj(crs)(lon, lat)
    x = xr.DataArray(x, dims=lon.dims)
    y = xr.DataArray(y, dims=lat.dims)
    tidx, ridx, cidx = xr.broadcast(tf['time'], y, x)
    tqf = qf.sel(TSTEP=tidx, ROW=ridx, COL=cidx, method=method)
    return tqf


def _cmaq2tempov(tf, qf, pvar):
    """
    Interpolate CMAQ to TEMPO (vertically).

    Arguments
    ---------
    tf : xr.Dataset
        TEMPO dataset which must have:
        - either midlev_pressure or surface_pressure for vertical interpolation
    qf : xr.Dataset or xr.DataArray
        CMAQ dataset or DataArray. Must have crs as an attribute if crs is not
        provided as a keyword
    pvar : xr.DataArray
        Mid-level pressure in Pa and will be used to vertically interpolate
        CMAQ to the TEMPO swt_level.

    Returns
    -------
    tqf : xr.Dataset or xr.DataArray
        Same type as input with dimension LAY converted to swt_level:
        - (mirror_step, xtrack, LAY) to (mirror_step, xtrack, level)
        - (ROW, COL, LAY) to (ROW, COL, level)
    """
    from ...utils import coord_interp
    PREShPa = pvar / 100.
    pres = _get_midp(tf)
    # interpolate to the TEMPO pressure coordinate.
    # ascending=True because the PRES [1000, ..., 50] variable is sorted
    # to [50, ..., 1000] using desceding LAY [0.999, ..., 0] coord
    tqf = coord_interp(
        pres, PREShPa, qf,
        indim='LAY', outdim='swt_level', ascending=True
    )
    return tqf


def _get_midp(ds, add=False):
    if 'midlev_pressure' in ds:
        pres = ds['midlev_pressure']
    else:
        psfc = ds['surface_pressure']
        pres = psfc * geos5_hybm_b + geos5_hybm_a
        if add:
            ds['midlev_pressure'] = pres
    return pres


def _get_delp(ds, add=False):
    if 'delta_pressure' in ds:
        delpres = ds['delta_pressure']
    else:
        psfc = ds['surface_pressure']
        bigpres = psfc * geos5_hybi_b[:-1] + geos5_hybi_a[:-1]
        lilpres = psfc * geos5_hybi_b[1:] + geos5_hybi_a[1:]
        outdims = psfc.dims + ('swt_level',)
        delpres = xr.DataArray(
            bigpres.data - lilpres.data,
            dims=outdims, attrs={'units': 'hPa'}
        )
        if add:
            ds['delta_pressure'] = delpres

    return delpres


def _get_sw(ds, TA):
    """
    Return temperature adjusted scattering weights (w') at PRES levels
    using the TEMPO L2 Trace Gas and Cloud Level 2 and 3 Data Products:
    User Guide Version 1.1 See pg 19. Where w_z is scattering_weights (sw),
    omega_z is the partical vertical column (pvcd), and alpha is calculated
    from temperature (TA)

    AMF = sum_z(w_z * alpha_z * omega_z) / sum_z(omega_z)
        = sum_z(w'_z          * omega_z) / sum_z(omega_z)

    alpha_z = 1 - 0.00316 * (T_z - T_0) + 3.39e-6 * (T_z - T_0)**2

    Arguments
    ---------
    TA : xr.DataArray
        Mid-level temperature in Kelvin regridded to l2 (mirror_step,
        xtrack, LAY)

    Returns
    -------
    sw : xr.DataArray
        Temperature adjusted scattering weights.
    """
    # Zero-scattering weights are not valid.
    sw = ds['scattering_weights'].where(lambda x: x > 0)  # unitless
    a = 0.00316
    b = 3.39e-6
    T0 = 220.
    dT = TA - T0
    alpha = 1 - a * dT + b * dT**2
    tasw = sw * alpha
    tasw.attrs.update(sw.attrs)
    tasw.attrs['description'] = 'temperature corrected scattering_weights'
    return tasw


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
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        _get_midp(sat.ds, add=True)
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

    def add_cmaq_amf(
        self, qf, qkey, extent='troposphere', amfkey=None,
        swopt='tempo', debug=False
    ):
        """
        Add an AMF using an alternate shape factor from the partial vertical
        column density (pvcd) calculated from CMAQ met and conc data
        interpolated to the TEMPO grid. This should be done before `to_level3`

        Arguments
        ---------
        qf : xr.Dataset
            Must contain DENS, TA, PRES, and qkey. If DENS is not present, then
            PRES/TA/Q can be used. If Q is missing, just PRES and TA.
        qkey : str
            Model variable to be used for shape factor.
        extent : str
            troposphere, stratosphere, total
        amfkey : str
            Name for output variable to be stored. Defaults to amf_cmaq_extent
            where extent is troposphere, stratosphere, or total
        swopt : str
            The scattering weight option determines if the scattering weight is
            interpolated to CMAQ or if CMAQ is interpolated to tempo. Options:
            'tempo' calculates the AMF on native TEMPO coordinate. 'cmaq'
            calculates the AMF on native CMAQ coordinate

        Returns
        -------
        amf : xr.DataArray
            Air mass factor with temperature corrected scattering weights
            applied to the shape factor from pvcd.

        Notes
        -----
        CMAQ met/conc data is often unreliable in the stratosphere. It is best
        to only use the troposphere.
        """
        from ...utils import coord_interp

        if amfkey is None:
            amfkey = f'amf_cmaq_{extent}'

        # Subset just the necessary variables for AMF calculation
        # - PRES for vertical interpolation
        # - TA for scattering weight correction.
        # - qkey for shape factor
        sqf = qf[['TA', 'PRES', qkey]]
        tf = self.ds
        # Pair in space and time with TEMPO using vertical interpolation
        if swopt == 'tempo':
            # Inteprolate CMAQ to TEMPO both horizontally and vertically
            tqf = _cmaq2tempo(tf, sqf, pkey='PRES', crs=qf.crs, droppvar=False)
            # calculate air area_density for later conversion of ppm to mole/m2
            delp_pa = _get_delp(tf) * 100.
            # P * A = m * a = n * M * g
            # P/Mw/g = n/A = delp_pa / MWAIR / GRAV
            #       [=] Pa * mol / kg * s**2 /m
            #       [=] kg / m / s**2 * mol / kg * s**2 / m
            #       [=] mol / m2
            #
            # copied from https://github.com/USEPA/CMAQ/blob/main/
            # CCTM/src/ICL/fixed/const/CONST.EXT
            MWAIR = 0.0289628
            GRAV = 9.80622
            tqf['mole_per_m2'] = delp_pa / MWAIR / GRAV
            taqf = tqf
            laykey = 'swt_level'
        elif swopt == 'cmaq':
            # calculate air area_density for later conversion of ppm to mole/m2
            sqf['mole_per_m2'] = qf.csp.mole_per_m2()
            # Interpolate NO2, TA, and PRES only horizontally
            tqf = _cmaq2tempo(tf, sqf, crs=qf.crs, droppvar=False)
            # Interpolate TA horizontally and vertically
            # - PRES will be used to interpolate Ta and then dropped
            # - TA used for scattering weight temperature correction
            #   on the TEMPO grid before interpolation to the CMAQ grid
            taqf = _cmaq2tempo(
                tf, tqf[['TA', 'PRES']], pkey='PRES', droppvar=True
            )
            laykey = 'LAY'
        else:
            raise KeyError(f'swopt must be tempo or cmaq; got {swopt}')

        # Get temperature corrected TEMPO scattering weights (sw_z * alpha_z)
        rsw = _get_sw(tf, taqf.TA)
        # Note for future:
        #
        # The scattering_weights variable is zero below clouds and masked
        # within _get_sw. The masked elements of gas_profile column
        # is measured only above the zero scattering weights. Regridding
        # magnifies challenges implicit in gridded data. The cloud is part way
        # between a layer. In that case, the sensitivity to the layer is
        # non-zero for the part above the cloud, and not-applicable for the
        # part under the cloud. There are options that can be illustrated by a
        # scattering weight vector with cloud masking (swc), without cloud
        # masking, and an opaque cloud at 825hPa.
        # profile:
        #                   cloud
        #                     |
        #                     v
        # P [hPa] = [990, 850, 700, 200, 0.015]
        # swc [1] = [0.0, 0.0, 2.0, 2.6, 2.5]
        # sws [1] = [0.3, 1.0, 2.0, 2.6, 2.5]
        #
        # and an interpolated sw at destination coordinates below:
        #              cloud
        #                |
        #                v
        # P [hPa] = [990, 800, 450]
        #
        # 1. Masking zeros before interpolating:
        #    swi [1] = [nan, nan, 2.3]
        # 2. Masking zeros after interpolating:
        #    swi [1] = [nan, 0.66, 2.3]
        # 3. Interpolating sws then masking (not an option because already 0)
        #    swi [1] = [nan, 1.33, 2.3]
        # 4. Integrate pvcd to TEMPO coordinate.
        #    Have to assume uniform distribution within cell.
        #
        # For now, interpolating after masking
        if swopt == 'cmaq':
            PREShPa = tqf.PRES / 100.
            pres = _get_midp(tf)
            # ascending=False because the midlev_pressure [1000, ..., 50]
            # variable is sorted [1000, ..., 50] using ascending swt_level
            # [0, ..., 72] coord
            rsw = coord_interp(
                PREShPa, pres, rsw,
                indim='swt_level', outdim='LAY', ascending=False
            )

        # Add additional filters for atmospheric focus
        if extent == 'troposphere':
            tppres = tf['tropopause_pressure']
            validsw = (~rsw.isnull()) & (tqf.PRES / 100 > tppres)
        elif extent == 'stratosphere':
            tppres = tf['tropopause_pressure']
            validsw = (~rsw.isnull()) & (tqf.PRES / 100 < tppres)
        elif extent == 'total':
            validsw = (~rsw.isnull())
        else:
            msg = f'extent can be troposphere/stratosphere/total; got {extent}'
            raise KeyError(msg)

        # remove unseen portions of the column
        pvcd = tqf[qkey].where(validsw) * 1e-6 * tqf['mole_per_m2']

        # AMF = sum_z(sw_z * alpha_z * omega_z) / sum_z(omega_z)
        amf = (pvcd * rsw).sum(laykey) / pvcd.sum(laykey)
        amf = amf.where(validsw.any(laykey))
        amf.attrs.update(
            long_name=amfkey, units='1',
            var_desc='air mass factor from alternative pvcd'
        )

        # add variable to dataset
        tf[amfkey] = amf.reset_coords(drop=True)
        if debug:
            for k in tqf.data_vars:
                tf['cmaq_' + k] = tqf[k]
            tf['cmaq_gas_profile'] = pvcd
        # if not already default output, add to default output
        if amfkey not in self._defaultkeys:
            self._defaultkeys = self._defaultkeys + (amfkey,)

        return amf

    @classmethod
    def cmaq_process(cls, qf, l3, qkey, l3key, extent='troposphere'):
        """
        Process CMAQ as though it were observed by TEMPO, which is simply
        based on the overpass time.

        Arguments
        ---------
        qf : xarray.Dataset
            CMAQ file that has composition (e.g., NO2) and met (DENS and TA or
            PRES/TA/Q). qf should already have the UTC time selected that
            matches l3.
        l3 : xarray.Dataset
            Output from to_level3, paths_to_level3, or cmr_to_level3 with
            as_dataset=True (the default).
        qkey : str
            CMAQ Key NO2 or HCHO.
        l3key : str
            Satellite key

        Returns
        -------
        overf : xr.DataArray
            An overpass file with satellite-like CMAQ.
        """
        from ...utils import coord_interp
        tpkey = 'tropopause_pressure'
        overf = qf.where(~l3[l3key].isnull())
        n_per_m2 = overf.csp.mole_per_m2(add=True)
        tgtvar = overf[qkey]
        tgtunit = tgtvar.units.strip().lower()
        if tgtunit.startswith('ppm'):
            pp2vmr = 1e-6
        elif tgtunit.startswith('ppb'):
            pp2vmr = 1e-9
        elif tgtunit.startswith('ppt'):
            pp2vmr = 1e-12

        nden = n_per_m2 * tgtvar * pp2vmr  # moles-air/m2 * moles-i/moles-air
        # interpolate TA from LAY to swt_level
        TAK = _cmaq2tempov(l3, qf['TA'], pvar=qf['PRES'])
        sw = _get_sw(l3, TAK)
        pmid = _get_midp(l3)
        # interpolate sw to LAY
        # ascending=False because the midlev_pressure [1000, ..., 50]
        # variable is sorted [1000, ..., 50] using ascending swt_level
        # [0, ..., 72] coord
        PREShPa = qf.PRES / 100.
        sw = coord_interp(
            PREShPa, pmid, sw,
            indim='swt_level', outdim='LAY', ascending=False
        )
        if extent == 'troposphere':
            ismine = PREShPa > l3[tpkey]
        elif extent == 'stratoposphere':
            ismine = PREShPa < l3[tpkey]
        elif extent == 'total':
            ismine = PREShPa > -999
        else:
            msg = f'extent can be troposphere/stratosphere/total; got {extent}'
            raise KeyError(msg)
        ismine = ismine & ~sw.isnull()
        nden = nden.where(ismine)
        sw = sw.where(ismine)
        vcd = nden.sum('LAY')
        amf = (sw * nden).sum('LAY') / vcd
        overf[f'{qkey}_{extent}'] = vcd
        overf[f'amf_{extent}'] = amf
        return overf


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

    @classmethod
    def cmaq_process(cls, qf, l3, qkey=None, l3key=None, extent=None):
        if qkey is None:
            qkey = 'NO2'
        if extent is None:
            extent = 'troposphere'
        if l3key is None:
            l3key = f'vertical_column_{extent}'

        opts = dict(qkey=qkey, l3key=l3key, extent=extent)
        return TEMPO_L2.cmaq_process(qf, l3, **opts)


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

    @classmethod
    def cmaq_process(cls, qf, l3, qkey=None, l3key=None, extent=None):
        if qkey is None:
            qkey = 'HCHO'
        if extent is None:
            extent = 'troposphere'
        if l3key is None:
            l3key = f'vertical_column_{extent}'

        opts = dict(qkey=qkey, l3key=l3key, extent=extent)
        return TEMPO_L2.cmaq_process(qf, l3, **opts)
