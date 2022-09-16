__all__ = ['OMI']

from ..core import satellite
from ...utils import walk_groups


class OMIL2(satellite):
    __doc__ = """
    Default OMI satellite processor.
    * bbox subsets the nTimes and nTimes_1  dimensions
    """

    @classmethod
    def _open_hierarchical_dataset(cls, path, bbox=None, **kwargs):
        import netCDF4
        import xarray as xr
        tmpf = netCDF4.Dataset(path)
        datakey, = [
            gk for gk in walk_groups(tmpf, '') if gk.endswith('Data Fields')
        ]
        geokey, = [
            gk for gk in walk_groups(tmpf, '')
            if gk.endswith('Geolocation Fields')
        ]
        del tmpf
        datads = xr.open_dataset(path, group=datakey, **kwargs)
        geods = xr.open_dataset(path, group=geokey, **kwargs)

        dimbylen = {
            4: 'nCorners'
        }
        dimbylen[datads.dims['phony_dim_0']] = 'nTimes'
        dimbylen[datads.dims['phony_dim_1']] = 'nXtrack'
        dimbylen[datads.dims['phony_dim_2']] = 'nPresLevels'
        dimbylen[datads.dims['phony_dim_0'] + 1] = 'nTimes_1'
        dimbylen[datads.dims['phony_dim_1'] + 1] = 'nXtrack_1'
        dimbylen[4] = 'nCorners'

        redatadims = {k: dimbylen.get(v, k) for k, v in datads.dims.items()}
        regeodims = {k: dimbylen.get(v, k) for k, v in geods.dims.items()}
        ds = xr.merge([datads.rename(**redatadims), geods.rename(**regeodims)])
        ds = cls.prep_dataset(ds, bbox=bbox, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def prep_dataset(cls, ds, bbox=None, path=None):
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            df = ds[['Latitude', 'Longitude']].to_dataframe()
            df = df.query(
                f'Latitude >= {swlat} and Latitude <= {nelat}'
                + f' and Longitude >= {swlon} and Longitude <= {nelon}'
            )
            times = df.index.get_level_values('nTimes').unique()
            if len(times) < 0:
                raise ValueError(f'{path} has no values in {bbox}')
            ds = ds.isel(
                nTimes=slice(times.min(), times.max() + 1)
            )
            if 'nTimes_1' in ds.dims:
                ds = ds.isel(nTimes_1=slice(times.min(), times.max() + 2))
        return ds

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a OMI OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: OMIL2
            Satellite processing instance
        """
        import xarray as xr

        ds = xr.open_dataset(path, **kwargs).reset_coords()
        if len(ds.dims) == 0:
            return cls._open_hierarchical_dataset(path, bbox=bbox, **kwargs)
        ds = cls.prep_dataset(ds, bbox=bbox, path=path)
        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat

    @classmethod
    def shorten_name(cls, key):
        key = key.replace(
            'ReferenceSectorCorrectedVerticalColumn', 'RefSctCor_VCD'
        )
        key = key.replace('SlantColumnAmount', 'SCD')
        key = key.replace('ColumnAmount', 'VCD')
        key = key.replace('Pressure', 'Press')
        key = key.replace('Scattering', 'Scat')
        key = key.replace('Altitude', 'Alt')
        key = key.replace('Spacecraft', 'Craft')
        key = key.replace('Wavelength', 'WvLen')
        key = key.replace('Registration', 'Reg')
        key = key.replace('Viewing', 'View')
        key = key.replace('Angle', 'Ang')
        key = key.replace('Pixel', 'Pix')
        key = key.replace('Measurement', 'Msrmt')
        key = key.replace('Radiance', 'Rad')
        key = key.replace('Latitude', 'Lat')
        key = key.replace('Longitude', 'Lon')
        key = key.replace('Check', 'Chk')
        key = key.replace('Fraction', 'Frac')
        key = key.replace('Configuration', 'Cfg')
        key = key.replace('Pointer', 'Ptr')
        key = key.replace('CloudRadianceFraction', 'CldRadFrac')
        key = key.replace('QualityFlags', 'QAFlag')
        key = key.replace('TerrainReflectivity', 'TerrainRefl')
        key = key.replace('ClimatologyLevels', 'ClimPresLevels')
        return key


class OMNO2(OMIL2):
    __doc__ = """
    OMNO2 satellite processor.
    * bbox subsets the nTimes and nTimes_1 dimensions
    * valid based three conditions
      * (VcdQualityFlags & 1) == 0
      * XTrackQualityFlags == 0
      * CloudFraction <= 0.3
    """
    _defaultkeys = ('ColumnAmountNO2Trop', 'AmfTrop')

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a TropOMI OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: OMINO2
            Satellite processing instance
        """
        omtmp = OMIL2.open_dataset(path, bbox=bbox, **kwargs)
        ds = omtmp.ds
        corners = {
            'll': 0, 'lu': 3, 'ul': 1, 'uu': 2,
        }
        for key, corner in corners.items():
            ds[f'{key}_x'] = ds['FoV75CornerLongitude'].sel(nCorners=corner)
            ds[f'{key}_y'] = ds['FoV75CornerLatitude'].sel(nCorners=corner)

        ds['valid'] = (
            ((ds['VcdQualityFlags'].astype('i') & 1) == 0)
            & (ds['XTrackQualityFlags'] == 0)
            & (ds['CloudFraction'] <= 0.3)
        )
        if not ds['valid'].any():
            raise ValueError(f'No valid pixels in {path} with {bbox}')

        sat = cls()
        sat.path = path
        sat.ds = ds
        sat.bbox = bbox
        return sat


class OMHCHO(OMIL2):
    __doc__ = """
    OMHCHO satellite processor.
    * bbox subsets the nTimes and nTimes_1 dimensions
    * valid based three conditions
      * MainDataQualityFlag == 0
      * (XtrackQualityFlagsExpanded & 1) == 0
      * AMFCloudFraction <= 0.3
    """
    _defaultkeys = (
        'ColumnAmount', 'ReferenceSectorCorrectedVerticalColumn',
        'AirMassFactor'
    )

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Arguments
        ---------
        path : str
            Path to a TropOMI OpenDAP-style file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: OMI
            Satellite processing instance
        """
        import xarray as xr
        omtmp = OMIL2.open_dataset(path, bbox=bbox, **kwargs)
        ds = omtmp.ds
        ds['valid'] = (
            (ds['MainDataQualityFlag'] == 0)
            & ((ds['XtrackQualityFlagsExpanded'].astype('i') & 1) == 0)
            & (ds['AMFCloudFraction'] <= 0.3)
        )
        corners = {
            'll': (slice(None, -1), slice(None, -1)),
            'ul': (slice(None, -1), slice(1, None)),
            'lu': (slice(1, None), slice(None, -1)),
            'uu': (slice(1, None), slice(1, None)),
        }
        for key, (tslice, xslice) in corners.items():
            ds[f'{key}_x'] = xr.DataArray(
                ds['PixelCornerLongitudes'].isel(
                    nTimes_1=tslice, nXtrack_1=xslice
                ).transpose('nTimes_1', 'nXtrack_1').values,
                dims=('nTimes', 'nXtrack')
            )
            ds[f'{key}_y'] = xr.DataArray(
                ds['PixelCornerLatitudes'].isel(
                    nTimes_1=tslice, nXtrack_1=xslice
                ).transpose('nTimes_1', 'nXtrack_1').values,
                dims=('nTimes', 'nXtrack')
            )

        if not ds['valid'].any():
            raise ValueError(f'No valid pixels in {path} with {bbox}')
        sat = cls()
        sat.ds = ds
        sat.bbox = bbox
        return sat
