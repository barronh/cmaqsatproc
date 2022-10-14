__all__ = ['goes_aod']

from ..core import satellite
from ...utils import EasyDataFramePolygon


class goes_aod(satellite):
    __doc__ = """
    goes_aod processor
    * valid if DQF < dqflt (default = 1)
    * pixel corners are interpolated in projected space.

    Note 1: Qualities are 0: high; 1: medium; 2: low; 3 not retrieved.
    Note 2: The default requires the highest quality only. This is based
            on experience. Because we are gridding, high and medium would get
            spatially mixed. This makes it hard to see if a monthly average is
            a bunch of isolated medium quality pixels with no repeat
            measurements or many of the same pixel. Thus, the higher quality
            requirement.
    """
    _defaultkeys = ('AOD', 'DQF')

    @classmethod
    def s3_to_level3(
        cls, date, satkey, grid, griddims=None, weighting='area', bbox=None,
        verbose=0, varkeys=None, as_dataset=True, link_kwargs=None, **kwargs
    ):
        """
        Thin wrapper around s3_links, open_datasets, and to_level3.

        For keyword documentation, see those functions.
        """
        from copy import copy
        if link_kwargs is None:
            link_kwargs = {}
        link_kwargs = copy(link_kwargs)
        link_kwargs.setdefault('resolution', 'H')
        if varkeys is None:
            varkeys = ()
        if bbox is None:
            grid.unary_union.envelope
        links = cls.s3_links(
            date, satkey, **link_kwargs
        )
        if verbose > 1:
            print(len(links), links)
        elif verbose > 0:
            print(len(links), 'first link', links[0])
        sat = cls.open_datasets(links, bbox=bbox, **kwargs)

        return sat.to_level3(
            *varkeys, grid=grid, griddims=griddims, weighting=weighting,
            as_dataset=as_dataset, verbose=verbose
        )

    @classmethod
    def s3_links(cls, date, satkey, product='ABI-L2-AODC', resolution='H'):
        """
        Return s3 links for date

        Arguments
        ---------
        date : str or date-like
            pandas.to_datetime will convert this into a date object.
        satkey : str
            goes16, goes17, or goes18... any goes satellite that has product
        product : str
            Usually ABI-L2-AODC, but could be ABI-L2-AODF
        resolution : str
            Choose, H, d, m or Y.

        Returns
        -------
        links : list
            All links within the folder based on resolution:
            - 'H': s3://noaa-{satkey}/ABI-L2-AODC/{date:%Y/%j/%H}/
            - 'd': s3://noaa-{satkey}/ABI-L2-AODC/{date:%Y/%j}/
            - 'Y': s3://noaa-{satkey}/ABI-L2-AODC/{date:%Y}/
            - 'm': iterates on all julian days in month
        """
        import s3fs
        import pandas as pd
        date = pd.to_datetime(date)
        if resolution == 'H':
            pattern = f's3://noaa-{satkey}/{product}/{date:%Y/%j/%H}/'
        elif resolution == 'd':
            pattern = f's3://noaa-{satkey}/{product}/{date:%Y/%j}/'
        elif resolution == 'Y':
            pattern = f's3://noaa-{satkey}/{product}/{date:%Y}/'
        elif resolution == 'm':
            # Use replace the day with 1 to get the start of the month. Then,
            # add the MonthEnd offset to get the end of the month. Next,
            # create a daily range for all days of the month. Finally, query
            # the s3 bucket for each day.
            sdate = date.replace(day=1)
            edate = sdate + pd.offsets.MonthEnd()
            dates = pd.date_range(sdate, edate, freq='d')
            paths = []
            for date in dates:
                dpaths = cls.s3_links(
                    date, satkey, product=product, resolution='d'
                )
                paths.extend(dpaths)
            return paths
        else:
            raise KeyError(
                f'resoluton must be H, d, m, or Y; got {resolution}'
            )

        fs = s3fs.S3FileSystem(anon=True)
        paths = fs.find(pattern)
        paths = ['s3://' + path for path in paths]

        return paths

    @classmethod
    def cmr_links(cls, *args, **kwds):
        raise NotImplementedError('Not implemented for GOES')

    @property
    def _defgeofunc(self):
        "Not wrapping dateline, because the data is in projected units"
        return lambda df: EasyDataFramePolygon(df, wrap=False)

    @classmethod
    def open_datasets(cls, paths, bbox, dqflt=1):
        """
        Similar to open_dataset, but assumes that all data is on a single grid.
        Therefore, it can mask pixels and average them before creating a "net"
        dataset.

        Basically,
        ds = xr.concat([open_dataset(path, bbox) for path in paths], 'scan')
        ds = ds.where(ds['valid']).mean('scan')

        Only variables valid, AOD, DQF, and goes_imager_projection are kept.

        Arguments
        ---------
        paths : list
            List of paths
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees

        Returns
        -------
        sat : satellite
        """
        import xarray as xr
        dss = []
        for path in paths:
            sat = cls.open_dataset(path, bbox, dqflt=dqflt)
            dss.append(
                sat.ds[['valid', 'AOD', 'DQF', 'goes_imager_projection']]
            )
        ds = xr.concat(dss, dim='scan')
        ds = ds.where(ds['valid']).mean('scan')
        ds.attrs.update(dss[0].attrs)
        ds['goes_imager_projection'] = dss[0]['goes_imager_projection']
        cls.prep_dataset(ds, bbox=bbox, dqflt=dqflt)
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        sat.bbox = bbox
        sat._crs = ds.attrs['crs']

        return sat

    @classmethod
    def open_dataset(cls, path, bbox=None, dqflt=1, **kwargs):
        """
        Open a GOES AOD dataset for satellite processing.

        Arguments
        ---------
        path : str
            Path to dataset
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
        dqflt : float
            Only AOD with data quality flags less than dqflt are considered
            valid
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat : satellite
            Satellite processing object.
        """
        import xarray as xr
        import copy
        kwargs = copy.copy(kwargs)
        kwargs.setdefault('decode_times', False)
        if 's3:' in path:
            import s3fs
            fs = s3fs.S3FileSystem(anon=True)
            fileObj = fs.open(path)
            ds = xr.open_dataset(fileObj, **kwargs)
        else:
            ds = xr.open_dataset(path, **kwargs)
        cls.prep_dataset(ds, bbox, dqflt=dqflt)
        sat = cls()
        sat.path = path
        sat.ds = ds.reset_coords()
        sat.bbox = bbox
        sat._crs = ds.attrs['crs']

        return sat

    @classmethod
    def prep_dataset(cls, ds, bbox=None, dqflt=1):
        """
        Prepare the dataset by adding a projection, applying valid spatial
        checks and requiring that the DQF variable have a value less than
        dqflt.
        """
        import numpy as np
        import pyproj

        ds['valid'] = ds['DQF'] < dqflt
        proj_info = ds['goes_imager_projection']
        lon_origin = float(proj_info.attrs['longitude_of_projection_origin'])
        H = float(proj_info.attrs['perspective_point_height'])
        r_eq = float(proj_info.attrs['semi_major_axis'])
        r_pol = float(proj_info.attrs['semi_minor_axis'])
        ds.attrs['crs'] = (
            f'+proj=geos +sweep=x +lon_0={lon_origin} +h={H} +x_0=0 +y_0=0'
            + f' +a={r_eq} +b={r_pol} +no_defs'
        )
        x_1d = ds['x'][:] * H
        y_1d = ds['y'][:] * H
        x, y = np.meshgrid(x_1d, y_1d)
        lon, lat = pyproj.Proj(ds.attrs['crs'])(x, y, inverse=True)
        ds['lon'] = (('y', 'x'), lon)
        ds['lat'] = (('y', 'x'), lat)
        x = ds['x'] * H
        y = ds['y'] * H
        xe = np.concatenate([
            x[:1], (x.values[:-1] + x.values[1:]) / 2, x[-1:]
        ])
        ye = np.concatenate([
            y[:1], (y.values[:-1] + y.values[1:]) / 2, y[-1:]
        ])
        l_y = ye[:-1]
        u_y = ye[1:]
        l_x = xe[:-1]
        u_x = xe[1:]
        ds['lu_x'] = ds['ll_x'] = (('x',), l_x)
        ds['uu_x'] = ds['ul_x'] = (('x',), u_x)
        ds['ul_y'] = ds['ll_y'] = (('y',), l_y)
        ds['uu_y'] = ds['lu_y'] = (('y',), u_y)
        ds['cn_x'] = x
        ds['cn_y'] = y
        if bbox is not None:
            swlon, swlat, nelon, nelat = bbox
            ds['valid'] = (
                ds['valid']
                & (lon >= swlon) & (lon <= nelon)
                & (lat >= swlat) & (lat <= nelat)
            )
