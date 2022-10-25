__all__ = ['satellite']

from ..utils import EasyDataFramePolygon, grouped_weighted_avg, rootremover
from ..utils import csp_version


class satellite:
    _crs = 4326
    _defaultkeys = ()
    _validkey = 'valid'
    _geokeys = (
        'll_x', 'll_y',
        'lu_x', 'lu_y',
        'uu_x', 'uu_y',
        'ul_x', 'ul_y',
        'cn_x', 'cn_y'
    )

    @property
    def _defgeofunc(self):
        """
        Function that takes a dataframe and returns geometries for each row.
        """
        return EasyDataFramePolygon

    @classmethod
    def cmr_links(cls, method='opendap', **kwds):
        """
        Use utils.getcmrlinks to get links from the NASA Common Metadata Repo

        Arguments
        ---------
        method : str
            Options are opendap, download, or s3
        kwds : mappable
            Passed through to utils.getcmrlinks. See getcmrlinks for valid
            keywords

        Returns
        -------
        links : list
            List of links to OpenDAP, http downloadable, or s3 links
        """
        from copy import copy
        from ..utils import getcmrlinks
        kwds = copy(kwds)
        down_f = (
            lambda x: (
                'opendap' not in x['href']
                and (
                    x['href'].endswith('he5')
                    or x['href'].endswith('nc')
                    or x['href'].endswith('h5')
                    or x['href'].endswith('hdf')
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
                    or x['href'].endswith('h5')
                    or x['href'].endswith('hdf')
                )
                and x['href'].startswith('s3')
            )
        )
        open_f = (
            lambda x: 'opendap' in x['href']
            and (
                not x['href'].endswith('html')
                or x['href'].endswith('.nc.html')
                or x['href'].endswith('.he5.html')
                or x['href'].endswith('.h5.html')
                or x['href'].endswith('.hdf.html')
            )
        )
        kwds.setdefault(
            'filterfunc', {
                'opendap': open_f, 's3': s3_f, 'download': down_f
            }[method]
        )

        rawlinks = sorted(getcmrlinks(**kwds))
        # MODIS requires extra processing because OpenDAP links from CMR are
        # to the html interface instead of to the file itself.
        #
        # Conceptually, this may apply to any/all datasets
        if method != 'opendap':
            links = rawlinks
        else:
            links = []
            for link in rawlinks:
                if link.endswith('.html'):
                    link = link[:-5]
                links.append(link)

        return links

    @classmethod
    def from_dataset(cls, ds, path='unknown'):
        """
        Create a satellite object from a dataset

        Arguments
        ---------
        ds : xarray.Dataset
            Satellite dataset

        Returns
        -------
        sat: satellite
            Satellite processing instance
        """
        sat = cls()
        sat.ds = ds
        sat.path = path
        sat.bbox = None
        return sat

    @classmethod
    def open_dataset(cls, path, bbox=None, **kwargs):
        """
        Create a satellite object from the dataset at a path. The path can be
        a local path or a remote path (OpenDAP or s3).

        Arguments
        ---------
        path : str
            Path to a satellite file
        bbox : iterable
            swlon, swlat, nelon, nelat in decimal degrees East and North
            of 0, 0
        kwargs : mappable
            Passed to xarray.open_dataset

        Returns
        -------
        sat: satellite
            Satellite processing instance
        """
        import xarray as xr
        sat = cls()
        sat.path = path
        if 's3:' in path:
            import s3fs
            fs = s3fs.S3FileSystem(anon=True)
            with fs.open(path) as fileObj:
                ds = xr.open_dataset(fileObj)
        else:
            ds = xr.open_dataset(path, *kwargs)
        sat.ds = ds
        if bbox is not None:
            import warnings
            warnings.warn(f'{cls} bbox not implemented; all cells returned')
        return sat

    @property
    def ds(self):
        return self._ds

    @ds.setter
    def ds(self, ds):
        self._ds = ds

    def to_dataframe(self, *varkeys, valid=True, geo=False, default_keys=False):
        """
        Transform the Dataset into a dataframe. This function is similar to
        xr.Dataset.to_dataframe, but includes geopandas support and filtering
        for valid pixels.

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
        if default_keys:
            varkeys = list(varkeys) + list(self._defaultkeys)
        validdims = self.ds[self._validkey].dims
        usevalid = len(varkeys) == 0
        for varkey in varkeys:
            for dk in self.ds[varkey].dims:
                if dk in validdims:
                    usevalid = True

        if usevalid and valid:
            keys = [self._validkey] + list(varkeys)
        else:
            keys = list(varkeys)

        df = self.ds[keys].to_dataframe()
        if valid and usevalid:
            df = df.query('valid == True')
        if not (geo is False):
            if hasattr(self, '_geodf'):
                gdf = self._geodf
            else:
                import geopandas as gpd
                if geo is True:
                    geo = self._defgeofunc
                gdf = self.to_dataframe(*self._geokeys, valid=valid)
                if gdf.shape[0] == 0:
                    raise ValueError('No valid pixels')
                gdf = gpd.GeoDataFrame(
                    gdf, geometry=geo(gdf), crs=self._crs
                ).drop(list(self._geokeys), axis=1)
                self._geodf = gdf
            df = gdf[['geometry']].join(df)
        if 'valid' in df.columns:
            if 'valid' not in varkeys:
                df = df.drop('valid', axis=1)
        for key in df.columns:
            if key in self.ds.variables:
                df[key].attrs.update(self.ds[key].attrs)
        return df

    @classmethod
    def add_weights(cls, intx, option='equal'):
        """
        Add Weights to intersection dataframe (intx) using predefined options:
        * 'equal' uses equal weights for all overlapping pixels
        * 'area' uses area of intersection for all overlapping pixels
        * other methods can be added by overriding this mehtod

        Arguments
        ---------
        intx : geopandas.GeoDataFrame
            Dataframe that has a geometry field that supports area

        Returns
        -------
        None
        """
        if option == 'equal':
            intx['weight'] = 1
        elif option == 'area':
            intx['weight'] = intx.geometry.area
            if (intx['weight'] == 0).all():
                import warnings
                warnings.warn(
                    'All areas are zero, which likely means that geometries'
                    + ' were points. Reverting to equal weights.'
                )
                intx['weight'] = 1
        else:
            raise KeyError(f'Unknown weighting option {option}')

        return None

    def shorten_name(self, k):
        """
        By default, nothing is done and k is returned. Override this method to
        make short names that are suitable for other purposes (e.g., IOAPI 16
        character names)
        """
        return k

    def to_level3(
        self, *varkeys, grid, griddims=None, weighting='area',
        as_dataset=True, verbose=0
    ):
        """
        Convert variables from L2 file to a custom L3 file.

        Arguments
        ---------
        varkeys : iterable
            See to_dataframe
        grid : geopandas.GeoDataFrame
            Defines the grid used as the L3 destination
        griddims : iterable
            Defaults to grid.index.names
        weighting : str
            Passed as option to self.add_weights

        Returns
        -------
        outputs : dict or xr.Dataset
            Dictionary of outputs by output dimensions
            or Dataset of outputs as a dataset
        """
        import pandas as pd
        import geopandas as gpd
        import time
        if len(varkeys) == 0:
            varkeys = [
                k for k in self._defaultkeys if k in list(self.ds.data_vars)
            ]
            if verbose > 0:
                print('defaults', varkeys)
        if isinstance(grid, gpd.GeoSeries):
            grid = gpd.GeoDataFrame(dict(), geometry=grid)
        if griddims is None:
            griddims = list(grid.index.names)
        if verbose > 1:
            print('griddims', griddims)
        dimsets = {}
        for key in varkeys:
            dims = self.ds[key].dims
            dimsets.setdefault(dims, []).append(key)
        if verbose > 1:
            print('dimsets', dimsets)
            t0 = time.time()

        if verbose > 0:
            print('Making geometry', flush=True)

        if verbose > 1:
            t0 = time.time()
        # Index is lost during overlay, and must be recreated
        geodf = self.to_dataframe(geo=True).to_crs(grid.crs)[['geometry']]
        geodf['source_area'] = geodf['geometry'].area
        # possible_matches_index = grid.sindex.intersection(geodf.total_bounds)
        # igrid = grid.iloc[possible_matches_index]
        if geodf.shape[0] == 0:
            raise ValueError('No valid pixels')
        if verbose > 1:
            print(
                f'Geometry has {geodf.shape[0]} rows'
                + f' ({time.time() - t0:.1f}s)'
            )
            t0 = time.time()
        myidx = list(geodf.index.names) + list(grid.index.names)
        if verbose > 0:
            print('Making overlay', flush=True)

        intx = gpd.overlay(
            geodf.reset_index(), grid.reset_index(), how='intersection',
            keep_geom_type=False, make_valid=True
        )
        intx.set_index(myidx, inplace=True)

        if verbose > 1:
            print(
                f'Overlay has {intx.shape[0]} rows '
                + f'({time.time() - t0:.1f}s)'
            )
            t0 = time.time()

        if verbose > 0:
            print('Adding weights', flush=True)

        self.add_weights(intx, option=weighting)
        justweight = intx[['weight']]
        overlays = {}

        for dimset, keys in dimsets.items():
            if verbose > 0:
                print('Working on', dimset, flush=True)
            outdims = griddims
            if verbose > 1:
                t0 = time.time()
                print(' - Making dataframe', flush=True)
            df = self.to_dataframe(*keys, geo=False)
            if verbose > 2:
                print(f'   - Dataframe has {df.shape[0]} rows', flush=True)
            if any([n in dimset for n in justweight.index.names]):
                notgeodims = set(dimset).difference(justweight.index.names)
                notgeodims = [dk for dk in dimset if dk in notgeodims]
                if verbose > 1:
                    print(' - Adding intx geometry', flush=True)
                    t0 = time.time()
                if len(notgeodims) > 0:
                    if verbose > 1:
                        print(f'   - Unstacking {notgeodims}', flush=True)
                    mjustweight = justweight.copy()
                    wkey = ('weight',) + ('',) * len(notgeodims)
                    wskey = ('weight_sum',) + ('',) * len(notgeodims)
                    mjustweight.columns = pd.MultiIndex.from_tuples(
                        [wkey], names=[None] + notgeodims
                    )
                    df = mjustweight.join(df.unstack(notgeodims))
                else:
                    wkey = 'weight'
                    # Adding on keyword, because it is significanly faster for
                    # some products. I would have thought the on would be
                    # inferred by names.
                    sharedidk = [
                        idk for idk in df.index.names
                        if idk in justweight.index.names
                    ]
                    df = justweight.join(df, on=sharedidk)

                if verbose > 2:
                    print(
                        f'   - Weighted inputs has {df.shape[0]} rows'
                        + f' ({time.time() - t0:.1f}s)',
                        flush=True
                    )
                    t0 = time.time()
                # Geometry objects are not weightable
                if verbose > 1:
                    print(' - Weighting', flush=True)
                gdf = grouped_weighted_avg(
                    df, df[wkey], outdims
                )
                if len(notgeodims) > 0:
                    ngdf = gdf.drop([wkey, wskey], axis=1).stack(notgeodims)
                    ngdf['weight'] = gdf[wkey]
                    ngdf['weight_sum'] = gdf[wskey]
                    gdf = ngdf
            else:
                gdf = self.to_dataframe(*keys, geo=False, valid=False)
                # tmpgrid = grid.drop('geometry', axis=1).copy()
                # tmpgrid['none'] = 1
                # gdf = tmpgrid.set_index('none', append=True).join(
                #     df.set_index('none', append=True)
                # ).droplevel('none')
                # gdf['weight_sum'] = 1

            if verbose > 2:
                print(
                    f'   - Weighted outputs has {gdf.shape[0]} rows'
                    f' ({time.time() - t0:.1f}s)',
                    flush=True
                )
            if verbose > 1:
                print(' - Adding attributes', flush=True)
            for key in gdf.columns:
                if key in self.ds.variables:
                    gdf[key].attrs.update(self.ds[key].attrs)

            overlays[tuple(dimset)] = gdf

        if as_dataset:
            import xarray as xr
            from datetime import datetime
            dss = {}
            for dimks, outdf in overlays.items():
                dss[dimks] = xr.Dataset.from_dataframe(outdf)
            outds = xr.merge(dss.values()).reindex(**{
                griddim: grid.index.get_level_values(griddim).unique()
                for griddim in griddims
            })
            for key in outds.variables:
                if key in self.ds.variables:
                    outds[key].attrs.update(self.ds[key].attrs)
            docstr = getattr(self, '__doc__', None)
            if docstr is None:
                docstr = ""
            outds.attrs['description'] = (
                docstr
                + '\n - {}'.format(getattr(self, 'path', '<path unknown>'))
            )
            outds.attrs['updated'] = datetime.now().strftime('%FT%H:%M:%S%z')
            outds.attrs['history'] = ''
            outds.attrs['crs'] = grid.crs.srs
            outds.attrs['cmaqsatproc_version'] = csp_version()
            comp = dict(zlib=True, complevel=1)
            for key, var in outds.data_vars.items():
                var.encoding.update(comp)

            return outds

        return overlays

    @classmethod
    def cmr_to_level3(
        cls, temporal, grid, griddims=None, weighting='area', bbox=None,
        verbose=0, varkeys=None, as_dataset=True, link_kwargs=None, **kwargs
    ):
        """
        Wrapper around cmr_links and paths_to_level3. For description of
        keywords, see those methods.
        """
        from copy import copy
        if link_kwargs is None:
            link_kwargs = {}
        link_kwargs = copy(link_kwargs)
        if bbox is None:
            grid.unary_union.envelope
        link_kwargs.setdefault('method', 'opendap')
        links = cls.cmr_links(
            temporal=temporal, bbox=bbox, **link_kwargs
        )
        if verbose > 1:
            print(len(links), links)
        elif verbose > 0:
            print(len(links), 'first link', links[0])
        return cls.paths_to_level3(
            links, grid=grid, griddims=griddims, weighting=weighting, bbox=bbox,
            verbose=verbose, varkeys=varkeys, as_dataset=as_dataset, **kwargs
        )

    @classmethod
    def paths_to_level3(
        cls, paths, grid, griddims=None, weighting='area', bbox=None,
        verbose=0, varkeys=None, as_dataset=True, **kwargs
    ):
        """
        Iteratively apply
        output = cls(path).to_level3(*args, path=path, **kwds)
        and then grouped_weighted_avg(output[dims]) for each dimset.

        For description of keywords, see to_level3.
        """
        import pandas as pd
        from copy import copy

        if varkeys is None:
            varkeys = ()

        withdata = {}
        nodata = {}
        if griddims is None:
            griddims = list(grid.index.names)

        for path in paths:
            if verbose > 0:
                print(path)
            try:
                sat = cls.open_dataset(path, bbox=bbox, **kwargs)
                withdata[path] = sat.to_level3(
                    *varkeys, grid=grid, griddims=griddims,
                    weighting=weighting, verbose=verbose, as_dataset=False
                )
            except Exception as e:
                nodata[path] = repr(e)
                if verbose > 0:
                    print(nodata[path])

        keys = sorted(withdata)
        dimdatasets = {}
        dimpaths = {}
        for path, output in withdata.items():
            for key, valdf in output.items():
                dimdatasets.setdefault(key, []).append(valdf)
                dimpaths.setdefault(key, []).append(path)
        outputs = {}
        for dimks, dimdfs in dimdatasets.items():
            keys = dimpaths[dimks]
            attrs = {
                key: copy(dimdfs[0][key].attrs) for key in dimdfs[0].columns
            }
            combinedf = pd.concat(
                dimdfs, keys=keys, names=['path']
            )
            outdims = list(dimdfs[0].index.names)
            if 'weight_sum' not in combinedf.columns:
                combinedf = combinedf.groupby(outdims).mean()
            else:
                combinedf = grouped_weighted_avg(
                    combinedf, combinedf['weight_sum'], outdims
                )

            for key in combinedf.columns:
                if key in attrs:
                    combinedf[key].attrs.update(attrs[key])
            outputs[dimks] = combinedf

        if as_dataset:
            import xarray as xr
            from datetime import datetime
            dss = {}
            for dimks, outdf in outputs.items():
                dss[dimks] = xr.Dataset.from_dataframe(outdf)

            outds = xr.merge(dss.values()).reindex(**{
                griddim: grid.index.get_level_values(griddim).unique()
                for griddim in griddims
            })
            for key in outds.variables:
                if key in sat.ds.variables:
                    outds[key].attrs.update(sat.ds[key].attrs)
            docstr = getattr(cls, '__doc__', None)
            if docstr is None:
                docstr = ""
            outds.attrs['description'] = (
                docstr
                + '\n - '.join(rootremover(paths, insert=True)[1])
            )
            outds.attrs['history'] = str(nodata)
            outds.attrs['updated'] = datetime.now().strftime('%FT%H:%M:%S%z')
            outds.attrs['crs'] = grid.crs.srs
            outds.attrs['cmaqsatproc_version'] = csp_version()
            return outds

        return outputs, nodata
