__all__ = ['satellite']

from ..utils import weight_vars, csp_formatwarnings


class satellite:
    __doc__ = """
`satellite` is the core satellite reader. It has the basic functionality
of converting a data store (i.e., HDF4, HDF5, or NetCDF file) to a dataframe
and georeferencing it.

A standard file is opened using the path to a data store. The path can be
any url that is readable by xarray. Commonly this will be local file paths
or OpenDAP urls to Level 2 or Level 3 satellite data.

The standard file then has methods to_dataframe and to_geodataframe. These
methods produce pandas or geopandas dataframes (like a csv file). The
geopandas has shapely Polygons or Points that identify the location of the
pixel.

A standard file then has general methods `weights` and `weighted`. `weights`
accepts a geopandas dataframe with polygons. Typically, the dataframe polygons
will be grid cells of a target domain. Weights calculates the contribution
weight for each pixel-gridcell pair. The typical weights algorithms include
equal (any overlapping pixel has an equal contribution) and fractional area
where the area overlap determines the contribution. `weighted` is a special
function that returns the weighted sum of pixels for each grid cell.

Conceptually, the weights method can be overwritten to include complex weights.
For example, the FOV brightness etc.
"""
    _varkeys2d = ()
    _varkeys3d = ()
    _zkey = None

    @property
    def short_description(self):
        import warnings
        desc = ""
        warnings.warn(
            f'short_description has not been implemented for this {type(self)}',
            RuntimeWarning
        )
        return desc

    @classmethod
    def from_paths(cls, paths, *keys, low_memory=True):
        """
        Open each data store (paths) and export to a dataframe all valid
        pixels for each path and then contruct a single object with access
        to all the pixels. This is particularly useful for multifile

        Arguments
        ---------
        paths : iterable
            Paths to individual data stores
        keys : list
            Variables to be exported. If none provided, standard export will
            be used.
        low_memory : bool
            Default True. Delete each object after it has been exported.
        """
        import gc
        dfs = []
        for path in paths:
            tmp = cls.open_path(path)
            df = tmp.export_dataframe(*keys)
            if low_memory:
                del tmp
                gc.collect()
            dfs.append(df)

        return cls.from_dataframes(dfs)

    @classmethod
    def open_path(cls, path, **kwds):
        """
        Initialize a satellite object where the datastore will be opened
        from the path. This is the same as teh typical __init__ call.
        """
        obj = cls(path, **kwds)
        return obj

    @classmethod
    def from_dataset(cls, ds, **kwds):
        """
        Initialize a satellite object and explicitly set the data store. This
        allows for a dataset to be created from other datasets and used as
        an input.
        """
        obj = cls(path='<inline>', **kwds)
        obj.ds = ds
        return obj

    @classmethod
    def from_dataframe(cls, df, **kwds):
        """
        Initialize a satellite object and explicitly set the dataframe
        that will be used to supersede a data store.

        Arguments
        ---------
        df : pandas.DataFrame
            Must have all columns needed by specific satelite reader.
        kwds : mappable
            Keyword arguments required for initialization
        """
        obj = cls(path='<inlinedf>', **kwds)
        obj.df = df
        return obj

    @classmethod
    def from_dataframes(cls, dfs, *keys):
        """
        Initialize a satellite object and explicitly set the dataframe
        that will be used to supersede a data store. The dataframe will
        be the concatenation of dataframes with attributes for each column
        and overall based on the attrs of the first df in dfs

        Arguments
        ---------
        dfs : list of pandas.DataFrame
            Each df must have all columns needed by specific satelite reader.
        kwds : mappable
            Keyword arguments required for initialization
        """
        import pandas as pd

        df = pd.concat(dfs, keys=range(len(dfs)), names=['df'])
        for key in df.columns:
            df[key].attrs.update(dfs[0][key].attrs)
        return cls.from_dataframe(df)

    def __init__(self, path=None, bbox=None, **attrs):
        """
        Arguments
        ---------
        path : str or Dataset
            If a str, xr.open_dataset to open the dataset as the ds property.
            If a Dataset, it will be used as the ds property.
        bbox : iterable
            Must contain swlon, swlat, nelon, nelat in decimal degrees. Will
            be used in valid_index to subset "valid" pixels to a domain. This
            is extremely good for speed becuase you often do not need all
            pixels of an orbit.
        attrs : mappable
            Used to set attrs of the
        """
        from copy import deepcopy

        self.path = path
        self._ds = None
        self._df = None
        self._bbox = bbox
        self._attrs = deepcopy(attrs)
        self._valididx = None
        self._geodf = None

    @property
    def required_args(self):
        return ()

    @property
    def valid_index(self):
        """
        This is a critical property that must be defined to filter good pixels.
        It should return a dataframe that has only valid rows and indices.

        Usually, this will be some sort of quality flag check
        """
        raise NotImplementedError('Not implemented')

    @property
    def geodf(self):
        """
        Create a base geopandas dataframe with a geometry for all valid_index
        elements.

        The geometry will typically be a polygon representing the pixel
        corners or a point representing the pixel centroid.
        """
        raise NotImplementedError('Not implemented')

    @property
    def ds(self):
        """
        The ds property contains a dataset, which must contain sufficient
        data for geodf to define a polygon. The ds may be superseded by
        a df if set.
        """
        import xarray as xr
        if self._ds is None:
            self._ds = xr.open_dataset(self.path)
        return self._ds

    @ds.setter
    def ds(self, ds):
        self._ds = ds

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, df):
        self._df = df

    @property
    def attrs(self):
        return self._attrs

    def to_dataframe(self, *keys, valid_only=True):
        """
        Export keys to a dataframe for only valid values.

        Arguments
        ---------
        keys : list
            List of keys to pull from the self.ds or self.df for all valid
            pixels.
        valid_only : bool
            Default True. Only export a dataframe where valid_index.

        Returns
        -------
        df : pandas.Dataframe
            Dataframe with columns keys
        """
        props = {}
        if self.df is None:
            df = self.ds[list(keys)].to_dataframe()
            for key in keys:
                if key in self.ds.data_vars:
                    props[key] = self.ds[key].attrs
        else:
            df = self.df[list(keys)]
            for key in keys:
                props[key] = self.df[key].attrs
        if valid_only:
            outdf = df.join(
                self.valid_index[[]], how='inner'
            )
        else:
            outdf = df

        for key in keys:
            if key in props:
                outdf[key].attrs.update(props[key])

        return outdf

    def to_geodataframe(self, *keys, geodf=None, valid_only=True):
        """
        Export keys to a dataframe for only valid values with geometries
        from geodf

        Arguments
        ---------
        keys : list
            List of keys to pull from the self.ds or self.df for all valid
            pixels.
        valid_only : bool
            Default True. Only export a dataframe where valid_index.

        Returns
        -------
        df : pandas.Dataframe
            Dataframe with columns keys with geometry
        """
        if geodf is None:
            geodf = self.geodf
        df = self.to_dataframe(*keys, valid_only=valid_only)
        outdf = geodf.join(df)
        for key in keys:
            outdf[key].attrs.update(df[key].attrs)

        return outdf

    def weights(self, othdf, option='fracarea', clip=None):
        """
        weights provides the weighting to apply to each L2 pixel. The default
        is the fractional area overlap of teh L2 pixel (tesselation) and the
        destination geometry.

        Arguments
        ---------
        othdf : geopandas.GeoDataFrame
            This must have a geometry column that defines the destination
            geometry.
        option : str
            choices fractional area overlap (fracarea) or equal (any overlap)

        Returns
        -------
        intx_df : geopandas.GeoDataFrame
            This is the same as the self.geodf.sjoin(othdf)
            where the fields are complemented by the destination geometry, the
            intersection area and the fractional intersection area. The final
            addition is weights, which defaults to the intx_fracarea.
        """
        import geopandas as gpd
        import warnings

        # Gets rid of self-intersections.
        mygeodf = self.geodf
        mygeodf = mygeodf[mygeodf.is_valid]
        if mygeodf.shape[0] < self.geodf.shape[0]:
            nlost = self.geodf.shape[0] - mygeodf.shape[0]
            warnings.warn(
                f'{nlost} geometries are not valid (cross pole or dateline)'
            )

        if clip is not None:
            mygeodf = gpd.clip(mygeodf, clip.to_crs(self.geodf.crs))
            if mygeodf.shape[0] == 0:
                intx_df = mygeodf.copy()
                if option == 'fracarea':
                    intx_df.insert(
                        column='intx_area', loc=len(intx_df.columns),
                        value=0
                    )
                    intx_df.insert(
                        column='intx_fracarea', loc=len(intx_df.columns),
                        value=0
                    )
                intx_df.insert(
                    column='weights', loc=len(intx_df.columns), value=0
                )
                return intx_df

        mygeodf_othcrs = mygeodf.to_crs(othdf.crs)

        # This is an uncertain assumption. Basically, I am having trouble with
        # pixels that cross the dateline, but do not encompass the dateline.
        # For example, a polygon that crosses the dateline may be interpretted
        # to include the entire world excluding the dateline.
        # The polygon shown below is 358 degrees wide
        # POLYGON ((179 0, -179 0, -179 2, 179 2)))
        # q1, q3 = area_othcrs.quantile([.25, .75])
        # iqr = q3 - q1
        # ub = q3 + 2 * iqr
        # area_othcrs = mygeodf_othcrs.area
        # isnotoutlier = area_othcrs <= ub
        # Exclude very large pixels
        # outlierfrac = 1 - isnotoutlier.mean()
        # if outlierfrac > 0.1:
        #     warnings.warn(
        #         f'{outlierfrac:%} were removed in an attempt to limit '
        #         + 'cross-dateline issues'
        #     )
        # mygeodf_othcrs_valid = mygeodf_othcrs[isnotoutlier]
        # Testing solution using _dateline
        mygeodf_othcrs_valid = mygeodf_othcrs
        intx_df = gpd.sjoin(mygeodf_othcrs_valid, othdf.reset_index())
        intx_df['dest_geometry'] = othdf.reset_index().iloc[
            intx_df.index_right
        ].geometry.values
        if option == 'fracarea':
            intx_df['intx_area'] = intx_df.geometry.intersection(
                intx_df['dest_geometry']
            ).area
            intx_df['intx_fracarea'] = (
                intx_df['intx_area'] / intx_df['dest_geometry'].area
            )
            intx_df['weights'] = intx_df['intx_fracarea']
        elif option == 'equal':
            intx_df['weights'] = 1
        else:
            import warnings
            warnings.warn('Using default weight of 1')
            intx_df['weights'] = 1
        # print(intx_df.shape)
        # import pdb; pdb.set_trace()
        return intx_df

    def weighted(
        self, *aggkeys, groupkeys, wgtdf=None, othdf=None, **weights_kw
    ):
        """
        Apply weights to all numeric fields and calculate output as:
        out = val * weights / sum(weights)

        weights will be returned as sum(weights)

        Arguments
        ---------
        aggkeys : list
            List of keys to aggregate at the cell othdf or wgtdf object.
        groupkeys : list
            List of keys to group by (e.g., ROW, COL)
        wgtdf : geopandas.GeoDataFrame
            This must have a weight column that defines the contribution
            to the groupby keys combination.
        othdf : geopandas.GeoDataFrame
            This must have a geometry column that defines the destination
            geometry. It will be used to calculate a wgtdf using the default
            weights options.
        weight_kw : mappable
            Passed through to to weights

        Returns
        -------
        intx_df : geopandas.GeoDataFrame
            This is the same as the self.geodf.sjoin(othdf)
            where the fields are complemented by the destination geometry, the
            intersection area and the fractional intersection area. The final
            addition is weights, which defaults to the intx_fracarea.
        """
        if wgtdf is None:
            wgtdf = self.weights(othdf, **weights_kw)
        if wgtdf.shape[0] == 0:
            aggdf = wgtdf.copy()
            for key in groupkeys:
                if key not in aggdf.columns:
                    aggdf.insert(column=key, loc=0, value=0)
            for key in aggkeys:
                if key not in aggdf.columns:
                    aggdf.insert(column=key, loc=0, value=0)
            return aggdf.set_index(groupkeys)

        withweights = self.to_geodataframe(*aggkeys, geodf=wgtdf)
        return weight_vars(withweights, groupkeys=groupkeys)

    def export_dataframe(self, *args):
        """
        Similar to to_dataframe, except it will be exported with sufficient
        information to construct polygons and filter data. This was designed
        primarily for creating merged objects (see from_dataframes and
        from_paths)
        """
        allargs = list(set(tuple(self.required_args) + tuple(args)))
        outdf = self.to_dataframe(*allargs)
        for key in allargs:
            outdf[key].attrs.update(outdf[key].attrs)

        return outdf

    def to_level3(
        self, grid, varkeys2d=None, varkeys3d=None, groupkeys=None,
        weights_kw=None, verbose=0
    ):
        """
        Convert level2 data to level3 data.

        For 2d and 3d variables:
            1. For each link, calculate pixel contributions to all grid cells.
              * This process is specific to each satellite reader.
              * requires a reader object with the weights method
            2. For each link, weight pixel values and produces grid cell values.
              * This process is specific to each satellite reader.
              * requires a reader object with the weighted method
            3. Combine weighted pixel values to produce overall value

        Arguments
        ---------
        grid : CMAQGrid
            Grid object.
        varkeys2d : iterable
            List of keys to process with first 2 dimensions in groupkeys
        varkeys3d : iterable
            List of keys to process with all dimensions in groupkeys
        groupkeys : list
            Usually ['ROW', 'COL', zdim] where zdim will depend on the dataset
        weights_kw : mappable
            Keywords for the weights method. For example, dict(option='equal')
            to make all intersecting pixels equally weighted.
        verbose : int
            Level of verbosity

        Returns
        -------
        out : dict
            '2d': 2d variables as a dataframe,
            '3d': 3d variables as a dataframe,
            'description': a short description of the processing,
            'history': a list of all warnings
        """
        import warnings
        if self.df is not None:
            allkeys = list(self.df.columns)
        else:
            allkeys = list(self.ds.data_vars)

        if varkeys2d is None:
            varkeys2d = [
                key for key in self._varkeys2d if key in allkeys
            ]
        if varkeys3d is None:
            varkeys3d = [
                key for key in self._varkeys3d if key in allkeys
            ]
        if groupkeys is None:
            if self._zkey is None:
                groupkeys = ['ROW', 'COL', ]
            else:
                groupkeys = ['ROW', 'COL', self._zkey]
        if weights_kw is None:
            weights_kw = {}
        with warnings.catch_warnings(record=True) as warning_list:
            outputs = {}
            outputs['description'] = (
                f'{self.path} with:\n{self.short_description}'
            )
            if verbose > 1:
                print(f'Found {self.valid_index.shape[0]} valid pixels')

            if verbose > 0:
                print('Calculating weights', flush=True)

            wgts = self.weights(
                grid.geodf, clip=grid.exterior, **weights_kw
            )
            if wgts.shape[0] == 0:
                warnings.warn(
                    f'No valid pixels for {self.path}'
                )
                outputs['history'] = csp_formatwarnings(warning_list)
                return outputs

            if len(varkeys2d) > 0:
                if verbose > 0:
                    print('Processing 2d data', flush=True)

                # Produce 2d output
                wgtd2d = self.weighted(
                    *varkeys2d, groupkeys=groupkeys[:2], wgtdf=wgts
                )
                outputs['2d'] = wgtd2d
                if verbose > 1:
                    print(f'{wgtd2d.shape[0]} unique cells')

            if len(varkeys3d) > 0:
                if verbose > 0:
                    print('Processing 3d data', flush=True)

                # Produce 3d output
                wgtd3d = self.weighted(
                    *varkeys3d, groupkeys=groupkeys, wgtdf=wgts
                )
                outputs['3d'] = wgtd3d
                if verbose > 1:
                    print(f'{wgtd3d.shape[0]} unique cells-layers')
            outputs['history'] = csp_formatwarnings(warning_list)
        return outputs
