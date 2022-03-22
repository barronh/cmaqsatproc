__all__ = ['satellite']


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

    def __init__(self, path=None, **kwds):
        from copy import deepcopy

        self.path = path
        self._ds = None
        self._df = None
        self._attrs = deepcopy(kwds)
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
            Default True. Only export a dataframe where valid_idx.

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

    def to_geodataframe(self, *keys, geodf=None):
        """
        Export keys to a dataframe for only valid values with geometries
        from geodf

        Arguments
        ---------
        keys : list
            List of keys to pull from the self.ds or self.df for all valid
            pixels.
        valid_only : bool
            Default True. Only export a dataframe where valid_idx.

        Returns
        -------
        df : pandas.Dataframe
            Dataframe with columns keys with geometry
        """
        if geodf is None:
            geodf = self.geodf
        df = self.to_dataframe(*keys)
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
        if clip is None:
            mygeodf = self.geodf
        else:
            mygeodf = gpd.clip(self.geodf, clip.to_crs(self.geodf.crs))
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

        intx_df = gpd.sjoin(mygeodf.to_crs(othdf.crs), othdf.reset_index())
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
                aggdf.insert(column=key, loc=0, value=0)
            for key in aggkeys:
                aggdf.insert(column=key, loc=0, value=0)
            return aggdf.set_index(groupkeys)

        withweights = self.to_geodataframe(*aggkeys, geodf=wgtdf)
        aggattrs = {
            key: withweights[key].attrs for key in aggkeys
        }
        notweighted = ['weights'] + list(groupkeys)
        weighted = withweights.drop(
            ['geometry', 'dest_geometry'], axis=1, inplace=False
        ).multiply(withweights.weights, axis=0)
        # overwrite weights with original values
        for key in notweighted:
            if key not in weighted.index.names:
                weighted[key] = withweights[key]
        # Groupby and sum
        aggweighted = weighted.groupby(groupkeys).sum()
        # Divide weighted values by the sum of weights
        agg = aggweighted.divide(aggweighted.weights, axis=0)
        # Overwrite with the sum of weights
        for key in notweighted:
            if key not in groupkeys:
                agg[key] = aggweighted[key]

        for key in aggkeys:
            agg[key].attrs.update(aggattrs[key])
        return agg

    def export_dataframe(self, *args):
        """
        Similar to to_dataframe, except it will be exported with sufficient
        information to construct polygons and filter data. This was designed
        primarily for creating merged objects (see from_dataframes and
        from_paths)
        """
        allargs = tuple(self.required_args) + tuple(args)
        outdf = self.to_dataframe(*allargs)
        for key in allargs:
            outdf[key].attrs.update(outdf[key].attrs)

        return outdf
