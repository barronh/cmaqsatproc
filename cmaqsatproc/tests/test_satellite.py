from .. import readers


def dataset_example():
    import xarray as xr
    ds = xr.Dataset(
        data_vars=dict(
            Flag=xr.DataArray(
                [[0, 1, 0], [0, 1, 0]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='none')
            ),
            Val=xr.DataArray(
                [[1, -999, 2], [3, -999, 4]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='DU')
            ),
            LL_Longitude=xr.DataArray(
                [[-90, -91, -92], [-90, -91, -92]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            LL_Latitude=xr.DataArray(
                [[40, 40, 40], [41, 41, 41]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            UL_Longitude=xr.DataArray(
                [[-90, -91, -92], [-90, -91, -92]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            UL_Latitude=xr.DataArray(
                [[41, 41, 41], [42, 42, 42]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            UU_Longitude=xr.DataArray(
                [[-91, -92, -93], [-91, -92, -93]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            UU_Latitude=xr.DataArray(
                [[41, 41, 41], [42, 42, 42]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            LU_Longitude=xr.DataArray(
                [[-91, -92, -93], [-91, -92, -93]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            LU_Latitude=xr.DataArray(
                [[40, 40, 40], [41, 41, 41]], dims=('nTime', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
        )
    )
    return ds


class satellite_example(readers.satellite):
    """
    Silly example designed for testing with the dataset from the
    dataset_example function.
    """
    @property
    def valid_index(self):
        if self._valididx is None:
            df = self.to_dataframe('Flag', valid_only=False)
            self._valididx = df.query('Flag == 0')
        return self._valididx

    @property
    def geodf(self):
        if self._geodf is None:
            import geopandas as gpd
            from shapely.geometry import Polygon
            df = self.to_dataframe(
                'LL_Latitude', 'LL_Longitude',
                'LU_Latitude', 'LU_Longitude',
                'UL_Latitude', 'UL_Longitude',
                'UU_Latitude', 'UU_Longitude',
            ).join(self.valid_index[[]], how='inner')
            self._geodf = gpd.GeoDataFrame(
                df[[]],
                geometry=df.apply(
                    lambda row: Polygon([
                        [row['LL_Longitude'], row['LL_Latitude']],
                        [row['LU_Longitude'], row['LU_Latitude']],
                        [row['UU_Longitude'], row['UU_Latitude']],
                        [row['UL_Longitude'], row['UL_Latitude']],
                        [row['LL_Longitude'], row['LL_Latitude']],
                    ]), axis=1
                ),
                crs=4326
            )

        return self._geodf


def test_satellite():
    import numpy as np

    ds = dataset_example()
    sat = satellite_example.from_dataset(ds)
    sat2 = satellite_example.from_dataframe(ds.to_dataframe())
    geodf = sat.geodf
    geodf2 = sat2.geodf
    assert(geodf.shape[0] == 4)
    assert((geodf == geodf2).all().all())
    assert(np.allclose(sat.to_dataframe('Val')['Val'].values, [1, 2, 3, 4]))


def test_satellite_weights():
    from .. import cmaq

    ds = dataset_example()
    sat = satellite_example.from_dataset(ds)

    cg = cmaq.CMAQGrid(None, '108US1')
    wgt = sat.weights(cg.geodf)
    wgtd1 = sat.weighted('Val', wgtdf=wgt, groupkeys=['ROW', 'COL'])
    assert(wgtd1['Val'].min() >= ds['Val'].min())
    assert(wgtd1['Val'].max() <= ds['Val'].max())
    wgtd2 = sat.weighted('Val', othdf=cg.geodf, groupkeys=['ROW', 'COL'])
    assert((wgtd1 == wgtd2).all().all())
