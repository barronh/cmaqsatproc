from .. import readers


def dataset_example():
    import xarray as xr
    ds = xr.Dataset(
        data_vars=dict(
            valid=xr.DataArray(
                [
                    [True, False, True],
                    [True, False, True]
                ], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='none')
            ),
            Val=xr.DataArray(
                [[1, -999, 2], [3, -999, 4]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='DU')
            ),
            ll_x=xr.DataArray(
                [[-90, -91, -92], [-90, -91, -92]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            ll_y=xr.DataArray(
                [[40, 40, 40], [41, 41, 41]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            lu_x=xr.DataArray(
                [[-90, -91, -92], [-90, -91, -92]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            lu_y=xr.DataArray(
                [[41, 41, 41], [42, 42, 42]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            uu_x=xr.DataArray(
                [[-91, -92, -93], [-91, -92, -93]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            uu_y=xr.DataArray(
                [[41, 41, 41], [42, 42, 42]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
            ul_x=xr.DataArray(
                [[-91, -92, -93], [-91, -92, -93]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_east')
            ),
            ul_y=xr.DataArray(
                [[40, 40, 40], [41, 41, 41]], dims=('nTimes', 'nXtrack'),
                attrs=dict(units='degrees_north')
            ),
        )
    )
    ds['cn_x'] = (
        ds['ll_x'] + ds['ul_x']
        + ds['lu_x'] + ds['uu_x']
    ) / 4
    ds['cn_y'] = (
        ds['ll_y'] + ds['ul_y']
        + ds['lu_y'] + ds['uu_y']
    ) / 4
    return ds


def test_satellite():
    import numpy as np
    import geopandas as gpd
    import warnings
    from shapely.geometry import box
    ds = dataset_example()
    sat = readers.satellite.from_dataset(ds)
    df = sat.to_dataframe('Val')
    assert (df.shape[0] == 4)
    assert (np.allclose(df['Val'].values, [1, 2, 3, 4]))
    df = sat.to_dataframe(geo=True)
    with warnings.catch_warnings(record=True):
        areas = df['geometry'].area
    assert ((areas == 1).all())
    gdf = gpd.GeoDataFrame(
        dict(ROW=[0], COL=[0], Val=[10 / 4.]),
        geometry=[box(-90, 40, -93, 42)], crs=4326
    ).set_index(['ROW', 'COL'])
    with warnings.catch_warnings(record=True):
        l3 = sat.to_level3('Val', grid=gdf[['geometry']])
    assert ((l3['Val'].values.ravel() == gdf['Val'].values).all())
