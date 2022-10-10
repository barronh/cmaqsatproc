from .. import readers


def modis_example_df():
    # Example Derived from
    # https://ladsweb.modaps.eosdis.nasa.gov/opendap/allData/61/MOD04_L2/2003/
    # 004/MOD04_L2.A2003004.1515.061.2017257143814.hdf
    # I exported all important fields and selected 3 Cell_Across_Swath and
    # 3 Cell_Along_Swath
    # Then, I add mask values
    import io
    import pandas as pd
    names = [
        'Cell_Across_Swath', 'Cell_Along_Swath', 'Longitude', 'Latitude',
        'dx', 'dy', 'Land_Ocean_Quality_Flag', 'Optical_Depth_Land_And_Ocean'
    ]
    df = pd.read_csv(
        io.StringIO("""
60,200,-69.599884,22.378624,0.097099304,0.089149475,3.0,0.134
60,201,-69.620316,22.289364,0.097042084,0.08916378,3.0,0.15
60,202,-69.641106,22.200296,0.096969604,0.08906746,3.0,0.157
61,200,-69.50296,22.364727,0.096767426,0.08914471,3.0,0.139
61,201,-69.523445,22.275476,0.09671402,0.08915424,3.0,0.139
61,202,-69.54431,22.186419,0.096637726,0.08905792,3.0,0.134
62,200,-69.40635,22.350815,0.09648132,0.08913803,3.0,0.126
62,201,-69.42689,22.26157,0.09642792,0.089146614,3.0,0.190
62,202,-69.44783,22.172522,0.09635162,0.089048386,3.0,0.130
"""),
        names=names
    ).set_index(['Cell_Across_Swath', 'Cell_Along_Swath'])
    return df


def modis_example_ds():
    import xarray as xr
    import numpy as np

    ds = xr.Dataset.from_dataframe(modis_example_df())
    Cell_Along_Swath_Edges = np.concatenate(
        [
            ds.Cell_Along_Swath[:1],
            (
                ds.Cell_Along_Swath[1:].values
                + ds.Cell_Along_Swath[:-1].values
            ) / 2,
            ds.Cell_Along_Swath[-1:]
        ], axis=0
    )
    Cell_Across_Swath_Edges = np.concatenate(
        [
            ds.Cell_Across_Swath[:1],
            (
                ds.Cell_Across_Swath[1:].values
                + ds.Cell_Across_Swath[:-1].values
            ) / 2,
            ds.Cell_Across_Swath[-1:]
        ], axis=0
    )
    lat_edges = ds.Latitude.interp(
        Cell_Along_Swath=Cell_Along_Swath_Edges,
        Cell_Across_Swath=Cell_Across_Swath_Edges
    )
    lon_edges = ds.Longitude.interp(
        Cell_Along_Swath=Cell_Along_Swath_Edges,
        Cell_Across_Swath=Cell_Across_Swath_Edges
    )
    corner_slices = {
        'll': (slice(None, -1), slice(None, -1)),
        'lu': (slice(None, -1), slice(1, None)),
        'ul': (slice(1, None), slice(None, -1)),
        'uu': (slice(1, None), slice(1, None)),
    }
    for corner, (xslice, aslice) in corner_slices.items():
        ds[f'{corner}_x'] = xr.DataArray(
            lon_edges.isel(
                Cell_Along_Swath=aslice, Cell_Across_Swath=xslice
            ),
            dims=ds.Longitude.dims,
            coords=ds.Longitude.coords,
        )
        ds[f'{corner}_y'] = xr.DataArray(
            lat_edges.isel(
                Cell_Along_Swath=aslice, Cell_Across_Swath=xslice
            ),
            dims=ds.Longitude.dims,
            coords=ds.Longitude.coords,
        )
    ds['cn_y'] = (
        ds['ll_y'] + ds['lu_y']
        + ds['ul_y'] + ds['uu_y']
    ) / 4
    ds['cn_x'] = (
        ds['ll_x'] + ds['lu_x']
        + ds['ul_x'] + ds['uu_x']
    ) / 4
    ds['valid'] = ds['Land_Ocean_Quality_Flag'] > 1
    outds = ds.drop_vars(['dx', 'dy'])
    return outds


def checksat(sat):
    import geopandas as gpd
    from shapely.geometry import box

    df = modis_example_df()
    inkey = 'Optical_Depth_Land_And_Ocean'
    satdf = sat.to_dataframe(inkey)
    assert ((df[inkey] == satdf[inkey]).all())
    gdf = gpd.GeoDataFrame(
        dict(ROW=[0], COL=[0], Val=[df[inkey].mean()]),
        geometry=[box(-179, 0, 179, 89)], crs=4326
    ).set_index(['ROW', 'COL']).to_crs(3785)
    l3 = sat.to_level3(inkey, grid=gdf[['geometry']])
    diff = (
        l3[inkey].values[0] - gdf['Val'].values[0]
    )
    pctdiff = diff / gdf['Val'] * 100
    assert ((pctdiff.abs() < 5).all())


def test_from_dataset():
    ds = modis_example_ds()
    sat = readers.modis.MOD04.from_dataset(ds)
    checksat(sat)


def test_open_dataset():
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdirname:
        outpath = os.path.join(tmpdirname, 'testomihcho.nc')
        ds = modis_example_ds()
        ds.to_netcdf(outpath)
        sat = readers.modis.MOD04.open_dataset(outpath)
        checksat(sat)
