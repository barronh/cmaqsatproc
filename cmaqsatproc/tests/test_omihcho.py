from .. import readers


def omi_example_df():
    """
    Example Derived from
    https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMHCHO.003/
    2018/152/OMI-Aura_L2-OMHCHO_2018m0601t1550-o73823_v003-2018m0601t230054.he5
    I exported all important fields and selected 3 nTimes and 3 nXtrack
    Then, I add mask values
    """
    import io
    import pandas as pd
    names = [
        'nTimes', 'nXtrack', 'll_x', 'lu_x', 'ul_x', 'uu_x', 'll_y', 'lu_y',
        'ul_y', 'uu_y', 'MainDataQualityFlag', 'ColumnUncertainty',
        'XtrackQualityFlagsExpanded', 'AMFCloudFraction', 'AirMassFactor',
        'ReferenceSectorCorrectedVerticalColumn'
    ]
    df = pd.read_csv(
        io.StringIO("""
990,15,-59.44,-59.04,-59.5,-59.1,46.76,46.85,46.87,46.97,0.0,89e14,0,0,1.1,79e14
990,16,-59.04,-58.66,-59.1,-58.72,46.85,46.94,46.97,47.06,0.0,104e14,0,0,1.14,55e14
990,17,-58.66,-58.29,-58.72,-58.35,46.94,47.03,47.06,47.14,0.0,61e14,0,0,1.24,13e15
991,15,-59.5,-59.1,-59.55,-59.16,46.87,46.97,46.99,47.09,0.0,151e14,0,0,1.09,90e14
991,16,-59.1,-58.72,-59.16,-58.77,46.97,47.06,47.09,47.18,0.0,124e14,0,0,1.12,0e14
991,17,-58.72,-58.35,-58.77,-58.4,47.06,47.14,47.18,47.26,0.0,58e14,0,0,1.18,95e14
992,15,-59.55,-59.16,-59.61,-59.21,46.99,47.09,47.11,47.2,0.0,77e14,0,0,1.05,166e14
992,16,-59.16,-58.77,-59.21,-58.83,47.09,47.18,47.2,47.29,0.0,10e15,0,0,1.02,-18e14
992,17,-58.77,-58.4,-58.83,-58.46,47.18,47.26,47.29,47.38,0.0,68e14,0,0,1.16,-78e14
"""),
        names=names
    ).set_index(['nTimes', 'nXtrack'])
    df['cn_y'] = df['Latitude'] = (
        df['ll_y'] + df['lu_y']
        + df['ul_y'] + df['uu_y']
    ) / 4
    df['cn_x'] = df['Longitude'] = (
        df['ll_x'] + df['lu_x']
        + df['ul_x'] + df['uu_x']
    ) / 4
    return df


def omi_example_ds():
    import xarray as xr
    import numpy as np

    ds = xr.Dataset.from_dataframe(omi_example_df())
    # Need to reconfigure to have nTimes_1 and nXtrack_1 which are
    # edge dimensions
    ds['PixelCornerLongitudes'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'll_x', 'ul_x', 'uu_x', 'lu_x'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    # Need to reconfigure to have nTimes_1 and nXtrack_1 which are
    # edge dimensions
    ds['PixelCornerLatitudes'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'll_y', 'ul_y', 'uu_y', 'lu_y'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    ds['PixelCornerLatitudes'] = lat_edges = xr.DataArray(
        np.zeros((ds.dims['nTimes'] + 1, ds.dims['nXtrack'] + 1)),
        dims=('nTimes_1', 'nXtrack_1')
    )
    lat_edges[:-1, :-1] = ds['ll_y'].values
    lat_edges[:-1, -1] = ds['lu_y'][:, -1].values
    lat_edges[-1, :-1] = ds['ul_y'][-1].values
    lat_edges[-1, -1] = ds['uu_y'][-1, -1].values
    ds['PixelCornerLongitudes'] = lon_edges = xr.DataArray(
        np.zeros((ds.dims['nTimes'] + 1, ds.dims['nXtrack'] + 1)),
        dims=('nTimes_1', 'nXtrack_1')
    )
    lon_edges[:-1, :-1] = ds['ll_x'].values
    lon_edges[:-1, -1] = ds['lu_x'][:, -1].values
    lon_edges[-1, :-1] = ds['ul_x'][-1].values
    lon_edges[-1, -1] = ds['uu_x'][-1, -1].values

    ds['valid'] = (
        (ds['MainDataQualityFlag'] == 0)
        & ((ds['XtrackQualityFlagsExpanded'].astype('i') & 1) == 0)
        & (ds['AMFCloudFraction'] <= 0.3)
    )
    outds = ds
    # .drop_vars([
    #     'll_y', 'ul_y', 'uu_y', 'lu_y',
    #     'll_x', 'ul_x', 'uu_x', 'lu_x'
    # ])
    return outds


def checksat(sat):
    import geopandas as gpd
    from shapely.geometry import box

    df = omi_example_df()

    inkey = 'ReferenceSectorCorrectedVerticalColumn'
    satdf = sat.to_dataframe(inkey)
    assert ((df[inkey] == satdf[inkey]).all())
    gdf = gpd.GeoDataFrame(
        dict(ROW=[0], COL=[0], Val=[df[inkey].mean()]),
        geometry=[box(-179, 0, 179, 89)], crs=4326
    ).set_index(['ROW', 'COL']).to_crs(3785)
    l3 = sat.to_level3(inkey, grid=gdf[['geometry']])
    diff = (
        l3[inkey].values.ravel() - gdf['Val'].values
    )
    pctdiff = diff / gdf['Val'] * 100
    assert ((pctdiff.abs() < 5).all())


def test_from_dataset():
    ds = omi_example_ds()
    sat = readers.omi.OMHCHO.from_dataset(ds)
    checksat(sat)


def test_open_dataset():
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdirname:
        outpath = os.path.join(tmpdirname, 'testomihcho.nc')
        ds = omi_example_ds()
        ds.to_netcdf(outpath)
        sat = readers.omi.OMHCHO.open_dataset(outpath)
        checksat(sat)
