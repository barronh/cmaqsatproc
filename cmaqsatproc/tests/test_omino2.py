from .. import readers


def omi_example_df():
    """
    Example Derived from
    https://aura.gesdisc.eosdis.nasa.gov/opendap/Aura_OMI_Level2/OMNO2.003/2018
    /152/OMI-Aura_L2-OMNO2_2018m0601t1550-o73823_v003-2019m0819t210955.he5
    I exported all important fields and selected 3 nTimes and 3 nXtrack
    Then, I add mask values
    """
    import io
    import pandas as pd

    names = [
        'nTimes', 'nXtrack', 'll_x', 'lu_x', 'ul_x', 'uu_x', 'll_y', 'lu_y',
        'ul_y', 'uu_y', 'CloudFraction', 'VcdQualityFlags',
        'XTrackQualityFlags', 'ColumnAmountNO2'
    ]
    df = pd.read_csv(
        io.StringIO("""
990,15,-59.44,-59.5,-59.04,-59.1,46.73,46.86,46.83,46.96,0.04,0.0,0.0,811e12
990,16,-59.04,-59.1,-58.65,-58.72,46.83,46.96,46.92,47.05,0.03,0.0,0.0,554e12
990,17,-58.65,-58.72,-58.28,-58.34,46.92,47.05,47.0,47.13,0.02,0.0,0.0,678e12
991,15,-59.49,-59.56,-59.09,-59.16,46.85,46.98,46.94,47.08,0.05,0.0,0.0,160e12
991,16,-59.09,-59.16,-58.71,-58.77,46.94,47.08,47.03,47.16,0.07,0.0,0.0,965e12
991,17,-58.71,-58.77,-58.34,-58.4,47.03,47.16,47.12,47.25,0.03,0.0,0.0,877e12
992,15,-59.55,-59.62,-59.15,-59.21,46.96,47.1,47.06,47.19,0.08,0.0,0.0,111e12
992,16,-59.15,-59.21,-58.77,-58.83,47.06,47.19,47.15,47.28,0.08,0.0,0.0,155e12
992,17,-58.77,-58.83,-58.39,-58.46,47.15,47.28,47.24,47.36,0.04,0.0,0.0,555e12
"""),
        names=names
    ).set_index(['nTimes', 'nXtrack'])
    df['cn_y'] = df['Latitude'] = (
        df['ll_y'] + df['lu_y'] + df['ul_y'] + df['uu_y']
    ) / 4
    df['cn_x'] = df['Longitude'] = (
        df['ll_x'] + df['lu_x'] + df['ul_x'] + df['uu_x']
    ) / 4
    return df


def omi_example_ds():
    import xarray as xr

    ds = xr.Dataset.from_dataframe(omi_example_df())
    ds['FoV75CornerLongitude'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'll_x', 'ul_x', 'uu_x', 'lu_x'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    ds['FoV75CornerLatitude'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'll_y', 'ul_y', 'uu_y', 'lu_y'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    ds['valid'] = (
        ((ds['VcdQualityFlags'].astype('i') & 1) == 0)
        & (ds['XTrackQualityFlags'] == 0)
        & (ds['CloudFraction'] <= 0.3)
    )
    # outds = ds.drop_vars([
    #     'll_y', 'ul_y', 'uu_y', 'lu_y', 'll_x', 'ul_x', 'uu_x', 'lu_x'
    # ])
    outds = ds
    return outds


def checksat(sat):
    import geopandas as gpd
    from shapely.geometry import box
    df = omi_example_df()
    satdf = sat.to_dataframe('ColumnAmountNO2')
    assert ((df['ColumnAmountNO2'] == satdf['ColumnAmountNO2']).all())
    gdf = gpd.GeoDataFrame(
        dict(ROW=[0], COL=[0], Val=[df['ColumnAmountNO2'].mean()]),
        geometry=[box(-179, 0, 179, 89)], crs=4326
    ).set_index(['ROW', 'COL']).to_crs(3785)
    l3 = sat.to_level3('ColumnAmountNO2', grid=gdf[['geometry']])
    diff = (
        l3['ColumnAmountNO2'].values.ravel() - gdf['Val'].values
    )
    pctdiff = diff / gdf['Val'] * 100
    assert ((pctdiff.abs() < 5).all())


def test_from_dataset():
    ds = omi_example_ds()
    sat = readers.omi.OMNO2.from_dataset(ds)
    checksat(sat)


def test_open_dataset():
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdirname:
        outpath = os.path.join(tmpdirname, 'testomino2.nc')
        ds = omi_example_ds()
        ds.to_netcdf(outpath)
        sat = readers.omi.OMNO2.open_dataset(outpath)
        checksat(sat)
