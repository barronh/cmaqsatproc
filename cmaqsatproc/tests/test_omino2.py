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
        'nTimes', 'nXtrack', 'LL_Longitude', 'LU_Longitude', 'UL_Longitude',
        'UU_Longitude', 'LL_Latitude', 'LU_Latitude', 'UL_Latitude',
        'UU_Latitude', 'CloudFraction', 'VcdQualityFlags',
        'XTrackQualityFlags', 'ColumnAmountNO2Std'
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
    df['Latitude'] = (
        df['LL_Latitude'] + df['LU_Latitude']
        + df['UL_Latitude'] + df['UU_Latitude']
    ) / 4
    df['Longitude'] = (
        df['LL_Longitude'] + df['LU_Longitude']
        + df['UL_Longitude'] + df['UU_Longitude']
    ) / 4
    return df


def omi_example_ds():
    import xarray as xr

    ds = xr.Dataset.from_dataframe(omi_example_df())
    ds['FoV75CornerLongitude'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'LL_Longitude', 'UL_Longitude', 'UU_Longitude', 'LU_Longitude'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    ds['FoV75CornerLatitude'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'LL_Latitude', 'UL_Latitude', 'UU_Latitude', 'LU_Latitude'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    outds = ds.drop_vars([
        'LL_Latitude', 'UL_Latitude', 'UU_Latitude', 'LU_Latitude',
        'LL_Longitude', 'UL_Longitude', 'UU_Longitude', 'LU_Longitude'
    ])
    return outds


def test_omi():
    df = omi_example_df()
    ds = omi_example_ds()

    sat = readers.omi.OMNO2.from_dataset(ds)
    df2 = sat.export_dataframe()
    assert((df2 == df.loc[:, df2.columns]).all().all())


def test_omi_valid():
    df = omi_example_df()
    ds = omi_example_ds()
    ds['VcdQualityFlags'][1, 1] = 1
    ds['CloudFraction'][1, 2] = 1
    ds['XTrackQualityFlags'][2, 0] = 1
    sat = readers.omi.OMNO2.from_dataset(ds)
    df2 = sat.export_dataframe()
    assert((df.shape[0] - 3) == df2.shape[0])


def test_satellite_weights():
    from .. import cmaq

    ds = omi_example_ds()
    sat = readers.omi.OMNO2.from_dataset(ds)

    cg = cmaq.CMAQGrid(None, '108US1')
    wgt = sat.weights(cg.geodf)
    wgtd1 = sat.weighted(
        'ColumnAmountNO2Std', wgtdf=wgt, groupkeys=['ROW', 'COL']
    )
    assert(
        wgtd1['ColumnAmountNO2Std'].min() >= ds['ColumnAmountNO2Std'].min()
    )
    assert(
        wgtd1['ColumnAmountNO2Std'].max() <= ds['ColumnAmountNO2Std'].max()
    )
    wgtd2 = sat.weighted(
        'ColumnAmountNO2Std', othdf=cg.geodf, groupkeys=['ROW', 'COL']
    )
    assert((wgtd1 == wgtd2).all().all())
