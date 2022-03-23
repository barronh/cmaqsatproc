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
        'nTimes', 'nXtrack', 'LL_Longitude', 'LU_Longitude', 'UL_Longitude',
        'UU_Longitude', 'LL_Latitude', 'LU_Latitude', 'UL_Latitude',
        'UU_Latitude', 'MainDataQualityFlag', 'ColumnUncertainty',
        'AirMassFactor', 'ReferenceSectorCorrectedVerticalColumn'
    ]
    df = pd.read_csv(
        io.StringIO("""
990,15,-59.44,-59.04,-59.5,-59.1,46.76,46.85,46.87,46.97,0.0,89e14,1.1,79e14
990,16,-59.04,-58.66,-59.1,-58.72,46.85,46.94,46.97,47.06,0.0,104e14,1.14,55e14
990,17,-58.66,-58.29,-58.72,-58.35,46.94,47.03,47.06,47.14,0.0,61e14,1.24,13e15
991,15,-59.5,-59.1,-59.55,-59.16,46.87,46.97,46.99,47.09,0.0,151e14,1.09,90e14
991,16,-59.1,-58.72,-59.16,-58.77,46.97,47.06,47.09,47.18,0.0,124e14,1.12,0e14
991,17,-58.72,-58.35,-58.77,-58.4,47.06,47.14,47.18,47.26,0.0,58e14,1.18,95e14
992,15,-59.55,-59.16,-59.61,-59.21,46.99,47.09,47.11,47.2,0.0,77e14,1.05,166e14
992,16,-59.16,-58.77,-59.21,-58.83,47.09,47.18,47.2,47.29,0.0,10e15,1.02,-18e14
992,17,-58.77,-58.4,-58.83,-58.46,47.18,47.26,47.29,47.38,0.0,68e14,1.16,-78e14
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
    # Need to reconfigure to have nTimes_1 and nXtrack_1 which are
    # edge dimensions
    ds['PixelCornerLongitudes'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'LL_Longitude', 'UL_Longitude', 'UU_Longitude', 'LU_Longitude'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    # Need to reconfigure to have nTimes_1 and nXtrack_1 which are
    # edge dimensions
    ds['PixelCornerLongitudes'] = xr.DataArray(
        xr.concat([
            ds[key]
            for key in (
                'LL_Latitude', 'UL_Latitude', 'UU_Latitude', 'LU_Latitude'
            )
        ], dim='nCorners').transpose('nTimes', 'nXtrack', 'nCorners')
    )
    outds = ds
    # .drop_vars([
    #     'LL_Latitude', 'UL_Latitude', 'UU_Latitude', 'LU_Latitude',
    #     'LL_Longitude', 'UL_Longitude', 'UU_Longitude', 'LU_Longitude'
    # ])
    return outds


def test_omi():
    df = omi_example_df()
    ds = omi_example_ds()

    sat = readers.omi.OMHCHO.from_dataset(ds)
    df2 = sat.export_dataframe()
    print(df.columns)
    print(df2.columns)
    assert((df2 == df.loc[:, df2.columns]).all().all())


def test_omi_valid():
    df = omi_example_df()
    ds = omi_example_ds()
    ds['MainDataQualityFlag'][1, 1] = 1
    ds['MainDataQualityFlag'][0, 1] = 1
    ds['MainDataQualityFlag'][2, 2] = 1
    sat = readers.omi.OMHCHO.from_dataset(ds)
    df2 = sat.export_dataframe()
    assert((df.shape[0] - 3) == df2.shape[0])


def test_satellite_weights():
    from .. import cmaq

    ds = omi_example_ds()
    sat = readers.omi.OMHCHO.from_dataset(ds)

    cg = cmaq.CMAQGrid(None, '108US1')
    wgt = sat.weights(cg.geodf)
    wgtd1 = sat.weighted(
        'ReferenceSectorCorrectedVerticalColumn', wgtdf=wgt,
        groupkeys=['ROW', 'COL']
    )
    assert(
        wgtd1['ReferenceSectorCorrectedVerticalColumn'].min()
        >= ds['ReferenceSectorCorrectedVerticalColumn'].min()
    )
    assert(
        wgtd1['ReferenceSectorCorrectedVerticalColumn'].max()
        <= ds['ReferenceSectorCorrectedVerticalColumn'].max()
    )
    wgtd2 = sat.weighted(
        'ReferenceSectorCorrectedVerticalColumn', othdf=cg.geodf,
        groupkeys=['ROW', 'COL']
    )
    assert((wgtd1 == wgtd2).all().all())
