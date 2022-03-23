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

    ds = xr.Dataset.from_dataframe(modis_example_df())
    outds = ds.drop_vars(['dx', 'dy'])
    return outds


def test_modis():
    df = modis_example_df()
    ds = modis_example_ds()

    sat = readers.modis.MOD04.from_dataset(ds)
    df2 = sat.export_dataframe()
    assert(
        (
            df2.loc[:, df.columns].drop(['dx', 'dy'], axis=1)
            == df.drop(['dx', 'dy'], axis=1)
        ).all().all()
    )


def test_modis_valid():
    df = modis_example_df()
    ds = modis_example_ds()
    ds['Land_Ocean_Quality_Flag'][1, 1] = 0
    ds['Land_Ocean_Quality_Flag'][0, 0] = 0
    ds['Land_Ocean_Quality_Flag'][2, 2] = 0
    sat = readers.modis.MOD04.from_dataset(ds)
    df2 = sat.export_dataframe()
    assert((df.shape[0] - 3) == df2.shape[0])


def test_satellite_weights():
    from .. import cmaq

    ds = modis_example_ds()
    sat = readers.modis.MOD04.from_dataset(ds)

    cg = cmaq.CMAQGrid(None, '108US1')
    wgt = sat.weights(cg.geodf)
    wgtd1 = sat.weighted(
        'Optical_Depth_Land_And_Ocean', wgtdf=wgt, groupkeys=['ROW', 'COL']
    )
    assert(
        wgtd1['Optical_Depth_Land_And_Ocean'].min()
        >= ds['Optical_Depth_Land_And_Ocean'].min()
    )
    assert(
        wgtd1['Optical_Depth_Land_And_Ocean'].max()
        <= ds['Optical_Depth_Land_And_Ocean'].max()
    )
    wgtd2 = sat.weighted(
        'Optical_Depth_Land_And_Ocean', othdf=cg.geodf,
        groupkeys=['ROW', 'COL']
    )
    assert((wgtd1 == wgtd2).all().all())
