from .. import utils


def _getmodislinks(**kwds):
    links = utils.getcmrlinks(
        short_name='MOD04_L2',
        temporal='2003-01-01T00:00:00Z/2003-01-01T23:59:59Z',
        **kwds
    )
    return links


def test_getcmrlinks_bbox():
    # Start by finding links for a day
    links = _getmodislinks(
        bbox=(-121, 21.6, -69.3, 52)
    )
    print(len(links))
    # Finds 100+ links. Depending on how you shape the polygon,
    assert(len(links) > 0)


def test_getcmrlinks_poly():
    from shapely import wkt
    poly = wkt.loads(
        'POLYGON ((-121.1 21.6, -69.3 20.6, -54.3 50.3,'
        + ' -93.6 56.6, -134.6 52.0, -121.1 21.6))'
    )
    links = _getmodislinks(
        poly=poly
    )
    # print(len(links))
    # Finds 120+ links. Depending on how you shape the polygon,
    assert(len(links) > 0)


def test_centertobox():
    wkt = utils.centertobox(0.5, 0.5, 1, 1).wkt
    assert(wkt == 'POLYGON ((1 0, 1 1, 0 1, 0 0, 1 0))')


def test_easypoly():
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    df = pd.DataFrame(dict(
        LL_Longitude=[0, -179.5, 0],
        LU_Longitude=[0, -179.5, 0],
        UU_Longitude=[1, 179.5, 1],
        UL_Longitude=[1, 179.5, 1],
        LL_Latitude=[0, 0, 89],
        LU_Latitude=[1, 1, -89],
        UU_Latitude=[1, 1, -89],
        UL_Latitude=[0, 0, 89]
    ))
    gdf = gpd.GeoDataFrame(
        df, geometry=df.apply(utils.EasyRowPolygon, axis=1), crs=4326
    )
    x, y = gdf.iloc[0].geometry.exterior.xy
    assert(np.allclose(x, [0, 0, 1, 1, 0]))
    assert(np.allclose(y, [0, 1, 1, 0, 0]))
    x, y = gdf.iloc[1].geometry.exterior.xy
    assert(np.allclose(x, [-179.5, -179.5, -180, -180, -179.5]))
    assert(np.allclose(y, [0, 1, 1, 0, 0]))
    x, y = gdf.iloc[2].geometry.exterior.xy
    assert(np.allclose(x, [0, 0, 1, 1, 0]))
    # Known failure. No reason a satellite should simultaneously see both poles
    assert(np.allclose(y, [89, -89, -89, 89, 89]))
    gdf = gpd.GeoDataFrame(
        df, geometry=utils.EasyDataFramePolygon(df), crs=4326
    )
    x, y = gdf.iloc[0].geometry.exterior.xy
    assert(np.allclose(x, [0, 0, 1, 1, 0]))
    assert(np.allclose(y, [0, 1, 1, 0, 0]))
    x, y = gdf.iloc[1].geometry.exterior.xy
    assert(np.allclose(x, [-179.5, -179.5, -180, -180, -179.5]))
    assert(np.allclose(y, [0, 1, 1, 0, 0]))
    x, y = gdf.iloc[2].geometry.exterior.xy
    assert(np.allclose(x, [0, 0, 1, 1, 0]))
    # Known failure. No reason a satellite should simultaneously see both poles
    assert(np.allclose(y, [89, -89, -89, 89, 89]))
