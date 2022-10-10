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
    assert (len(links) > 0)


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
    assert (len(links) > 0)


def test_easypoly():
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    df = pd.DataFrame(dict(
        ll_x=[0, -179.5, 0],
        lu_x=[0, -179.5, 0],
        uu_x=[1, 179.5, 1],
        ul_x=[1, 179.5, 1],
        ll_y=[0, 0, 89],
        lu_y=[1, 1, -89],
        uu_y=[1, 1, -89],
        ul_y=[0, 0, 89]
    ))
    opts = [dict(lowmem=True, wrap=False), dict(lowmem=False, wrap=True)]
    for opt in opts:
        gdf = gpd.GeoDataFrame(
            df, geometry=utils.EasyDataFramePolygon(df, **opt), crs=4326
        )
        x, y = gdf.iloc[0].geometry.exterior.xy
        assert (np.allclose(x, [0, 0, 1, 1, 0]))
        assert (np.allclose(y, [0, 1, 1, 0, 0]))
        x, y = gdf.iloc[1].geometry.exterior.xy
        if opt['wrap']:
            assert (np.allclose(x, [-179.5, -179.5, -180, -180, -179.5]))
        else:
            assert (np.allclose(x, [-179.5, -179.5, 179.5, 179.5, -179.5]))
        assert (np.allclose(y, [0, 1, 1, 0, 0]))
        x, y = gdf.iloc[2].geometry.exterior.xy
        assert (np.allclose(x, [0, 0, 1, 1, 0]))
        # Known failure. No reason a polar orbiting satellite should
        # simultaneously see both poles
        assert (np.allclose(y, [89, -89, -89, 89, 89]))


def test_rootremover():
    stem, short_list = utils.rootremover(
        ['/really/long/testing/1', '/really/long/testing/2']
    )
    assert (stem == '/really/long/testing')
    assert (short_list[0] == '{root}/1')
    assert (short_list[1] == '{root}/2')
