from .. import utils


def _getmodislinks(**kwds):
    links = utils.getcmrlinks(
        short_name='MOD04_L2',
        temporal='2003-01-01T00:00:00Z/2003-01-01T23:59:59Z',
        **kwds
    )
    return links


def _test_getcmrlinks_bbox():
    # Start by finding links for a day
    links = _getmodislinks(
        bbox=(-121, 21.6, -69.3, 52)
    )
    print(len(links))
    # Finds 100+ links. Depending on how you shape the polygon,
    assert(len(links) > 0)


def _test_getcmrlinks_poly():
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
    assert(wkt == 'POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))')
