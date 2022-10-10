from .. import cmaq


def test_cmaqgrid():
    import numpy as np
    cg = cmaq.open_griddesc(gdpath=None, GDNAM='108US3')
    bbox = cg.csp.bbox()
    assert (np.allclose(
        bbox,
        (
            -145.769499654764, 10.529366393866985,
            -40.46308679254433, 63.835501598090225
        )
    ))


def test_cmaqgrid_geodf():
    cg = cmaq.open_griddesc(gdpath=None, GDNAM='108US3')
    gdf = cg.csp.geodf
    assert (gdf.shape == (3000, 1))
    wkt = gdf.loc(axis=0)[0.5, 0.5].geometry.wkt
    assert (wkt == 'POLYGON ((1 0, 1 1, 0 1, 0 0, 1 0))')
    wkt = gdf.loc(axis=0)[49.5, 59.5].geometry.wkt
    assert (wkt == 'POLYGON ((60 49, 60 50, 59 50, 59 49, 60 49))')
