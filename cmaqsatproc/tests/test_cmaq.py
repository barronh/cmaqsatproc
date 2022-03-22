from .. import cmaq
import warnings


def test_cmaqgrid():
    import numpy as np
    cg = cmaq.CMAQGrid(None, '108US1')
    bbox = cg.bbox()
    assert(np.allclose(
        bbox,
        (
            -145.769499654764, 10.529366393866985,
            -40.46308679254433, 63.835501598090225
        )
    ))

def test_cmaqgrid_template():
    cg = cmaq.CMAQGrid(None, '108US1')
    tmpl = cg.template(
        ['NO2', 'SO2'], vglvls=(1, 0), vgtop=5000, ntsteps=1, VGTYP=7,
        SDATE=1970015, STIME=140000
    )
    assert('NO2' in tmpl.variables)
    assert('SO2' in tmpl.variables)
    assert(tmpl.VGTOP == 5000)
    assert(len(tmpl.dimensions['LAY']) == 1)
    assert(len(tmpl.dimensions['VAR']) == 2)
    assert(len(tmpl.dimensions['TSTEP']) == 1)
    assert(tmpl.SDATE == 1970015)
    assert(tmpl.STIME == 140000)


def test_cmaqgrid_geodf():
    cg = cmaq.CMAQGrid(None, '108US1')
    gdf = cg.geodf
    assert(gdf.shape == (3000, 1))
    wkt = gdf.loc(axis=0)[0, 0].geometry.wkt
    assert(wkt == 'POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))')
    wkt = gdf.loc(axis=0)[49, 59].geometry.wkt
    assert(wkt == 'POLYGON ((59 49, 60 49, 60 50, 59 50, 59 49))')


def test_cmaq_to_ioapi():
    import pandas as pd
    import numpy as np
    cg = cmaq.CMAQGrid(None, '108US1')
    df = pd.DataFrame.from_dict(dict(
        ROW=[0, 49], COL=[0, 59], NAMEISTOOLONGNO2=[0, .5]
    )).set_index(['ROW', 'COL'])
    df['NAMEISTOOLONGNO2'].attrs['units'] = 'DU'
    df['NAMEISTOOLONGNO2'].attrs['long_name'] = 'NO2'
    df['NAMEISTOOLONGNO2'].attrs['var_desc'] = 'NO2'
    renamer = {
        'NAMEISTOOLONGNO2': 'NO2'
    }
    iof = cg.to_ioapi(
        df, TSTEP=0, LAY=0, ROW='ROW', COL='COL', rename=renamer
    )
    assert('NO2' in iof.variables)
    vno2 = iof.variables['NO2']
    assert(vno2[0, 0, 0, 0] == 0)
    assert(vno2[0, 0, 49, 59] == 0.5)
    nmasked = iof.variables['NO2'].mask.sum()
    assert((np.prod(vno2.shape) - 2) == nmasked)
    assert(vno2.units.strip() == 'DU')