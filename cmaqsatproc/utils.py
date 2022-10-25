__all__ = [
    'getcmrlinks', 'getcmrgranules', 'rootremover', 'grouped_weighted_avg',
    'walk_groups', 'EasyDataFramePoint', 'EasyDataFramePolygon', 'cdconvert',
    'csp_version'
]


_ckeys = ['ll', 'lu', 'uu', 'ul']
_xcrnrkeys = [f'{ckey}_x' for ckey in _ckeys]
_ycrnrkeys = [f'{ckey}_y' for ckey in _ckeys]


def csp_version():
    from . import __version__ as _csp_version
    return _csp_version


def walk_groups(f, gkey, outputs=None):
    if outputs is None:
        outputs = []
    for key in sorted(f.groups.keys()):
        outputs.append(f'{gkey}/{key}')
        walk_groups(f[key], f'{gkey}/{key}', outputs)
    return outputs


def grouped_weighted_avg(values, weights, by):
    numerator = values.multiply(
        weights, axis=0
    ).groupby(by).sum()
    denominator = weights.groupby(by).sum()
    outdf = numerator.divide(denominator, axis=0)
    outdf['weight_sum'] = denominator
    return outdf


def getcmrgranules(
    short_name, temporal, bbox=None, poly=None, verbose=0, **kwds
):
    """
    Return all links from the Common Metadata Repository for the product
    granules with short_name, date_range, and optionally a bounding box.

    Arguments
    ---------
    short_name : str
        NASA short_name that is recognized by NASA Common Metadata Repository
    temporal : str
        NASA temporal that is recognized by NASA Common Metadata Repository
        i.e, YYYY-MM-DDTHH:MM:SSZ or YYYY-MM-DDTHH:MM:SSZ/YYYY-MM-DDTHH:MM:SSZ
        or YYYY-MM-DDTHH:MM:SSZ/P01D (or P01M or P01Y) etc.
    bbox : list
        Longitude/latitude bounding edges as floats (wlon,slat,elon,nlat)
    poly : shapely.geometry.Polygon
        The exterior will be converted to floats and ordered x1,y1,...,xN,yN
    filterfunc : function
        Takes a link dictionary from CMR and returns True if it should be
        retained

    Returns
    -------
    jr : dictionary
        json result as a dictionary
    """
    from copy import deepcopy
    import requests

    opts = deepcopy(kwds)
    opts.setdefault('page_size', '1000')
    opts.setdefault('pretty', 'false')
    opts['short_name'] = short_name
    opts['temporal'] = temporal

    if bbox is not None:
        opts['bounding_box'] = '{:f},{:f},{:f},{:f}'.format(*bbox)

    if poly is not None:
        import numpy as np
        from shapely.geometry import polygon
        pgstr = ','.join([
            f'{x:f},{y:f}'
            for x, y in zip(*np.asarray(polygon.orient(poly, 1).exterior.xy))
        ])
        opts['polygon'] = pgstr

    if verbose:
        print(opts)
    r = requests.get(
        'https://cmr.earthdata.nasa.gov/search/granules.json',
        params=opts
    )
    jr = r.json()
    return jr


def getcmrlinks(*args, filterfunc=None, **kwds):
    """
    Return all links from the Common Metadata Repository for the product
    granules with *args and **kwds are passed thru to getcmrgranules

    Arguments
    ---------
    *args, **kwds :
        See getcmrgranuels.
    filterfunc : function
        Takes a link dictionary from CMR and returns True if it should be
        retained

    Returns
    -------
    links : list
        List of all links where filterfunc is None or just links for which
        filterfunc returns True.
    """
    jr = getcmrgranules(*args, **kwds)
    entries = jr['feed']['entry']
    links = []
    for entry in entries:
        for link in filter(filterfunc, entry['links']):
            href = link['href']
            links.append(href)

    return links


def EasyDataFramePoint(df, xkey='cn_x', ykey='cn_y'):
    """
    Thin wrapper on geopandas.points_from_xy. Create points from a rows with
    centers named xkey and ykey.

    Arguments
    ---------
    df : pandas.DataFrame
        Must contain x, y

    Returns
    -------
    points : list
        List of shapely.geometry.Polygons
    """
    import geopandas as gpd
    points = gpd.points_from_xy(df[xkey].values, df[ykey].values)
    return points


def EasyDataFramePolygon(df, wrap=True, progress=False, lowmem=False):
    """
    Create polygons from a row with corners of a pixel specificied using
    columns ll_x, ll_y ... uu_x, uu_y.

    The wrap functionality prevents polygons from straddling the dateline.

    Arguments
    ---------
    df : pandas.DataFrame
        Must contain ll, lu, uu, ul for x and y (e.g., ll_x, ll_y)
    wrap : bool
        If True (default), each polygon that crosses the dateline will be
        truncated to the Western portion.

    Returns
    -------
    polys : list
        List of shapely.geometry.Polygons
    """
    from shapely.geometry import Polygon
    import numpy as np
    x = df[_xcrnrkeys].copy()
    y = df[_ycrnrkeys].copy()
    dx = x.max(axis=1) - x.min(axis=1)
    # Conceivably apply the same approach to the pole
    # dy = y.max(axis=1) - y.min(axis=1)
    if wrap:
        newx = x.where(~((dx > 90).values[:, None] & (x > 0).values), -180)
        newy = y  # .where(~((dy > 45).values[:, None] & (y < 0).values), 90)
    else:
        newx = x
        newy = y

    polys = []
    # Faster to convert all to values arrays first
    if lowmem:
        i = 0
        for idx, x in newx.iterrows():
            if progress:
                print(f'\r{i/newx.shape[0]:8.3%}', end='')
            y = newy.loc[idx]
            polys.append(Polygon(np.asarray([x, y]).T))
            i += 1
    else:
        xys = np.stack([newx.values, newy.values], axis=2)
        for i, xy in enumerate(xys):
            if progress:
                print(f'\r{i/xys.shape[0]:8.3%}', end='')
            polys.append(Polygon(xy))
    if progress:
        print('\r100.000%')

    return polys


def rootremover(strlist, insert=False):
    """
    Find the longest common root and replace it with {root}

    Arguments
    ---------
    strlist : list
        List of strings from which to find a common root and replace with
        '{root}'
    insert : bool
        If true, insert f'root: {root}' at the beginning of the short list.

    Return
    ------
    stem, short_list
        List with each element of strlist where the longest common root has
        been removed. If insert, then the root is inserted
    """
    import os

    stem = os.path.dirname(strlist[0])
    while not all([stem in l for l in strlist]) or False:
        oldstem = stem
        stem = os.path.dirname(stem)
        if oldstem == stem:
            short_strlist = strlist
            break
    else:
        short_strlist = [
            l.replace(stem, '{root}')
            for l in strlist
        ]
        if insert:
            short_strlist.insert(0, f'root: {stem}')

    return stem, short_strlist


def coord_interp(
    coordout, coordin, varin, verbose=0, interp='numpy',
    outdim='LAY', indim='LAY', ascending=True, **kwds
):
    """
    Arguments
    ---------
    coordout : xarray.DataArray
        n-dimensional (must have outdim) where values of coordout are the
        coordinate
    coordin : xarray.DataArray
        n-dimensional (must have indim) where values of coordin are the
        coordinate
    varin : xarray.DataArray
        n-dimensional (must have indim) where values will be interpolated
    **kwds : mappable
        supplied as keywords to numpy interp

    Returns
    -------
    out : xarray.DataArray
        varin interpolated to coordout
    """
    import xarray as xr
    import numpy as np

    if interp.lower() == 'numpy':
        def interp1d(data, x, xi, **kwds):
            return np.interp(xi, x, data, **kwds)
    else:
        def interp1d(data, x, xi, **kwds):
            from scipy import interpolate
            f = interpolate.interp1d(x, data, fill_value='extrapolate')
            return f(xi)
    tempdimname = 'temp_dim_name'
    X = coordout.rename(**{outdim: tempdimname})
    XP = coordin.sortby(indim, ascending=ascending)
    FP = varin.sortby(indim, ascending=ascending)
    interped = xr.apply_ufunc(
        interp1d,
        FP,
        XP,
        X,
        input_core_dims=[[indim], [indim], [tempdimname]],
        output_core_dims=[[tempdimname]],
        exclude_dims=set((indim,)),
        vectorize=True,
        kwargs=kwds
    )
    out = interped.rename(**{tempdimname: outdim})
    return out


def cdconvert(inval, inunit, outunit):
    """
    Converts between du, mole m**-2, and molecules cm**-2
    """

    inunit = inunit.strip().lower().replace('^', '**')
    outunit = outunit.strip().lower().replace('^', '**')

    if inunit == outunit:
        return inval
    elif inunit == 'du':
        if outunit == 'mole m**-2':
            return inval / 1e5 * 101325 / 8.314 / 273.15
        else:
            return cdconvert(
                cdconvert(inval, inunit, 'mole m**-2'), 'mole m**-2', outunit
            )
    elif inunit == 'molecules cm**-2':
        if outunit == 'mole m**-2':
            return inval / 6.022e23 * 1e4
        else:
            return cdconvert(
                cdconvert(inval, inunit, 'mole m**-2'), 'mole m**-2', outunit
            )
    elif inunit == 'mole m**-2':
        if outunit == 'du':
            return inval * (1 / 1e5 * 101325 / 8.314 / 273.15)**-1
        elif outunit == 'molecules cm**-2':
            return inval * 6.022e23 / 1e4
        else:
            raise KeyError(
                f'unknown outunit ({outunit}); use du, mole m**-2, or'
                + ' molecules cm**-2'
            )
    else:
        raise KeyError(
            f'unknown inunit ({inunit}); use du, mole m**-2, or'
            + ' molecules cm**-2'
        )
