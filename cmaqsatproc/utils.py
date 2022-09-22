__all__ = [
    'getcmrlinks', 'getcmrgranules', 'centertobox', 'cornertobox',
    'EasyRowPolygon', 'weight_vars', 'rootremover', 'csp_formatwarning',
    'csp_formatwarnings', 'grouped_weighted_avg', 'row_to_poly', 'walk_groups'
]
import warnings


_ckeys = ['ll', 'lu', 'uu', 'ul']
_xcrnrkeys = [f'{ckey}_x' for ckey in _ckeys]
_ycrnrkeys = [f'{ckey}_y' for ckey in _ckeys]


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


def row_to_poly(row, wrap=False):
    from shapely.geometry import Polygon
    import numpy as np
    coords = np.asarray([
        (row['ll_x'], row['ll_y']),
        (row['lu_x'], row['lu_y']),
        (row['uu_x'], row['uu_y']),
        (row['ul_x'], row['ul_y']),
        (row['ll_x'], row['ll_y']),
    ])
    if wrap:
        minx = coords[:, 0].min()
        maxx = coords[:, 0].mix()
        dx = maxx - minx
        if dx > 90:
            x = coords[:, 0]
            coords[:, 0] = np.where(x > 0, -180, x)

    return Polygon(coords)


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


def centertobox(xc, yc, width, height):
    """
    Return a polygon using a fixed width (2*hdx) and height (2*hdy)
    where the corners are offset from the center by hdx and hdy

    Arguments
    ---------
    xc, yc : scalar
        x,y-coordinates for the center of the polygon
    width, height : scalar
        Width and height (Half on each side)

    Returns
    -------
    poly : shapely.geometry.Polygon
        Polygon with exterior points at
        x = xc - hdx, x + hdx, x + hdx, x - hdx, x - hdx
        y = yc - hdy, y + hdy, y + hdy, y - hdy, y - hdy
    """
    from shapely.geometry import box

    hdx = width / 2
    hdy = height / 2
    return box(xc - hdx, yc - hdy, xc + hdx, yc + hdy)


def cornertobox(xc, yc, width, height):
    """
    Return a polygon using a fixed width (2*hdx) and height (2*hdy)
    where the corners are offset from the center by hdx and hdy

    Arguments
    ---------
    xc, yc : scalar
        x,y-coordinates for the corner of the polygon
    width, height : scalar
        Width and height

    Returns
    -------
    poly : shapely.geometry.Polygon
        Polygon with exterior points at
        x = xc, x + width, x + width, x, x
        y = yc, y + height, y + height, y, y
    """
    from shapely.geometry import box
    return box(xc, yc, xc + width, yc + height)


def weight_vars(withweights, *aggkeys, groupkeys=('ROW', 'COL')):
    """
    Arguments
    ---------
    withweights : pandas.DataFrame
    aggkeys : keys
    groupkeys : iterable

    Returns
    -------
    agg : pandas.DataFrame
    """
    notweighted = ['weights'] + list(groupkeys)
    aggattrs = {
        key: withweights[key].attrs
        for key in withweights.columns
    }
    dropkeys = ['geometry', 'dest_geometry']
    dropkeys = [k for k in dropkeys if k in withweights.columns]
    weighted = withweights.drop(
        dropkeys, axis=1, inplace=False
    ).multiply(withweights.weights, axis=0)
    # overwrite weights with original values
    for key in notweighted:
        if key not in weighted.index.names:
            weighted[key] = withweights[key]
    # Groupby and sum
    aggweighted = weighted.groupby(groupkeys).sum()
    # Divide weighted values by the sum of weights
    agg = aggweighted.divide(aggweighted.weights, axis=0)
    # Overwrite with the sum of weights
    for key in notweighted:
        if key not in groupkeys:
            agg[key] = aggweighted[key]

    for key in agg.columns:
        if key in aggattrs:
            agg[key].attrs.update(aggattrs[key])

    return agg


def EasyRowPolygon(row, wrap=True):
    """
    Create polygons from a row with corners of a pixel specificied using
    columns LL_Longitude, LL_Latitude... UU_Longitude, UU_Latitude.

    The wrap functionality prevents polygons from straddling the dateline.

    Arguments
    ---------
    row : pandas.DataFrame row
        Must contain LL, LU, UU, UL for Longitude and Latitude (e.g.,
        LL_Longitude, LL_Latitude)
    wrap : bool
        If True (default), each polygon that crosses the dateline will be
        truncated to the Western portion.

    Returns
    -------
    poly : shapely.geometry.Polygon
    """
    from shapely.geometry import Polygon
    import numpy as np

    x = row[_xcrnrkeys].values
    y = row[_ycrnrkeys].values
    if wrap:
        dx = x.max() - x.min()
        if dx > 90:
            x = np.where(x > 0, -180, x)

    # Conceivably apply the same approach to the pole
    # dy = y.max() - y.min()
    # if dy > 45:
    #     y = np.where(y < 0, 90, y)
    return Polygon(np.asarray([x, y]).T)


def EasyDataFramePoint(df, wrap=True, progress=False, lowmem=False):
    """
    Create polygons from a row with corners of a pixel specificied using
    columns LL_Longitude, LL_Latitude... UU_Longitude, UU_Latitude.

    The wrap functionality prevents polygons from straddling the dateline.

    (Same functionality as EasyRowPolygon, but intended to be faster due to
    the ability to rapidly apply wrapping to multiple rows at a time.)

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
    points = gpd.points_from_xy(df['x'].values, df['y'].values)
    return points


def EasyDataFramePolygon(df, wrap=True, progress=False, lowmem=False):
    """
    Create polygons from a row with corners of a pixel specificied using
    columns ll_x, ll_y ... uu_x, uu_y.

    The wrap functionality prevents polygons from straddling the dateline.

    (Same functionality as EasyRowPolygon, but intended to be faster due to
    the ability to rapidly apply wrapping to multiple rows at a time.)

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


def csp_formatwarning(message, category, filename, lineno, line=None):
    """
    Simplify warnings that come from cmaqsatproc. All others will remain in
    standard format.
    """
    if 'cmaqsatproc' in filename:
        filename = 'cmaqsatproc' + filename.split('cmaqsatproc')[-1]
        if category in (UserWarning, RuntimeWarning):
            category = 'Note:'
        else:
            category = str(category).split("'")[1] + ':'
        warnstr = f'{filename}:{lineno}:{category} {message}\n'
    else:
        warnstr = warnings.formatwarning(
            message=message, category=category,
            filename=filename, lineno=lineno,
            line=line
        )
    return warnstr


def csp_formatwarnings(warning_list, asstr=True):
    """
    Apply csp_formatwarnings to a list of warnings. If asstr is true,
    concatenate the results.
    """
    fmt_warning_list = []
    for w in warning_list:
        fmt_warning_list.append(
            csp_formatwarning(
                message=w.message, category=w.category,
                filename=w.filename, lineno=w.lineno, line=w.line
            )
        )
    if asstr:
        return ''.join(fmt_warning_list)
    else:
        return fmt_warning_list


def csp_showwarning(message, category, filename, lineno, file=None, line=None):
    """
    Could be used to overwrite warnings.showwarning. For any other modules,
    warnings will still be shown as normal. For cmaqsatproc warnings, a special
    format is applied
    """
    import sys
    if file is None:
        file = sys.stdout
    warnstr = csp_formatwarning(message, category, filename, lineno, line=line)
    file.write(warnstr, flush=True)


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
