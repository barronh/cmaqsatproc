__all__ = ['getcmrlinks', 'centertobox']


def getcmrlinks(short_name, temporal, bbox=None, poly=None, verbose=0, **kwds):
    """
    Return all links from the Common Metadata Repository
    for the product with short_name, date_range, and optionally
    a bounding box.

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
    entries = jr['feed']['entry']
    links = []
    for entry in entries:
        for link in entry['links']:
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
    from shapely.geometry import Polygon

    hdx = width / 2
    hdy = height / 2
    return Polygon([
        [xc - hdx, yc - hdy],
        [xc + hdx, yc - hdy],
        [xc + hdx, yc + hdy],
        [xc - hdx, yc + hdy],
        [xc - hdx, yc - hdy],
    ])
