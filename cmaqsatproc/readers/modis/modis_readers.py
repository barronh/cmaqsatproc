__all__ = ['gdalhdf4_to_xrdataset', 'gdal_to_xrdataarray']


def gdal_to_xrdataarray(path):
    import gdal
    import xarray as xr
    import numpy as np
    ds0 = gdal.Open(path)
    data = ds0.ReadAsArray()
    attrs = ds0.GetMetadata_Dict()
    attrs['crs_wkt'] = ds0.GetProjection()
    attrs['geo_transform'] = ds0.GetGeoTransform()
    if 'valid_range' in attrs:
        vmin, vmax = eval(attrs['valid_range'])
        data = np.ma.masked_greater(np.ma.masked_less(data, vmin), vmax)
    try:
        offset = np.asarray(attrs.get('add_offset', 0), dtype=data.dtype)
        scale = np.asarray(attrs.get('scale_factor', 1), dtype=data.dtype)
    except Exception:
        offset = np.asarray(attrs.get('add_offset', 0), dtype='d')
        scale = np.asarray(attrs.get('scale_factor', 1), dtype='d')
    finaldata = (
        offset
        + np.ma.masked_values(data, -28672)
        * scale
    )
    return xr.DataArray(
        finaldata,
        dims=('nCandidate', 'y', 'x'),
        attrs=attrs,
        name=path.split(':')[-1]
    )


def gdalhdf4_to_xrdataset(path, unpack=False):
    import pyproj
    import numpy as np
    from osgeo import gdal
    import xarray as xr
    ds = gdal.Open(path)
    sds = ds.GetSubDatasets()
    outds = xr.merge([
        gdal_to_xrdataarray(sdn)
        for sdn, defn in sds
        if (
            sdn.endswith('AOD_QA')
            or sdn.endswith('AOD_Uncertainty')
            or sdn.endswith('Optical_Depth_055')
        )
    ])
    attrs = outds.Optical_Depth_055.attrs
    nc = int(attrs['DATACOLUMNS'])
    nr = int(attrs['DATAROWS'])
    proj = pyproj.Proj(attrs['crs_wkt'])
    wx, sy = proj(
        float(attrs['WESTBOUNDINGCOORDINATE']),
        float(attrs['SOUTHBOUNDINGCOORDINATE'])
    )
    ex, ny = proj(
        float(attrs['EASTBOUNDINGCOORDINATE']),
        float(attrs['NORTHBOUNDINGCOORDINATE'])
    )
    GT = attrs['geo_transform']
    xbnds = np.linspace(wx, ex, nc + 1)
    ybnds = np.linspace(sy, ny, nr + 1)
    x = (xbnds[:-1] + xbnds[1:]) / 2
    y = (ybnds[:-1] + ybnds[1:]) / 2
    X_pixel = np.arange(nc)
    Y_line = np.arange(nr)
    x = GT[0] + X_pixel * GT[1] + Y_line * GT[2]
    y = GT[3] + X_pixel * GT[4] + Y_line * GT[5]
    outds.coords['x'] = x
    outds.coords['y'] = y
    outds['crs'] = xr.DataArray(
        0,
        dims=(),
        attrs=dict(crs_wkt=attrs['crs_wkt'])
    )
    outds['QualityLevel'] = xr.DataArray(
        (
            (outds['AOD_QA'].values.astype('<H') >> 8)
            & np.array(int('1111', 2), dtype='<H')
        ),
        dims=('nCandidate', 'y', 'x'),
        name='QualityLevel',
        attrs=dict(
            notes=(
                "see https://amt.copernicus.org/articles/11/5741/2018/"
                + "amt-11-5741-2018.pdf"
            ),
            long_name='QualityLevel', units='none',
            description=(
                '0: "best", 1: "Water sediment"; 2: "n/a", 3: "one cloud",'
                + ' 4: ">1 clouds", 5: "no retrieval",'
                + ' 6: "no retrieval near", 7: "Climatology AOD (0.02)",'
                + ' 8: "No retrieval glint",'
                + ' 9: "Retrieval low due to glint", 10: "AOD near coast",'
                + ' 11: "– Land, research quality: AOD retrieved but CM is'
                + 'possibly cloudy"'
            )
        )
    )
    if unpack:
        outds['AOD_QA_Bits'] = xr.DataArray(
            np.unpackbits(
                outds.AOD_QA.values.view('uint8'), bitorder='little'
            ).reshape(outds.AOD_QA.shape + (16,))[..., ::-1],
            dims=('nCandidate', 'y', 'x', 'bits'), name='QABits',
            attrs=dict(
                units='None',
                description=(
                    'count from right 15,14,13...,0\n'
                    + outds.AOD_QA.attrs['data description']
                )
            )
        )

    outds.coords['x_bounds'] = xbnds
    outds.coords['y_bounds'] = ybnds
    return outds
