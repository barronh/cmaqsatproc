__all__ = ['Pipeline', 'modis_pipeline']


from .pipeline_core import Pipeline


def get_readers(short_name):
    from .. import readers
    return {
        'MOD04_3K': readers.modis.MOD04,
        'MOD04_L2': readers.modis.MOD04,
        'MYD04_3K': readers.modis.MOD04,
        'MYD04_L2': readers.modis.MOD04,
        'OMHCHO': readers.omi.OMHCHO,
        'OMNO2': readers.omi.OMNO2,
        'S5P_L2__NO2___': readers.tropomi.TropOMI,
        'S5P_L2__HCHO__': readers.tropomi.TropOMI,
        'S5P_L2__CH4___': readers.tropomi.TropOMI,
        'S5P_L2__CO____': readers.tropomi.TropOMI,
    }[short_name]


def _omi_filter(links):
    return [l for l in links if 'opendap' in l and l.endswith('.he5')]


def _tropomi_filter(links):
    links = [l for l in links if l.endswith('.nc')]
    if any(['opendap' in l.lower() for l in links]):
        links = [l for l in links if 'opendap' in l.lower()]
    else:
        links = [l.replace('/data/', '/opendap/hyrax/') for l in links]
    return links


def _modis_filter(links):
    return [
        l[:-5] for l in links if 'opendap' in l and l.endswith('.hdf.html')
    ]


_link_filters = {
    'MOD04_3K': _modis_filter,
    'MOD04_L2': _modis_filter,
    'MYD04_3K': _modis_filter,
    'MYD04_L2': _modis_filter,
    'OMHCHO': _omi_filter,
    'OMNO2': _omi_filter,
    'S5P_L2__NO2___': _tropomi_filter,
    'S5P_L2__HCHO__': _tropomi_filter,
    'S5P_L2__CH4___': _tropomi_filter,
    'S5P_L2__CO____': _tropomi_filter,
}

_renamer2ds = {
    'S5P_L2__NO2___': {
        'PRODUCT_nitrogendioxide_tropospheric_column': 'VCDTROPNO2',
        'weights': 'weights'
    },
    'S5P_L2__HCHO__': {
        'PRODUCT_formaldehyde_tropospheric_vertical_column': 'VCDTROPHCHO',
        'weights': 'weights'
    },
    'S5P_L2__CH4___': {
        'PRODUCT_methane_mixing_ratio': 'VMRCH4',
        'weights': 'weights'
    },
    'S5P_L2__CO____': {
        'PRODUCT_carbonmonoxide_total_column': 'VCDTOTCO',
        'weights': 'weights'
    },
    'MOD04_3K': {
        'Optical_Depth_Land_And_Ocean': 'AOD_550', 'weights': 'weights'
    },
    'MOD04_L2': {
        'Optical_Depth_Land_And_Ocean': 'AOD_550', 'weights': 'weights'
    },
    'MYD04_3K': {
        'Optical_Depth_Land_And_Ocean': 'AOD_550', 'weights': 'weights'
    },
    'MYD04_L2': {
        'Optical_Depth_Land_And_Ocean': 'AOD_550', 'weights': 'weights'
    },
    'OMHCHO': {
        'AirMassFactor': 'AMF', 'ColumnUncertainty': 'VCD_SIGMA',
        'ReferenceSectorCorrectedVerticalColumn': 'VCD', 'weights': 'weights'
    },
    'OMNO2': {
        'ColumnAmountNO2Trop': 'VCDTROPNO2',
        'ColumnAmountNO2Strat': 'VCDSTRATNO2',
        'ColumnAmountNO2': 'VCDNO2',
        'SlantColumnAmountNO2': 'SCDNO2',
        'AmfTrop': 'AMFTROP', 'AmfStrat': 'AMFSTRAT'
    }
}

_renamer3ds = {
    # 'OMHCHO': {
    #     'ClimatologyLevels': 'PRES', 'ScatteringWeights': 'SCATWGT'
    # }
}


def get_pipeline(
    short_name, cmaqgrid, persist=True, output3d=False,
    link_filter=None, reader=None, renamer2d=None, renamer3d=None, local=False
):
    """
    Pipelines perform a series of steps, but creating them can be difficult.
    They can be difficult because there are an infinite set of possibile
    configurations. However, there are only a few common configurations.

    This function pairs short_name with the standard options associated with
    it.

    Arguments
    ---------
    short_name : str
        CMR short_name of satellite product
    cmaqgrid : str
        GRID name
    persist : bool
        Should data be output as IOAPI?
    output3d : bool
        Should 3d data be output as well
    link_filter : function
        If provided, override default link_filter
    reader : function
        If provided, override default reader
    renamer2d : function
        If provided, override default renamer2d
    renamer3d : function
        If provided, override default renamer3d
    local : bool
        If local, then files have been downloaded using `wget -r`

    Returns
    -------
    pipe : cmaqsatproc.pipeline.Pipeline
        Configured for standard application to short_name satellite product
    """
    if persist:
        outtmpl2d = None
        outtmpl3d = None
    else:
        outtmpl2d = False
        outtmpl3d = False

    if reader is None:
        reader = get_readers(short_name)

    if link_filter is not None:
        tmp_link_filter = link_filter
    else:
        tmp_link_filter = _link_filters[short_name]

    if local:
        def mylinkfilter(links):
            return [l.replace('https://', '') for l in tmp_link_filter(links)]
    else:
        link_filter = tmp_link_filter

    if renamer2d is None:
        renamer2d = _renamer2ds.get(short_name, None)

    if renamer2d is None:
        varkeys2d = None
    else:
        varkeys2d = [key for key in list(renamer2d) if key != 'weights']

    if output3d:
        if renamer3d is None:
            renamer3d = _renamer3ds.get(short_name, None)
    else:
        renamer3d = None

    if renamer3d is None:
        varkeys3d = None
    else:
        varkeys3d = list(renamer3d)

    return Pipeline(
        short_name=short_name, cmaqgrid=cmaqgrid,
        link_filter=link_filter, reader=reader,
        varkeys2d=varkeys2d, renamer2d=renamer2d, outtmpl2d=outtmpl2d,
        varkeys3d=varkeys3d, renamer3d=renamer3d, outtmpl3d=outtmpl3d,
    )
