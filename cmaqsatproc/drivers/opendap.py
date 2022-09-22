import os
import pandas as pd
from ..import readers
from ..import cmaq


def opendap(
    reader, gdnam, gdpath, SDATE, EDATE, outformat, variables, geo, outpath,
    overwrite
):
    if outformat == 'l3':
        opendap_l3(reader, gdnam, gdpath, SDATE, EDATE, outpath, overwrite)
    elif outformat == 'csv':
        opendap_csv(
            reader, gdnam, gdpath, SDATE, EDATE, variables, geo, outpath,
            overwrite
        )
    else:
        raise KeyError('format must be l3 or csv')


def opendap_l3(reader, gdnam, gdpath, SDATE, EDATE, outpath, overwrite):
    """
    Write a Level3 NetCDF file to disk.

    Arguments
    ---------
    reader : str
        Key for any satellite reader (e.g., 'OMNO2'). See
        cmaqsatproc.reader.reader_dict for options.
    gdnam : str
        Name of grid (either known or in gdpath)
    gdpath : str or None
        if str, a path to a GRIDDESC file.
    SDATE : str
        Any date-like starting date for CMR query
    EDATE : str
        Any date-like starting date for CMR query
    outpath : str
        path to write out result; If None, default {reader}_{gdnam}.csv
    overwrite : bool
        Overwrite existing files?

    Results
    -------
    None :
        Outputs are written to disk.
    """

    if outpath is None:
        outpath = f'{reader}_{SDATE}_{EDATE}_{gdnam}.nc'

    if os.path.exists(outpath):
        if overwrite:
            os.remove(outpath)
        else:
            raise IOError(f'{outpath} exists; use -O or --overwrite to force')

    cg = cmaq.CMAQGrid(GDNAM=gdnam, gdpath=gdpath)

    sdate = pd.to_datetime(SDATE)
    edate = pd.to_datetime(EDATE)

    satreader = readers.reader_dict[reader]

    # Create Level3 from OpenDAP
    timerange = f'{sdate:%FT%H:%M:%S}Z/{edate:%FT%H:%M:%S}Z'
    outputs = satreader.cmr_to_level3(
        temporal=timerange, bbox=cg.bbox(),
        grid=cg.geodf[['geometry']], as_dataset=True
    )

    outputs.to_netcdf(outpath)


def opendap_csv(
    reader, gdnam, gdpath, SDATE, EDATE, variables, geo, outpath, overwrite
):
    """
    Write a Level2 csv file to disk.

    Arguments
    ---------
    reader : str
        Key for any satellite reader (e.g., 'OMNO2'). See
        cmaqsatproc.reader.reader_dict for options.
    gdnam : str
        Name of grid (either known or in gdpath)
    gdpath : str or None
        if str, a path to a GRIDDESC file.
    SDATE : str
        Any date-like starting date for CMR query
    EDATE : str
        Any date-like starting date for CMR query
    variables : list
        Variables to write out.
    geo : bool
        Write geometry as a WKT polygon
    outpath : str
        path to write out result; If None, default {reader}_{gdnam}.csv
    overwrite : bool
        Overwrite existing files?

    Results
    -------
    None :
        Outputs are written to disk.
    """
    myvariables = []
    for varkey in variables:
        if ',' in varkey:
            myvariables.extend(varkey.split(','))
        else:
            myvariables.append(varkey)

    satreader = readers.reader_dict[reader]
    if len(myvariables) < 1:
        raise ValueError(
            'opendap csv requires variables; choose one or more of: '
            + ', '.join(satreader._defaultkeys)
        )

    if outpath is None:
        outpath = f'{reader}_{SDATE}_{EDATE}_{gdnam}.csv'

    if os.path.exists(outpath):
        if overwrite:
            os.remove(outpath)
        else:
            raise IOError(f'{outpath} exists; use -O or --overwrite to force')

    cg = cmaq.CMAQGrid(GDNAM=gdnam, gdpath=gdpath)

    sdate = pd.to_datetime(SDATE)
    edate = pd.to_datetime(EDATE)

    # Create Level3 from OpenDAP
    timerange = f'{sdate:%FT%H:%M:%S}Z/{edate:%FT%H:%M:%S}Z'
    links = satreader.cmr_links(
        temporal=timerange, bbox=cg.bbox(),
    )
    dfs = []
    if not geo:
        myvariables = list(satreader._geokeys) + list(myvariables)
    for link in links:
        sat = satreader.open_dataset(link, bbox=cg.bbox())
        df = sat.to_dataframe(*myvariables, geo=geo)
        dfs.append(df)
    paths = [link.split('/')[-1] for link in links]

    df = pd.concat(dfs, keys=paths, names=['path'])
    df.to_csv(outpath)


def add_opendap_parser(subparsers):
    opendapparser = subparsers.add_parser(
        'opendap',
        help=(
            'Create custom Level2 csv (csv) or Level3 (l3) file using data from'
            + ' remote OpenDAP files. Requires .netrc and .dodsrc'
            + ' configuration (https://disc.gsfc.nasa.gov/data-access).'
        )
    )
    opendapparser.add_argument(
        '-O', '--overwrite', default=False, action='store_true',
        help='--outpath will be removed before running the command'
    )
    opendapparser.add_argument(
        '--outpath', default=None,
        help='Defaults to {reader}_YYYY-MM-DD_YYYY-MM-DD_{GDNAM}.{suffix}'
    )
    opendapparser.add_argument(
        '--gdpath', default=None, metavar='GRIDDESC',
        help='If provided, points to an IOAPI GRIDDESC file'
    )
    opendapparser.add_argument(
        '-v', '--variables', default=[], action='append',
        help=(
            'Required for csv format; use comma separated variables or multiple'
            + ' --variables options.'
        )
    )
    opendapparser.add_argument(
        '--geo', default=False, action='store_true',
        help=(
            'Only used by csv format; include polygon as WKT? If not, include'
            + ' pixel corners'
        )
    )
    opendapparser.add_argument(
        'outformat', choices=('l3', 'csv'),
        help='Create Level2 csv or custom Level3 file directly from OpenDAP'
    )
    opendapparser.add_argument(
        'reader', choices=tuple(readers.reader_dict.keys()), metavar='READER',
        help=(
            'Choose a satellite READER: '
            + ', '.join(readers.reader_dict.keys())
        )
    )
    opendapparser.add_argument(
        'gdnam', metavar='GDNAM',
        help='Name CMAQ Grid, either known by cmaqsatproc or in GRIDDESC'
    )
    opendapparser.add_argument(
        'SDATE', metavar='START_DATE',
        help='Start date CMR query'
    )
    opendapparser.add_argument(
        'EDATE', metavar='END_DATE',
        help='End date CMR query'
    )
