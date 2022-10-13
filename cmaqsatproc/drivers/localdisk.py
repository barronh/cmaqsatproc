import os
import pandas as pd
from ..import readers
from ..import cmaq


def localdisk(
    reader, gdnam, gdpath, inpaths, outformat, variables, geo, outpath,
    overwrite
):
    if outformat == 'l3':
        localdisk_l3(reader, gdnam, gdpath, inpaths, outpath, overwrite)
    elif outformat == 'csv':
        localdisk_csv(
            reader, gdnam, gdpath, inpaths, variables, geo, outpath, overwrite
        )
    else:
        raise KeyError('format must be l3 or csv')


def localdisk_l3(reader, gdnam, gdpath, inpaths, outpath, overwrite):
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
    inpaths : list
        List of files on disk
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
        outpath = f'{reader}_{gdnam}.nc'

    if os.path.exists(outpath):
        if overwrite:
            os.remove(outpath)
        else:
            raise IOError(f'{outpath} exists; use -O or --overwrite to force')

    cg = cmaq.CMAQGrid(GDNAM=gdnam, gdpath=gdpath)

    satreader = readers.reader_dict[reader]

    # Create Level3 from localdisk
    outputs = satreader.paths_to_level3(
        inpaths, bbox=cg.bbox(),
        grid=cg.geodf[['geometry']], as_dataset=True
    )

    outputs.to_netcdf(outpath)


def localdisk_csv(
    reader, gdnam, gdpath, inpaths, variables, geo, outpath, overwrite
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
    inpaths : list
        List of files on disk
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
            'localdisk csv requires variables; choose one or more of: '
            + ', '.join(satreader._defaultkeys)
        )

    if outpath is None:
        outpath = f'{reader}_{gdnam}.csv'

    if os.path.exists(outpath):
        if overwrite:
            os.remove(outpath)
        else:
            raise IOError(f'{outpath} exists; use -O or --overwrite to force')

    cg = cmaq.CMAQGrid(GDNAM=gdnam, gdpath=gdpath)

    # Create Level3 from localdisk
    if not geo:
        myvariables = list(satreader._geokeys) + list(myvariables)
    dfs = []
    for link in inpaths:
        sat = satreader.open_dataset(link, bbox=cg.bbox())
        df = sat.to_dataframe(*myvariables, geo=geo)
        dfs.append(df)
    paths = [link.split('/')[-1] for link in inpaths]

    df = pd.concat(dfs, keys=paths, names=['path'])
    df.to_csv(outpath)


def add_localdisk_parser(subparsers):
    localdiskparser = subparsers.add_parser(
        'localdisk',
        help=(
            'Create custom Level2 csv (csv) or Level3 (l3) file using data'
            + ' from files on your local disk.'
        )
    )
    localdiskparser.add_argument(
        '-O', '--overwrite', default=False, action='store_true',
        help='--outpath will be removed before running the command'
    )
    localdiskparser.add_argument(
        '--outpath', default=None,
        help='Defaults to {reader}_{GDNAM}.{suffix}'
    )
    localdiskparser.add_argument(
        '--gdpath', default=None, metavar='GRIDDESC',
        help='If provided, points to an IOAPI GRIDDESC file'
    )
    localdiskparser.add_argument(
        '-v', '--variables', default=[], action='append',
        help=(
            'Required for csv format; use comma separated variables or'
            + ' multiple --variables options.'
        )
    )
    localdiskparser.add_argument(
        '--geo', default=False, action='store_true',
        help=(
            'Only used by csv format; include polygon as WKT? If not, include'
            + ' pixel corners'
        )
    )
    localdiskparser.add_argument(
        'outformat', choices=('l3', 'csv'),
        help='Create Level2 csv or custom Level3 file directly from localdisk'
    )
    localdiskparser.add_argument(
        'reader', choices=tuple(readers.reader_dict.keys()), metavar='READER',
        help='Choose a reader from ' + ', '.join(readers.reader_dict.keys())
    )
    localdiskparser.add_argument(
        'gdnam', metavar='GDNAM',
        help=(
            'Name CMAQ Grid, either known by cmaqsatproc or in GRIDDESC'
            + ' (e.g., 12US1, 36US2, 108NHEMI1, global_1x1, global_0pt1deg,'
            + ' 108US1, US_1deg, etc)'
        )
    )
    localdiskparser.add_argument(
        'inpaths', nargs='+',
        help='In paths'
    )
