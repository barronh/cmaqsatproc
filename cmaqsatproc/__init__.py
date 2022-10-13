__all__ = [
    'utils', 'readers', 'cmaq', 'drivers', 'open_ioapi', 'open_griddesc',
    'reader_dict', 'print_reader_list'
]

from . import utils
from . import readers
from . import cmaq
from . import drivers

__version__ = '0.2.0'

reader_dict = readers.reader_dict
open_ioapi = cmaq.open_ioapi
open_griddesc = cmaq.open_griddesc


def print_reader_list():
    """
    Print to the screen a short description reader_dict keys and a short
    description of the reader they are associated with.
    """
    for key, reader in reader_dict.items():
        print(f'{key}: ', reader.__doc__)
    print('For more help on any reader, use `help(reader_dict[key])`,')
    print('where key is the reader name. e.g., help(reader_dict["OMNO2"])')
