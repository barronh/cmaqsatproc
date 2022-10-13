__all__ = [
    'satellite', 'goes', 'iasi', 'modis', 'omi', 'omps', 'tropomi', 'viirs',
    'reader_dict'
]


from .core import satellite
from . import goes
from . import iasi
from . import modis
from . import omi
from . import omps
from . import tropomi
from . import viirs
import inspect

reader_dict = {}
for submod in [goes, iasi, modis, omi, omps, tropomi, viirs]:
    readers = getattr(submod, '__all__', [])
    for _readerkey in readers:
        _reader = getattr(submod, _readerkey)
        if inspect.isclass(_reader) and issubclass(_reader, satellite):
            long_name = '.'.join(
                submod.__name__.split('.')[2:] + [_reader.__name__]
            )
            short_name = _reader.__name__
            reader_dict[short_name] = _reader
            reader_dict[long_name] = _reader

reader_dict['TropOMINO2'] = reader_dict['S5P_L2__NO2___']
reader_dict['TropOMICO'] = reader_dict['S5P_L2__CO____']
reader_dict['TropOMIHCHO'] = reader_dict['S5P_L2__HCHO__']
reader_dict['TropOMICH4'] = reader_dict['S5P_L2__CH4___']


def print_reader_list():
    """
    Print to the screen a short description reader_dict keys and a short
    description of the reader they are associated with.
    """
    for key, reader in reader_dict.items():
        print(f'{key}: ', reader.__doc__)
    print('For more help on any reader, use `help(reader_dict[key])`,')
    print('where key is the reader name. e.g., help(reader_dict["OMNO2"])')
