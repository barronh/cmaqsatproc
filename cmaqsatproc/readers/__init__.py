__all__ = [
    'satellite', 'omi', 'modis', 'tropomi', 'iasi', 'viirs', 'reader_dict'
]


from .core import satellite
from . import omi
from . import modis
from . import tropomi
from . import iasi
from . import viirs
import inspect

reader_dict = {}
for submod in [omi, modis, tropomi, iasi, viirs]:
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
