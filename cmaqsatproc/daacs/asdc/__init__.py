__all__ = ['s3access_asdc']
from .._utils import s3access


class s3access_asdc(s3access):
    def __init__(
        self, url='https://data.asdc.earthdata.nasa.gov/s3credentials',
        buffer='2min', verbose=0
    ):
        super().__init__(credurl=url, buffer=buffer, verbose=verbose)
