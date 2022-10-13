__all__ = ['opendap', 'localdisk', 'parser', 'run_from_args']

from . import opendap
from . import localdisk
import argparse


parser = argparse.ArgumentParser(prog='cmaqsatproc')
subparsers = parser.add_subparsers(
    dest='command', title='subcommands',
    description='Valid subcommands are show below:',
    help='For help on subcommands run %(prog)s subcommand -h'
)
opendap.add_opendap_parser(subparsers)
localdisk.add_localdisk_parser(subparsers)


def parse_args(args, run=True, noexit=True):
    """
    args : list
        Like argparse.ArgumentParser.parse_args (use '-h' for more details)
    run : bool
        If True, run the commands. Otherwise simply return the kwargs.
    noexit : bool
        By default, do not exits on error. If running from CLI, noexit should
        be False
    """
    from .opendap import opendap
    from .localdisk import localdisk
    _ = opendap
    _ = localdisk
    try:
        args = parser.parse_args(args)
        kwargs = vars(args)
        cmdname = kwargs.pop('command')
        if cmdname is None:
            parser.print_usage()
        else:
            if run:
                cmdfunc = eval(cmdname.replace('-', '_'))
                cmdfunc(**kwargs)
            else:
                return kwargs
    except SystemExit as e:
        if not noexit:
            raise e
        else:
            print(repr(e))
