import argparse
from .drivers.opendap import add_opendap_parser, opendap
from .drivers.localdisk import add_localdisk_parser, localdisk

parser = opendap
parser = localdisk
parser = None

parser = argparse.ArgumentParser(prog='cmaqsatproc')
subparsers = parser.add_subparsers(
    dest='command', title='subcommands',
    description='valid subcommands',
    help='For help on subcommands run %(prog)s subcommand -h'
)
add_opendap_parser(subparsers)
add_localdisk_parser(subparsers)

# Define CMAQ (q), MCIP (m) and output (out) paths
args = parser.parse_args()

kwargs = vars(args)
cmdname = kwargs.pop('command')
if cmdname is None:
    parser.print_usage()
else:
    cmdfunc = eval(cmdname.replace('-', '_'))
    cmdfunc(**kwargs)
