from __future__ import division, print_function

import os, sys
import logging
log = logging.getLogger('wemd_cli')

import argparse

# We must prefer to load the wemd package over this script;
# this is only a problem for a program called "wemd.py"
try:
    wemd_lib_path = os.environ['WEMDLIB']
except KeyError:
    wemd_lib_path = os.path.dirname(os.path.realpath(__file__)) + '/../src'
log.debug('prepending %r to sys.path' % wemd_lib_path)
sys.path.insert(0, wemd_lib_path)

import wemd
from wemd import rc

def cmd_init(args):
    from wemd.util.config_dict import ConfigDict
    sim_config_src = ConfigDict()
    try:
        sim_config_src.read_config_file(args.sim_config_file)
    except IOError as e:
        sys.stderr.write('cannot open simulation configuration file: %s\n' % e)
        sys.exit(rc.EX_ENVIRONMENT_ERROR)
        
    

    

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--rcfile', metavar='RCFILE', dest='run_config_file',
                    help='use RCFILE as the WEMD run-time configuration file (default: %s)' 
                          % rc.RC_DEFAULT_FILENAME)
parser.add_argument('-v', '--version', action='version', version='WEMD version %s' % wemd.version)

subparsers = parser.add_subparsers()

parser_init =    subparsers.add_parser('init', help='initialize a new simulation')
parser_init.add_argument('sim_config_file', metavar='SIM_CONFIG',
                         help='simulation configuration file')
parser_init.set_defaults(func=cmd_init)

parser_run =     subparsers.add_parser('run', help='start/continue a simulation')
parser_status =  subparsers.add_parser('status', help='report simulation status')
parser_console = subparsers.add_parser('console', help='console with current simulation loaded')


# Parse command line arguments
args = parser.parse_args()

# Read runtime configuration file
runtime_config = rc.read_config(args.run_config_file)

# Apply logging options
rc.configure_logging(runtime_config)

# Branch to appropriate function
args.func(args)
