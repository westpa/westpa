from __future__ import division, print_function; __metaclass__ = type
from core import WESTTool
import westpa
    
class WESTDataReader(WESTTool):
    '''Tool for reading data from WEST-related HDF5 files. Coordinates finding
    the main HDF5 file from west.cfg or command line arguments, caching of certain
    kinds of data (eventually), and retrieving auxiliary data sets from various
    places.'''
    
    def __init__(self):
        super(WESTDataReader,self).__init__()
        self.data_manager = westpa.rc.get_data_manager() 
        self.we_h5filename = None
        
    def add_args(self, parser):
        group = parser.add_argument_group('WEST input data options')
        group.add_argument('-W', '--west-data', dest='we_h5filename', metavar='WEST_H5FILE',
                           help='''Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in west.cfg).''')
        
    def process_args(self, args):
        if args.we_h5filename:
            self.data_manager.we_h5filename = self.we_h5filename = args.we_h5filename
        else:
            self.we_h5filename = self.data_manager.we_h5filename
        
    def open(self, mode='r'):
        self.data_manager.open_backing(mode)
        
    def close(self):
        self.data_manager.close_backing()
        
    def __getattr__(self, key):
        return getattr(self.data_manager, key)
