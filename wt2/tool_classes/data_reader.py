from core import WEMDTool

import wemd

class WEMDDataReader(WEMDTool):
    def __init__(self):
        super(WEMDDataReader,self).__init__()
        self.data_manager = wemd.rc.get_data_manager() 
        self.we_h5filename = None
        
    def add_args(self, parser):
        group = parser.add_argument_group('WEMD input data options')
        group.add_argument('-W', '--wemd-data', dest='we_h5filename', metavar='WEMD_H5FILE',
                           help='''Take WEMD data from WEMD_H5FILE (default: read from the HDF5 file specified in wemd.cfg).''')
        
    def process_args(self, args):
        if args.we_h5filename:
            self.data_manager.we_h5filename = self.we_h5filename = args.we_h5filename
        else:
            self.we_h5filename = self.data_manager.we_h5filename
        
    def open(self, mode):
        self.data_manager.open_backing(mode)
        
    def close(self):
        self.data_manager.close_backing()
