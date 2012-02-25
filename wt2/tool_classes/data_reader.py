from core import WEMDTool

import wemd

class WEMDDataReader(WEMDTool):
    def __init__(self):
        super(WEMDDataReader,self).__init__()
        self.data_manager = wemd.rc.get_data_manager() 
        self.wt2_wemd_h5name = None
        
    def add_args(self, parser):
        group = parser.add_argument_group('WEMD input data options')
        group.add_argument('-W', '--wemd-data', dest='wt2_wemd_h5name', metavar='WEMD_H5FILE',
                           help='''Take WEMD data from WEMD_H5FILE (default: read from the HDF5 file specified in wemd.cfg).''')
        
    def process_args(self, args):
        self.data_manager.wemd_h5name = self.wt2_wemd_h5name = args.wt2_wemd_h5name
