from __future__ import print_function, division; __metaclass__ = type

import wemd, h5py

class WEMDTool:
    def __init__(self):
        self.wt2_config_required = False
    
    def add_args(self, parser):
        '''Add arguments specific to this tool to the given argparse parser.'''
        wemd.rc.add_args(parser)
    
    def process_args(self, args):
        '''Take argparse-processed arguments associated with this tool and deal
        with them appropriately (setting instance variables, etc)'''
        wemd.rc.process_args(args, config_required = self.wt2_config_required)
    
class HDF5Storage(WEMDTool):
    '''Class/mixin for storing data in an HDF5 file.'''
    def __init__(self):
        self.wt2_analysis_h5filename = None
        self.analysis_h5file = None
        
    def add_args(self, parser):
        group = parser.add_argument_group('analysis storage options')
        group.add_argument('-A', '--analysis-storage', dest='wt2_analysis_h5filename', metavar='HDF5FILE',
                           help='''Store results in HDF5FILE.''')
    
    def process_args(self, args):
        self.wt2_analysis_h5filename = args.wt2_analysis_h5filename
        
    def open_analysis_h5file(self, mode=None):
        '''Open the analysis HDF5 file with the given ``mode``.'''
        self.analysis_h5file = h5py.File(self.wt2_analysis_h5filename, mode=mode)
    
    def require_analysis_group(self, path, replace=False):
        '''Ensure that the given group exists in the analysis HDF5 file, optionally
        replacing (deleting and recreating) it.'''
        if not self.analysis_h5file:
            self.open_analysis_h5file()
            
        if replace:
            try:
                del self.analysis_h5file[path]
            except KeyError:
                pass
            return self.analysis_h5file.create_group(path)
        else:
            return self.analysis_h5file.require_group(path)
                
        