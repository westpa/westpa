from __future__ import print_function, division; __metaclass__=type
import os, argparse
import numpy
import wemd

import logging, warnings
log = logging.getLogger('wbf_import')

from wemdtools.aframe import (WEMDAnalysisTool,ExtDataReaderMixin,BFDataManager)


class WBFImporter(ExtDataReaderMixin,BFDataManager,WEMDAnalysisTool):
    def __init__(self):
        super(WBFImporter,self).__init__()
        self.bf_mode = True
        
    def import_trajectories(self):
        '''Import one or more trajectories into the brute force HDF5 file.'''
        
        for filename in self.ext_input_filenames:
            traj_id, traj_group = self.create_traj_group()
            wemd.rc.pstatus('Reading "{:s}" and storing as trajectory {:d}...'.format(filename, traj_id))
            if self.is_npy(filename):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.npy_to_h5dataset(numpy.load(filename, 'r'), traj_group, 'pcoord', usecols=self.ext_input_usecols)
            else:
                self.text_to_h5dataset(filename, traj_group, 'pcoord', dtype=numpy.float64, usecols=self.ext_input_usecols)
            self.update_traj_index(traj_id, traj_group['pcoord'].len(), source_data=os.path.abspath(filename))
        
    def main(self):
        parser = argparse.ArgumentParser('wbf_import', description='''\
        Import brute force trajectory progress coordinate data into an HDF5 file for subsequent analysis. 
        ''')
        wemd.rc.add_args(parser)
        self.add_args(parser)
        args = parser.parse_args()
        wemd.rc.process_args(args, config_required=False)
        
        self.process_args(args)
        self.require_bf_h5file()
        self.import_trajectories()

