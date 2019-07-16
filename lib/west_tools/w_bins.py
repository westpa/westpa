# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


import logging, sys
import numpy 
log = logging.getLogger('w_bins')

from westtools import WESTTool, WESTDataReader, BinMappingComponent
import westpa

from westtools.binning import write_bin_info

class WBinTool(WESTTool):
    prog='w_bins'
    description = '''\
Display information and statistics about binning in a WEST simulation, or
modify the binning for the current iteration of a WEST simulation.        
-------------------------------------------------------------------------------
'''    
    def __init__(self):
        super(WBinTool,self).__init__()
        self.subcommand = None
        self.data_reader = WESTDataReader() 
        self.binning = BinMappingComponent()
        self.args = None
        self.n_iter = None
                
    # Interface for command-line tools
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        
        subparsers = parser.add_subparsers(help='available commands')
        
        info_parser = subparsers.add_parser('info', help='Display information about binning.')
        info_parser.add_argument('-n', '--n-iter', type=int, 
                                 help='''Consider initial points of segment N_ITER (default: current iteration).''')
        info_parser.add_argument('--detail', action='store_true',
                                 help='''Display detailed per-bin information in addition to summary
                                 information.''')
        self.binning.add_args(info_parser)
        info_parser.set_defaults(func=self.cmd_info)
        
        rebin_parser = subparsers.add_parser('rebin',help='Rebuild current iteration with new binning.')
        rebin_parser.add_argument('--confirm', action='store_true', 
                                  help='''Commit the revised iteration to HDF5; without this option, the effects of the
                                  new binning are only calculated and printed.''')
        rebin_parser.add_argument('--detail', action='store_true',
                                  help='''Display detailed per-bin information in addition to summary
                                     information.''')
        self.binning.add_args(rebin_parser, suppress=['--bins-from-file'])
        self.binning.add_target_count_args(rebin_parser)
        rebin_parser.set_defaults(func=self.cmd_rebin)
            
    def process_args(self, args):

        self.data_reader.process_args(args)
        self.data_reader.open(mode='r+')
        self.n_iter = getattr(args,'n_iter',None) or self.data_reader.current_iteration
        
        # we cannot read bin information during rebins
        # interesting note: '==' is required here; 'is' fails
        if args.func == self.cmd_rebin:
            self.binning.target_counts_required = True
        else:
            self.binning.set_we_h5file_info(self.n_iter, self.data_reader)    
        
        self.binning.process_args(args)
        
        self.args = args
        self.subcommand = args.func
        
    def go(self):
        self.subcommand()
        
    def cmd_info(self):
        mapper = self.binning.mapper
        
        # Get target states and their assignments
        target_states = self.data_reader.get_target_states(self.n_iter)
        n_target_states = len(target_states)
        
        iter_group = self.data_reader.get_iter_group(self.n_iter)
        
        # bin initial pcoords for iteration n_iter
        initial_pcoords = iter_group['pcoord'][:,0,:]
        assignments = mapper.assign(initial_pcoords)
        del initial_pcoords
        
        print('Bin information for iteration {:d}'.format(self.n_iter))
        
        # Get bin counts and weights
        weights = iter_group['seg_index']['weight']
        
        write_bin_info(mapper, assignments, weights, n_target_states, detailed=self.args.detail)
            
    def cmd_rebin(self):
        mapper = self.binning.mapper
        assert mapper is not None    
        if self.n_iter == 1:
            sys.stderr.write('rebin is not supported for the first iteration; reinitialize with w_init instead\n')
            sys.exit(1)
        n_target_states = len(self.data_reader.get_target_states(self.n_iter))
        we_driver = westpa.rc.get_we_driver()
        data_manager = self.data_reader.data_manager
        
        segments = data_manager.get_segments(self.n_iter,load_pcoords=True)
        last_iter_segments = data_manager.get_segments(self.n_iter-1,load_pcoords=False,load_auxdata=False)
                
        # Bin on this iteration's initial points
        # We don't have to worry about recycling because we are binning on
        # initial points rather than final points, so recycling has already
        # occurred for this iteration.
        # We do need initial states, in case we merge a newly-created walker out of existence
        #avail_initial_states = {state.state_id: state
        #                        for state in data_manager.get_unused_initial_states(n_iter = self.n_iter)}
        avail_initial_states = data_manager.get_unused_initial_states(n_iter = self.n_iter)
        used_initial_states = data_manager.get_segment_initial_states(segments)
        we_driver.new_iteration(initial_states=avail_initial_states,
                                bin_mapper=mapper, bin_target_counts=self.binning.bin_target_counts)
        we_driver.used_initial_states = {state.state_id: state for state in used_initial_states}
        we_driver.assign(segments,initializing=True)
        we_driver.rebin_current(parent_segments=last_iter_segments)
        
        weights = numpy.array([segment.weight for segment in we_driver.next_iter_segments])
        assignments = numpy.fromiter(we_driver.next_iter_assignments,dtype=int,count=len(weights))
        write_bin_info(mapper, assignments, weights, n_target_states, detailed=self.args.detail)
        
        if self.args.confirm:
            data_manager.prepare_iteration(self.n_iter, list(we_driver.next_iter_segments))
            
            # manually update endpoint statuses only
            endpoint_types = sorted([(segment.seg_id, segment.endpoint_type) for segment in last_iter_segments])
            last_iter_group = data_manager.get_iter_group(self.n_iter-1)
            last_iter_index = last_iter_group['seg_index'][...]
            last_iter_index['endpoint_type'] = [pair[1] for pair in endpoint_types]
            last_iter_group['seg_index'][...] = last_iter_index
            
            data_manager.save_iter_binning(self.n_iter, self.binning.mapper_hash, self.binning.mapper_pickle,
                                           we_driver.bin_target_counts)
            data_manager.update_initial_states(we_driver.all_initial_states)
            data_manager.flush_backing()
            
            
            
            
            

if __name__ == '__main__':
    WBinTool().main()
