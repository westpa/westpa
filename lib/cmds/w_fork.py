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



import os, sys, logging, argparse
import numpy
 
log = logging.getLogger('w_fork')

import westpa
from west import Segment
from west.states import InitialState
from west.data_manager import n_iter_dtype, seg_id_dtype

parser = argparse.ArgumentParser('w_fork',description='''\
Prepare a new weighted ensemble simulation from an existing one at a particular
point. A new HDF5 file is generated. In the case of executable propagation,
it is the user's responsibility to prepare the new simulation directory
appropriately, particularly making the old simulation's restart data from the
appropriate iteration available as the new simulations initial state data; a
mapping of old simulation segment to new simulation initial states is
created, both in the new HDF5 file and as a flat text file, to aid in this.
Target states and basis states for the new simulation are taken from those 
in the original simulation.  
''')

westpa.rc.add_args(parser)
parser.add_argument('-i', '--input', dest='input_h5file',
                    help='''Create simulation from the given INPUT_H5FILE (default: read from
                            configuration file.''')
parser.add_argument('-I', '--iteration', dest='n_iter', type=int,
                    help='''Take initial distribution for new simulation from iteration N_ITER
                            (default: last complete iteration).''')
parser.add_argument('-o', '--output', dest='output_h5file', default='forked.h5',
                    help='''Save new simulation HDF5 file as OUTPUT (default: %(default)s).''')
parser.add_argument('--istate-map', default='istate_map.txt',
                    help='''Write text file describing mapping of existing segments to new initial
                            states in ISTATE_MAP (default: %(default)s).''')
parser.add_argument('--no-headers', action='store_true',
                    help='''Do not write header to ISTATE_MAP''')
args = parser.parse_args()
westpa.rc.process_args(args)

# Open old HDF5 file
dm_old = westpa.rc.new_data_manager()
if args.input_h5file:
    dm_old.we_h5filename = args.input_h5file
dm_old.open_backing(mode='r')

# Get iteration if necessary
n_iter = args.n_iter or dm_old.current_iteration-1

# Create and open new HDF5 file
dm_new = westpa.rc.new_data_manager()
dm_new.we_h5filename = args.output_h5file
dm_new.prepare_backing()
dm_new.open_backing()

# Copy target states
target_states = dm_old.get_target_states(n_iter)
dm_new.save_target_states(target_states, n_iter)

# Copy basis states
basis_states = dm_old.get_basis_states(n_iter)
dm_new.create_ibstate_group(basis_states, n_iter=1)


# Transform old segments into initial states and new segments
# We produce one initial state and one corresponding
# new segment for each old segment. Further adjustment
# can be accomplished by using w_binning.
old_iter_group = dm_old.get_iter_group(n_iter)
old_index = old_iter_group['seg_index'][...]
old_pcoord_ds = old_iter_group['pcoord']
n_segments = old_pcoord_ds.shape[0]
pcoord_len = old_pcoord_ds.shape[1]
pcoord_ndim = old_pcoord_ds.shape[2]
old_final_pcoords = old_pcoord_ds[:,pcoord_len-1,:]

istates = dm_new.create_initial_states(n_segments, n_iter=1)
segments = []
state_map_dtype = numpy.dtype([('old_n_iter', n_iter_dtype),
                               ('old_seg_id', seg_id_dtype),
                               ('new_istate_id', seg_id_dtype)])
state_map = numpy.empty((n_segments,),dtype=state_map_dtype)
state_map['old_n_iter'] = n_iter


for (iseg, (index_row, pcoord)) in enumerate(zip(old_index, old_final_pcoords)):
    istate = istates[iseg]
    istate.iter_created = 0
    istate.iter_used = 1
    istate.istate_type = InitialState.ISTATE_TYPE_RESTART
    istate.istate_status = InitialState.ISTATE_STATUS_PREPARED
    istate.pcoord = pcoord
    
    segment = Segment(n_iter=1, seg_id=iseg, weight=index_row['weight'],
                      parent_id =-(istate.state_id+1),
                      wtg_parent_ids = [-(istate.state_id+1)], 
                      status=Segment.SEG_STATUS_PREPARED)
    segment.pcoord = numpy.zeros((pcoord_len, pcoord_ndim), dtype=pcoord.dtype)
    segment.pcoord[0] = pcoord
    segments.append(segment)
    state_map[iseg]['old_seg_id'] = iseg
    state_map[iseg]['new_istate_id'] = istate.state_id
    
dm_new.update_initial_states(istates, n_iter=0)
dm_new.prepare_iteration(n_iter=1, segments=segments)

# Update current iteration and close both files
dm_new.current_iteration = 1
dm_new.close_backing()
dm_old.close_backing()

# Write state map
istate_map_file = open(args.istate_map, 'wt')
if not args.no_headers:
    istate_map_file.write('# mapping from previous segment IDs to new initial states\n')
    istate_map_file.write('# generated by w_fork\n')
    istate_map_file.write('# column 0: old simulation n_iter\n')
    istate_map_file.write('# column 1: old simulation seg_id\n')
    istate_map_file.write('# column 2: new simulation initial state ID\n')

for row in state_map:
    istate_map_file.write('{old_n_iter:20d}    {old_seg_id:20d}    {new_istate_id:20d}\n'
                          .format(old_n_iter=int(row['old_n_iter']),
                                  old_seg_id=int(row['old_seg_id']),
                                  new_istate_id=int(row['new_istate_id'])))
    
    

