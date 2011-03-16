from __future__ import division, print_function
import os
import math
import numpy

import logging
log = logging.getLogger('w_sqlite2hdf5')

import wemd
from wemd import Segment

parser = wemd.rc.common_arg_parser(description = '''\
Port a pre-0.5 SQLite DB to the current HDF5 format
''')
parser.add_argument('sqlite_db', help='SQLite DB to read')
parser.add_argument('h5file', help='HDF5 file to write')
parser.add_argument('maxiter', help='Read up to iteration MAXITER', type=int)
args = parser.parse_args()

from wemd.util.config_dict import ConfigDict
wemd.rc.config_logging(args, 'w_sqlite2hdf5')
#runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config = ConfigDict()
runtime_config.update_from_object(args)
sim_manager = wemd.rc.load_sim_manager(runtime_config)

# Create the HDF5 file
runtime_config['data.h5file'] = args.h5file
if os.path.exists(args.h5file): os.unlink(args.h5file)
sim_manager.load_data_manager()
sim_manager.data_manager.prepare_backing()

# Open the SQLite file
import wemdtools.convert.sqlite2hdf5, wemdtools.convert.sqlite2hdf5.sqlcore 
from wemdtools.convert.sqlite2hdf5 import sqlcore
from wemdtools.convert.sqlite2hdf5.sqlcore import SQLAlchemyDataManager
runtime_config['data.db.url'] = 'sqlite:///{}'.format(args.sqlite_db)
sql_dm = SQLAlchemyDataManager(sim_manager)

# A mapping of SQL seg_id -> (n_iter, hdf5_seg_id)  
idmap = dict()
for n_iter in xrange(1,args.maxiter+1):
    # get all old segments
    old_segments = sql_dm.get_segments(sqlcore.OldSegment.n_iter == n_iter, load_parents = True, load_p_parent = True)
    print('{} segments in iteration {}'.format(len(old_segments), n_iter))
    
    new_segments = []
    pcoord_ndim = old_segments[0].pcoord.shape[1]
    pcoord_len = old_segments[0].pcoord.shape[0]
    pcoord_dtype = old_segments[0].pcoord.dtype
    
    we_iter = sql_dm.get_we_sim_iter(n_iter)
    
    total_weight = 0.0
    min_prob = 1.0
    max_prob = 0.0
    
    for (seg_id, old_segment) in enumerate(old_segments):
        idmap[old_segment.seg_id] = (n_iter, seg_id)
        
        total_weight += old_segment.weight
        min_prob = min(min_prob or 1.0, old_segment.weight)
        max_prob = max(max_prob, old_segment.weight)

        if not old_segment.p_parent:
            p_parent_id = -1
            parent_ids = [-1]
        else:
            (prev_iter, p_parent_id) = idmap[old_segment.p_parent.seg_id]
            parent_ids = set(idmap[parent.seg_id][1] for parent in old_segment.parents)
            
        pcoord = old_segment.pcoord
        if pcoord.ndim < 2:
            pcoord = numpy.expand_dims(pcoord, 1)
        
        new_segment = Segment(seg_id = seg_id, n_iter = n_iter,
                              status = old_segment.status,
                              endpoint_type = old_segment.endpoint_type,
                              walltime = old_segment.walltime,
                              cputime = old_segment.cputime,
                              weight = old_segment.weight,
                              pcoord = pcoord,
                              p_parent_id = p_parent_id,
                              parent_ids = parent_ids,
                              )
        new_segments.append(new_segment)
        
    sim_manager.data_manager.prepare_iteration(n_iter, new_segments, pcoord_ndim, pcoord_len, pcoord_dtype)
    iter_summary = sim_manager.data_manager.get_iter_summary(n_iter)
    iter_summary['norm'] = total_weight
    iter_summary['n_particles'] = len(new_segments)
    iter_summary['cputime'] = we_iter.cputime
    iter_summary['walltime'] = we_iter.walltime
    iter_summary['min_seg_prob'] = min_prob    
    iter_summary['max_seg_prob'] = max_prob
    iter_summary['seg_dyn_range'] = math.log(max_prob/min_prob)
    
    populations = we_iter.data['bins_population']
    min_bin_prob = populations[populations!=0].min()
    max_bin_prob = populations.max()
    iter_summary['min_bin_prob'] = min_bin_prob
    iter_summary['max_bin_prob'] = max_bin_prob
    iter_summary['bin_dyn_range'] = math.log(max_bin_prob/min_bin_prob)
    iter_summary['target_flux']  = we_iter.data['recycled_population'] or 0.0
    sim_manager.data_manager.update_iter_summary(n_iter, iter_summary)
    
    # Delete everything we were working with, since the garbage collector seems
    # pretty bad at discerning when we're done with nested OldSegments.
    del iter_summary, populations, pcoord, new_segment, new_segments, old_segment, old_segments
    
sim_manager.data_manager.current_iteration = n_iter
