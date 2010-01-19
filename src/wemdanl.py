import os, sys
from itertools import izip
from optparse import OptionParser
import wemd
from wemd import Segment, WESimIter

from wemd.util.wetool import WECmdLineMultiTool
from wemd.environment import *
from wemd.sim_managers import make_sim_manager

from logging import getLogger
log = getLogger(__name__)

class WEMDAnlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDAnlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('mapsim',
                        'trace trajectories',
                        self.cmd_mapsim, True)
        self.add_rc_option(cop)
        
    def add_iter_param(self, parser):
        parser.add_option('-i', '--iteration', dest='we_iter', type='int',
                          help = 'use trajectories as of iteration number '
                               + 'WE_ITER (default: the last complete '
                               + ' iteration)'
                         )
        
    def get_sim_manager(self):
        self.sim_manager = make_sim_manager(self.runtime_config)
        return self.sim_manager
        
    def get_sim_iter(self, we_iter):
        if we_iter is None:
            dbsession = self.sim_manager.data_manager.require_dbsession()
            we_iter = dbsession.execute(
'''SELECT MAX(n_iter) FROM we_iter 
     WHERE NOT EXISTS (SELECT * FROM segments 
                         WHERE segments.n_iter = we_iter.n_iter 
                           AND segments.status != :status)''',
                           params = {'status': Segment.SEG_STATUS_COMPLETE}).fetchone()[0]
        
        self.we_iter = dbsession.query(WESimIter).get([we_iter])
        return self.we_iter
    
    def cmd_mapsim(self, args):
        # Walk tree and convert trajectory data to HDF5 array
        parser = self.make_parser()
        self.add_iter_param(parser)
        parser.add_option('--squeeze', dest='squeeze', action='store_true',
                          help='eliminate duplicate records at segment '
                              +'boundaries (default)')
        parser.add_option('--nosqueeze', dest='squeeze', action='store_false',
                          help='do not treat segment boundaries in any special '
                              +'way')
        parser.add_option('-o', '--output', dest='output_file',
                          help='store compiled output in OUTPUT_FILE '
                              +'(default: trajectories.h5)')
        parser.add_option('-t', '--timestep', dest='timestep', type='float',
                          help='manually specify dynamics timestep')
        parser.add_option('--timestep-units', dest='timestep_units',
                          help='store given timestep units in HDF5 file')
        parser.add_option('-T', '--segment-length', '--tau', dest='tau',
                          type='float',
                          help='manually specify segment length')
        parser.add_option('--segment-length-units', '--tau-units', 
                          dest='tau_units',
                          help='store given segment length units in HDF5 file')
        parser.set_defaults(squeeze = True,
                            output_file = 'trajectories.h5')
        (opts, args) = parser.parse_args(args)
        
        self.get_sim_manager()
        we_iter = self.get_sim_iter(opts.we_iter)
        
        import numpy, h5py, sqlalchemy
        from sqlalchemy.sql import select, bindparam
        schema = self.sim_manager.data_manager.get_schema()
        dbsession = self.sim_manager.data_manager.require_dbsession()
        n_trajs = dbsession.query(Segment).filter(Segment.n_iter == we_iter.n_iter).count()
        model_segment = dbsession.query(Segment).filter(Segment.n_iter == we_iter.n_iter).first()
        
        self.output_stream.write('%d live trajectories in iteration %d\n'
                                 % (n_trajs, we_iter.n_iter))
        self.output_stream.write('progress coordinate is %s for each segment\n'
                                 % ('x'.join(str(x) for x in model_segment.pcoord.shape),))
        
        n_md_steps = model_segment.pcoord.shape[0]
        ndim_pcoord = model_segment.pcoord.ndim - 1 
        try:
            dt = model_segment.data['timestep']
        except KeyError:
            dt = opts.timestep
            if dt is not None:
                self.output_stream.write('timestep %g taken from command line\n'
                                         % dt)
        else:
            self.output_stream.write('timestep %g read from database\n' % dt)

        h5file = h5py.File(opts.output_file)
        dynamics_group = h5file.create_group('Dynamics')
        we_group = h5file.create_group('WE')
                    
        seg_paths = dynamics_group.create_dataset('SegPaths', 
                                                  shape=(n_trajs,
                                                         we_iter.n_iter),
                                                  dtype=numpy.uint64)
        if opts.squeeze:
            # Omit first data point from all segments except the first
            pcoords_shape = ( n_trajs, 
                              (n_md_steps-1)*(we_iter.n_iter-1) + n_md_steps,) \
                            + model_segment.pcoord.shape[1:]
        else:
            pcoords_shape = ( n_trajs, n_md_steps * we_iter.n_iter,) \
                            + model_segment.pcoord.shape[1:]
        # pcoords[trajectory][t][...]
        pcoords = dynamics_group.create_dataset('ProgCoord',
                                                shape=pcoords_shape,
                                                dtype=model_segment.pcoord.dtype)
        
        
        # Store some metadata
        if dt is not None:
            dynamics_group.attrs['timestep'] = pcoords.attrs['timestep'] = dt
        else:
            log.warning('no timestep specified')
            
        if opts.timestep_units:
            dynamics_group.attrs['timestep'] = pcoords.attrs['timestep_units'] = opts.timestep_units
            
        if opts.tau:
            we_group.attrs['tau'] = opts.tau
        if opts.tau_units:
            we_group.attrs['tau_units'] = opts.tau_units
            
        # Do some fast select magic without reconstituting true Segment objects
        sel = select([schema.segmentsTable.c.seg_id,
                      schema.segmentsTable.c.p_parent_id,
                      schema.segmentsTable.c.pcoord,
                      schema.segmentsTable.c.data],
                     schema.segmentsTable.c.n_iter == bindparam('n_iter'))
        
        for n_iter in xrange(we_iter.n_iter, 0, -1):
            i_iter = n_iter - 1
                        
            rows = dbsession.execute(sel, params={'n_iter': n_iter}).fetchall()
            rows_by_seg_id = dict((row.seg_id, row) for row in rows)
            
            if n_iter == we_iter.n_iter:
                # If we ever set up automatic handling of everything in the 
                # "data" field, here's where we'd do it
                
                # Set up the last point to start up the induction
                seg_paths[:, i_iter] = rows_by_seg_id.keys()
                
            seg_ids = seg_paths[:, i_iter][...]
                
            if i_iter > 0:
                seg_paths[:, i_iter-1] = [rows_by_seg_id[seg_id].p_parent_id 
                                          for seg_id in seg_ids]
                if opts.squeeze:
                    seg_pc_lb = 1
                    pc_lb = i_iter * (n_md_steps - 1) + 1
                    pc_ub = pc_lb + (n_md_steps-1)
                else:
                    seg_pc_lb = 0
                    pc_lb = i_iter * n_md_steps
                    pc_ub = pc_lb + n_md_steps
            else:
                seg_pc_lb = 0
                pc_lb = 0
                pc_ub = n_md_steps
            
            pcoords[:, pc_lb:pc_ub] = \
                [rows_by_seg_id[seg_id].pcoord[seg_pc_lb:]
                 for seg_id in seg_ids
                ]
                
        h5file.close()
            
        
if __name__ == '__main__':
    WEMDAnlTool().run()

