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
                        'trace all trajectories into an HDF5 file',
                        self.cmd_mapsim, True)
        cop.add_command('transanl',
                        'analyze transitions',
                        self.cmd_transanl, True)
        cop.add_command('tracetraj',
                        'trace an individual trajectory and report',
                        self.cmd_tracetraj, True)
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
        dbsession = self.sim_manager.data_manager.require_dbsession()
        if we_iter is None:
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
            dt = model_segment.data['dt']
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
                    
        seg_paths_ds = dynamics_group.create_dataset('SegPaths', 
                                                  shape=(n_trajs,
                                                         we_iter.n_iter),
                                                  dtype=numpy.uint64)
        seg_paths = numpy.empty(shape=seg_paths_ds.shape, dtype=numpy.uint64)
        if opts.squeeze:
            # Omit first data point from all segments except the first
            pcoords_shape = ( n_trajs, 
                              (n_md_steps-1)*(we_iter.n_iter-1) + n_md_steps,) \
                            + model_segment.pcoord.shape[1:]    
        else:
            pcoords_shape = ( n_trajs, n_md_steps * we_iter.n_iter,) \
                            + model_segment.pcoord.shape[1:]

        # pcoords[trajectory][t][...]
        pcoords_ds = dynamics_group.create_dataset('ProgCoord',
                                                shape=pcoords_shape,
                                                dtype=model_segment.pcoord.dtype)
        pcoords = numpy.empty(shape=pcoords_shape, dtype=model_segment.pcoord.dtype)
        
        weights_shape = pcoords_shape[0:2]
        weights_ds = dynamics_group.create_dataset('Weight',
                                                shape=weights_shape,
                                                dtype=numpy.float64)
        weights = numpy.empty(shape=weights_shape, dtype=numpy.float64)
        
        # Store some metadata
        if dt is not None:
            dynamics_group.attrs['timestep'] = pcoords_ds.attrs['timestep'] = dt
        else:
            log.warning('no timestep specified')
            
        if opts.timestep_units:
            dynamics_group.attrs['timestep'] = pcoords_ds.attrs['timestep_units'] = opts.timestep_units
            
        if opts.tau:
            we_group.attrs['tau'] = opts.tau
        if opts.tau_units:
            we_group.attrs['tau_units'] = opts.tau_units
            
        # Do some fast select magic without reconstituting true Segment objects
        sel = select([schema.segmentsTable.c.seg_id,
                      schema.segmentsTable.c.p_parent_id,
                      schema.segmentsTable.c.weight,
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
                
            seg_ids = seg_paths[:, i_iter]
                
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
                
            for (i_seg, seg_id) in enumerate(seg_ids):
                row = rows_by_seg_id[seg_id]
                pcoords[i_seg, pc_lb:pc_ub] = row.pcoord[seg_pc_lb:]
                weights[i_seg, pc_lb:pc_ub] = row.weight
            
            #pcoords[:, pc_lb:pc_ub] = \
            #    [rows_by_seg_id[seg_id].pcoord[seg_pc_lb:]
            #     for seg_id in seg_ids
            #    ]
        
        seg_paths_ds[...] = seg_paths
        weights_ds[...] = weights
        pcoords_ds[...] = pcoords     
        h5file.close()
        
    def cmd_transanl(self, args):
        # Find transitions
        parser = self.make_parser()
        self.add_iter_param(parser)
        parser.add_option('-a', '--analysis-config', dest='anl_config',
                          help='use ANALYSIS_CONFIG as configuration file '
                              +'(default: analysis.cfg)')
        parser.set_defaults(anl_config = 'analysis.cfg')
        (opts,args) = parser.parse_args(args)
        
        from wemd.util.config_dict import ConfigDict, ConfigError
        transcfg = ConfigDict({'data.squeeze': True,
                               'data.timestep.units': 'ps'})
        transcfg.read_config_file(opts.anl_config)
        transcfg.require('regions.edges')
        transcfg.require('regions.names')
        
        region_names = transcfg.get_list('regions.names')
        region_edges = [float(v) for v in transcfg.get_list('regions.edges')]
        if len(region_edges) != len(region_names) + 1:
            self.error_stream.write('region names and boundaries do not match\n')
            self.exit(EX_ERROR)
            
        try:
            translog = open(transcfg['transitions.log'], 'wt')
        except KeyError:
            translog = None 
        
        regions = []
        for (irr,rname) in enumerate(region_names):
            regions.append((rname, (region_edges[irr], region_edges[irr+1])))
            
        from wemd.analysis.transitions import OneDimTransitionEventFinder
        
        self.get_sim_manager()
        final_we_iter = self.get_sim_iter(opts.we_iter)
        
        import numpy, h5py, sqlalchemy
        from sqlalchemy.sql import select, bindparam
        schema = self.sim_manager.data_manager.get_schema()
        dbsession = self.sim_manager.data_manager.require_dbsession()
        
        # Do some fast select magic without reconstituting true Segment objects
        sel = select([schema.segmentsTable.c.seg_id,
                      schema.segmentsTable.c.p_parent_id,
                      schema.segmentsTable.c.weight,
                      schema.segmentsTable.c.pcoord,
                      schema.segmentsTable.c.data],
                      schema.segmentsTable.c.n_iter == bindparam('n_iter'))
                
        n_md_steps = None
        ndim_pcoord = None
        dt = None
        
        event_durations = numpy.empty((0,2), numpy.float64)
        event_counts = numpy.zeros((len(regions), len(regions)), numpy.uint64)

        for max_iter in xrange(1, final_we_iter.n_iter+1):      
            n_trajs = dbsession.query(Segment).filter(Segment.n_iter == max_iter).count()
            model_segment = dbsession.query(Segment).filter(Segment.n_iter == max_iter).first()
            
            self.output_stream.write('%d live trajectories in iteration %d\n'
                                     % (n_trajs, max_iter))
        
            if n_md_steps is None:
                self.output_stream.write('progress coordinate is %s for each segment\n'
                                         % ('x'.join(str(x) for x in model_segment.pcoord.shape),))
                n_md_steps = model_segment.pcoord.shape[0]
                ndim_pcoord = model_segment.pcoord.ndim - 1 
                try:
                    dt = model_segment.data['dt']
                except KeyError:
                    dt = opts.timestep
                    if dt is not None:
                        self.output_stream.write('timestep %g taken from command line\n'
                                                 % dt)
                else:
                    self.output_stream.write('timestep %g read from database\n' % dt)

            seg_paths = numpy.empty(shape=(n_trajs, max_iter), 
                                    dtype=numpy.uint64)
                        
            if transcfg.get_bool('data.squeeze'):
                # Omit first data point from all segments except the first
                pcoords_shape = ( n_trajs, 
                                  (n_md_steps-1)*(max_iter-1) + n_md_steps,) \
                                + model_segment.pcoord.shape[1:]    
            else:
                pcoords_shape = ( n_trajs, n_md_steps * max_iter,) \
                                + model_segment.pcoord.shape[1:]
    
            # pcoords[trajectory][t][...]
            pcoords = numpy.empty(shape=pcoords_shape, dtype=model_segment.pcoord.dtype)            
            weights_shape = pcoords_shape[0:2]
            weights = numpy.empty(shape=weights_shape, dtype=numpy.float64)
                            
            for n_iter in xrange(max_iter, 0, -1):
                i_iter = n_iter - 1
                            
                rows = dbsession.execute(sel, params={'n_iter': n_iter}).fetchall()
                rows_by_seg_id = dict((row.seg_id, row) for row in rows)
                
                if n_iter == max_iter:
                    # Set up the last point to start up the induction
                    seg_paths[:, i_iter] = rows_by_seg_id.keys()
                    
                seg_ids = seg_paths[:, i_iter]
                    
                if i_iter > 0:
                    seg_paths[:, i_iter-1] = [rows_by_seg_id[seg_id].p_parent_id 
                                              for seg_id in seg_ids]
                    if transcfg.get_bool('data.squeeze'):
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
                    
                for (i_seg, seg_id) in enumerate(seg_ids):
                    row = rows_by_seg_id[seg_id]
                    pcoords[i_seg, pc_lb:pc_ub] = row.pcoord[seg_pc_lb:]
                    weights[i_seg, pc_lb:pc_ub] = row.weight
            
            for itraj in xrange(0, seg_paths.shape[0]):
                trans_finder = OneDimTransitionEventFinder(regions,
                                                           pcoords[itraj],
                                                           dt = dt,
                                                           traj_id=seg_paths[itraj,-1],
                                                           weights = weights[itraj],
                                                           transition_log = translog)
                trans_finder.identify_regions()
                trans_finder.identify_transitions()
                event_counts += trans_finder.event_counts

                try:
                    tfed = trans_finder.event_durations[2,0]
                except KeyError:
                    pass
                else:
                    event_durations.resize((event_durations.shape[0] + tfed.shape[0],2))
                    event_durations[-tfed.shape[0]:,:] = tfed[:,:]

        if event_durations.shape[0]:
            (ed_mean, ed_norm) = numpy.average(event_durations[:,0],
                                               weights = event_durations[:,1],
                                               returned = True)    
            self.output_stream.write('Number of events:    %d\n' % event_durations.shape[0])
            self.output_stream.write('ED average:          %g\n' % event_durations[:,0].mean())
            self.output_stream.write('ED weighted average: %g\n' % ed_mean)
            self.output_stream.write('ED min:              %g\n' % event_durations[:,0].min())
            self.output_stream.write('ED median:           %g\n' % event_durations[event_durations.shape[0]/2,0])
            self.output_stream.write('ED max:              %g\n' % event_durations[:,0].max())
        
            
    def cmd_tracetraj(self, args):
        parser = self.make_parser('TRAJ_ID')
        self.add_iter_param(parser)
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        
        self.get_sim_manager()
        we_iter = self.get_sim_iter(opts.we_iter)
        
        seg = self.sim_manager.data_manager.get_segment(we_iter, long(args[0]))
        segs = [seg]
        while seg.p_parent and seg.n_iter:
            seg = seg.p_parent
            segs.append(seg)
        traj = list(reversed(segs))
        
        for seg in traj:
            self.output_stream.write('%5d  %10d  %10d  %20.16g  %20.16g\n'
                                     % (seg.n_iter,
                                        seg.seg_id,
                                        seg.p_parent_id or 0,
                                        seg.weight,
                                        seg.pcoord[-1]))
        
if __name__ == '__main__':
    WEMDAnlTool().run()

