import os, sys
from itertools import izip
from optparse import OptionParser
import wemd
from wemd import Segment, WESimIter

from wemd.util.wetool import WECmdLineMultiTool
from wemd.environment import *
from wemd.sim_managers import make_sim_manager

class WEMDAnlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDAnlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('mapsim',
                        'trace trajectories',
                        self.cmd_mapsim, True)
        cop.add_command('simcat',
                        'show progress ofa  simulation',
                        self.cmd_simcat, True)
        cop.add_command('trajcat',
                        'show progress of a trajectory',
                        self.cmd_trajcat, True)
        cop.add_command('lstraj', 
                        'list and show information about trajectories',
                        self.cmd_lstraj, True)
        cop.add_command('genpc',
                        'extract progress coordinates for analysis',
                        self.cmd_genpc)
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
        live_segments = self.sim_manager.data_manager.get_segments(we_iter)
        model_segment = live_segments[0]
        n_md_steps = model_segment.pcoord.shape[0]
        ndim_pcoord = model_segment.pcoord.ndim - 1
        try:
            dt = model_segment.data['timestep']
        except KeyError:
            dt = opts.timestep

        self.output_stream.write('using WE iteration %d (%d live trajectories)\n' 
                                 % (we_iter.n_iter, len(live_segments)))
        self.output_stream.write('pcoord is %d-dimensional with %d frames per segment\n'
                                 % (ndim_pcoord, n_md_steps))
        
        import numpy, h5py, sqlalchemy
        from sqlalchemy.sql import select, bindparam
        schema = self.sim_manager.data_manager.get_schema()
        dbsession = self.sim_manager.data_manager.require_dbsession()

        h5file = h5py.File(opts.output_file)
        dynamics_group = h5file.create_group('Dynamics')
        we_group = h5file.create_group('WE')
                    
        seg_paths = dynamics_group.create_dataset('SegPaths', 
                                                  shape=(len(live_segments),
                                                         we_iter.n_iter),
                                                  dtype=numpy.uint64)
        if opts.squeeze:
            # Omit first data point from all segments except the first
            pcoords_shape = ( len(live_segments), (n_md_steps-1)*(we_iter.n_iter-1) + n_md_steps,) + model_segment.pcoord.shape[1:]
        else:
            pcoords_shape = ( len(live_segments), n_md_steps * we_iter.n_iter,) + model_segment.pcoord.shape[1:]
        # pcoords[trajectory][t][...]
        pcoords = dynamics_group.create_dataset('ProgCoord',
                                                shape=pcoords_shape,
                                                dtype=model_segment.pcoord.dtype)
        if dt is not None:
            pcoords.attrs['timestep'] = dt
        
        # Do some fast select magic without reconstituting true Segment objects
        sel = select([schema.segmentsTable.c.n_iter,
                      schema.segmentsTable.c.seg_id,
                      schema.segmentsTable.c.p_parent_id,
                      schema.segmentsTable.c.pcoord])
        
        n_iter = we_iter.n_iter
        segments = live_segments
        seg_ids = list(sorted([segment.seg_id for segment in live_segments]))
        segments_by_seg_id = dict((segment.seg_id, segment) for segment in live_segments)
        parent_ids_by_seg_id = dict((segment.seg_id, segment.p_parent.seg_id)
                                    for segment in segments)
        
        # We do not record WE iteration 0 (the bootstrap iteration), so
        # seg_paths[:, n] is iteration n+1
        if opts.squeeze:
            pcoord_row_offset = 1
            n_pcoord_rows = n_md_steps - 1
        else:
            pcoord_row_offset = 0
            n_pcoord_rows = n_md_steps
        
        seg_paths[:, n_iter-1] = seg_ids
        for (i_seg, seg_id) in enumerate(seg_ids):
            pc_lb = (n_iter-1)*n_pcoord_rows+pcoord_row_offset
            pc_ub = pc_lb + n_pcoord_rows
            pcoords[i_seg, pc_lb:pc_ub] = segments_by_seg_id[seg_id].pcoord[pcoord_row_offset:]
            
        seg_paths[:, n_iter-2] = [parent_ids_by_seg_id[seg_id] for seg_id in seg_ids]
        i_iter = n_iter - 2
        MAXCHUNK=768
        while i_iter >= 1:
            seg_ids = [long(seg_id) for seg_id in seg_paths[:, i_iter]]
            parent_ids = []
            for b_seg in xrange(0, len(seg_ids), MAXCHUNK):
                subset_ids = seg_ids[b_seg:b_seg+MAXCHUNK]
                sel2 = sel.where((schema.segmentsTable.c.n_iter == i_iter+1)
                                 & schema.segmentsTable.c.seg_id.in_(subset_ids))
                rows = dbsession.execute(sel2).fetchall()
                
                parent_ids_by_seg_id = dict((row.seg_id, row.p_parent_id) for row in rows)
                parent_ids.extend([parent_ids_by_seg_id[seg_id] for seg_id in subset_ids])
                                    
            seg_paths[:, i_iter-1] = parent_ids            
            i_iter -= 1
            
    def cmd_simcat(self, args):
        parser = self.make_parser()
        (opts, args) = parser.parse_args(args)
        
        self.get_sim_manager()
        sim_iters = self.sim_manager.dbsession.query(WESimIter).order_by(WESimIter.we_iter)
        last_iter = self.get_sim_iter(None)
        nbins = last_iter.data['bins_nparticles'].shape[0]
        
        npart_field = '%8d'
        pop_field = '%20.15g'
        
        rpt_fields = '%-8d    %-8d    %12.6g    %12.6g    %10g    %10g'
        hdr_text =   '%-8s    %-8s    %12s    %12s    %10s    %10s' \
                      % ('we_iter', 'n_part', 'norm', 'err_norm', 'cputime',
                         'walltime')
        self.output_stream.write(hdr_text + '\n')
        for sim_iter in sim_iters[:]:
            self.output_stream.write(rpt_fields % (sim_iter.we_iter,
                                                   sim_iter.n_particles,
                                                   sim_iter.norm,
                                                   sim_iter.norm - 1,
                                                   round(sim_iter.cputime or 0),
                                                   round(sim_iter.walltime or 0)))
            try:
                npart_text = '    '.join(npart_field % n for n in sim_iter.data['bins_nparticles'].flat)
                pop_text   = '    '.join(pop_field % p for p in sim_iter.data['bins_population'].flat)
            except TypeError:
                self.output_stream.write('\n')
            else:
                self.output_stream.write('    %s    %s\n' % (npart_text, pop_text))

    def cmd_lstraj(self, args):
        parser = self.make_parser()
        self.add_iter_param(parser)
        (opts, args) = parser.parse_args(args)
        
        self.get_sim_manager()
        self.get_sim_iter(opts.we_iter)
        
        col_fmt = {'seg_id': '%-8d',
                   'we_iter': '%-4d',
                   'data_ref': '%-24s',
                   'weight': '%-22.17e',
                   'cputime': '%-8.3f',
                   'walltime': '%-8.3f'}
        hdr_fmt = dict((k, v[:-1] + 's') for k,v in col_fmt.iteritems())
        cols = ('data_ref', 'seg_id', 'we_iter', 'weight', 
                'cputime', 'walltime')
        self.output_stream.write('#' + '    '.join(hdr_fmt[col] % col for col in cols) + '\n')
        for final_segment in self.sim_manager.dbsession.query(Segment)\
                                 .filter(Segment.we_iter == self.sim_iter.we_iter):
            self.output_stream.write('    '.join(col_fmt[col] 
                                              % getattr(final_segment, col) 
                                              for col in cols))
            self.output_stream.write('\n')
            
            
    def cmd_trajcat(self, args):
        parser = self.make_parser('SEG_ID',
                                  description = 
'''List and show information about trajectories. Emits information about the
trajectory ending with segment SEG_ID.''')
        (opts, args) = parser.parse_args(args)
        
        try:
            seg_id = int(args[0])
        except ValueError:
            self.error_stream.write('invalid segment ID %r specified\n' % args[0])
            self.exit(EX_USAGE_ERROR)
        except IndexError:
            self.error_stream.write('A segment ID must be provided\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
            
        self.get_sim_manager()
        
        segment = self.sim_manager.dbsession.query(Segment).get([seg_id])
        if segment is None:
            self.error_stream.write('no such segment (%r)\n' % seg_id)
            self.exit(EX_USAGE_ERROR)
        
        col_fmt = {'seg_id': '%-8d',
                   'we_iter': '%-4d',
                   'data_ref': '%-24s',
                   'weight': '%-22.17e',
                   'cputime': '%-8.3f',
                   'walltime': '%-8.3f'}
        hdr_fmt = dict((k, v[:-1] + 's') for k,v in col_fmt.iteritems())
        cols = ('data_ref', 'seg_id', 'we_iter', 'weight', 
                'cputime', 'walltime')

        
        segments = [segment]
        parent = segment.p_parent
        while parent and parent.data_ref:
            segments.append(parent)
            parent = parent.p_parent

        self.output_stream.write('#' + '    '.join(hdr_fmt[col] % col for col in cols) + '\n')
        for segment in reversed(segments):
            self.output_stream.write('    '.join(col_fmt[col] 
                                              % (getattr(segment, col) or 0) 
                                              for col in cols))
            self.output_stream.write('    %s' % ','.join(repr(pc) for pc in segment.pcoord))
            self.output_stream.write('\n')
                    
    def cmd_genpc(self, args):
        parser = self.make_parser(description = 
'''Extract and compile progress coordinate information for further analysis''')
        self.add_iter_param(parser)
        parser.add_option('-s', '--script', dest='script',
                          help='run SCRIPT on each trajectory segment to '
                              +'extract progress coordinate data')
        parser.add_option('--squeeze', dest='squeeze', action='store_true',
                          help='when obtaining data by script, '
                              + 'omit initial point from non-starting '
                              +'segment')
        parser.add_option('-o', '--output', dest='output_file',
                          default='wedist.h5',
                          help='store compiled output in OUTPUT_FILE '
                              +'(default: wedist.h5')
        parser.add_option('-t', '--timestep', dest='timestep', type='float',
                          help='timestep (MD or WE)')
        (opts, args) = parser.parse_args(args)
        
        self.get_sim_manager()
        self.get_sim_iter(opts.we_iter)
        
        import h5py, numpy
        
        # Get all living trajectories
        # Loop over them, one WE iteration at a time
        # Collect PC data and condense to binary format
        q = self.sim_manager.dbsession.query(Segment)\
                         .filter(Segment.we_iter == self.sim_iter.we_iter)
        ndim = len(q.first().pcoord.shape)
        n_we_iter = self.sim_iter.we_iter
        ntraj = q.count()
        dt = opts.timestep or 1
        
        
        h5file = h5py.File(opts.output_file, 'w')
        grp = h5file.create_group('pcoords')
        pcoords = grp.create_dataset('pcoord', 
                                     (n_we_iter, ntraj, ndim), numpy.float64,
                                      maxshape=(None, ntraj, ndim),
                                      chunks=(n_we_iter,1,1))
        weights = grp.create_dataset('weight', 
                                     (n_we_iter, ntraj), numpy.float64,
                                     maxshape=(None, ntraj),
                                     chunks=(n_we_iter,1))
        pcoords.attrs['timestep'] = dt
        weights.attrs['timestep'] = dt
        
        trajectories = []
        for (itraj, final_segment) in enumerate(q[:]):
            print "Tracing trajectory %d of %d" % (itraj+1, ntraj)
            segments = [final_segment]
            parent = final_segment.p_parent
            while parent and parent.data_ref:
                segments.append(parent)
                parent = parent.p_parent
            for (iseg, seg) in enumerate(reversed(segments)):
                pcoords[iseg, itraj, ...] = seg.pcoord
                weights[iseg, itraj] = seg.weight
                
        h5file.close()
        
if __name__ == '__main__':
    WEMDAnlTool().run()

