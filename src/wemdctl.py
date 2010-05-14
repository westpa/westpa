from __future__ import division
import wemd
from wemd import Segment, WESimIter
from wemd.sim_managers import make_sim_manager
from wemd.util.wetool import WECmdLineMultiTool
from wemd.rc import EX_USAGE_ERROR, EX_ENVIRONMENT_ERROR, EX_STATE_ERROR

import logging
log = logging.getLogger(__name__)

import os, sys

class WEMDCtlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDCtlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('init', 'initialize a new WE simulation',
                        self.cmd_init, True)
        cop.add_command('status', 'show simulation status',
                        self.cmd_status, True)
        cop.add_command('rebin', 'apply new bin limits to a WE simulation',
                        self.cmd_rebin, True)
        cop.add_command('truncate', 'truncate a WE simulation',
                        self.cmd_truncate, True)
        cop.add_command('console', 'open an interactive Python console',
                        self.cmd_console, True)
        cop.add_command('update_schema', 
                        'update a database to the most recent schema',
                        self.cmd_update_schema, True)
        self.add_rc_option(cop)
        
    def cmd_init(self, args):
        parser = self.make_parser('SIM_CONFIG', 
                                  'Configure a new WEMD simulation according '
                                  +'to the directives in SIM_CONFIG')
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            self.error_stream.write('a simulation configuration file is required\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        from wemd.util.config_dict import ConfigDict
        sim_config = ConfigDict()
        try:
            sim_config.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            self.error_stream.write('cannot open simulation config file: %s\n' % e)
            self.exit(EX_ENVIRONMENT_ERROR)
        
        sim_manager = make_sim_manager(self.runtime_config)
        sim_manager.initialize_simulation(sim_config)
        sim_manager.save_state()
        self.exit()
        
    def cmd_rebin(self, args):
        parser = self.make_parser('SIM_CONFIG', 
                                  'apply new bin limits to a WE simulation '
                                  +' according to the directives in SIM_CONFIG')
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            self.error_stream.write('a simulation configuration file is required\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        from wemd.util.config_dict import ConfigDict
        sim_config = ConfigDict()
        try:
            sim_config.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            self.error_stream.write('cannot open simulation config file: %s\n' 
                                    % e)
            self.exit(EX_ENVIRONMENT_ERROR)
            
        from sqlalchemy.sql import select, delete
        from wemd.data_manager import schema
            
        sim_manager = make_sim_manager(self.runtime_config)
        sim_manager.restore_state()
        
        dbsession = sim_manager.data_manager.require_dbsession()
        n_iter = sim_manager.we_driver.current_iteration
        qsegs = dbsession.query(Segment).filter(Segment.n_iter == n_iter)
        
        dbsession.begin()
        n_total = qsegs.count()
        n_prepared = qsegs.filter(Segment.status == Segment.SEG_STATUS_PREPARED).count()
        
        try:
            if n_prepared and n_prepared != n_total:
                self.error_stream.write('%d segments are incomplete\n' 
                                        % n_prepared)
                self.error_stream.write('complete or truncate current iteration'
                                        +' before rebinning\n')
                self.exit(EX_STATE_ERROR)
            
            dbsession.execute(delete(schema.segmentLineageTable,
                                     schema.segmentLineageTable.c.seg_id.in_(
                                         select([schema.segmentsTable.c.seg_id],
                                                schema.segmentsTable.c.n_iter == n_iter))))
            qsegs.delete()
            dbsession.flush()
            
            sim_manager.load_we_driver(sim_config)
            sim_manager.we_driver.current_iteration = n_iter-1
            sim_manager.we_iter = sim_manager.data_manager.get_we_sim_iter(n_iter-1)
            dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).delete()
            sim_manager.run_we()
            dbsession.flush()
            sim_manager.save_state()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
            self.exit(0)
    
    def cmd_truncate(self, args):
        parser = self.make_parser('[N_ITER]',
                                  'truncate the simulation record just prior to'
                                  + ' N_ITER '
                                  +'(thus, WE will run on the results of '
                                  +'iteration N_ITER-1). If no N_ITER is '
                                  +'given, then delete the current iteration.')
        parser.add_option('-A', '--after', 
                          dest='trunc_after', action='store_true',
                          help='truncate after N_ITER, instead of before')
        (opts, args) = parser.parse_args(args)
        
        sim_manager = make_sim_manager(self.runtime_config)
        sim_manager.restore_state()
        dbsession = sim_manager.data_manager.require_dbsession()
        
        from sqlalchemy.sql.functions import max as max_
        from sqlalchemy.sql import select, delete
        from wemd.data_manager.schema import segmentsTable, segmentLineageTable
        max_we_iter = dbsession.query(max_(WESimIter.n_iter)).scalar()
        
        try:
            n_iter = int(args[0])
        except IndexError:
            n_iter = sim_manager.we_driver.current_iteration
        except ValueError:
            sys.stderr.write('invalid iteration number %r\n' % args[0])
            self.exit(EX_USAGE_ERROR)
        
        if opts.trunc_after:
            n_iter += 1
            
        if n_iter > max_we_iter:
            sys.stderr.write(('iteration number (%d) exceeds those found '
                              +'in the simulation (max = %d)')
                              % (n_iter, max_we_iter))
            self.exit(EX_USAGE_ERROR)
            
        dbsession.begin()
        try:
            # Clean up lineage table
            dbsession.execute(delete(segmentLineageTable,
                                     segmentLineageTable.c.seg_id.in_(
                                         select([segmentsTable.c.seg_id],
                                                segmentsTable.c.n_iter >= n_iter))))
            # Delete segments
            dbsession.query(Segment).filter(Segment.n_iter >= n_iter).delete()
            
            # Delete iteration record
            dbsession.query(WESimIter).filter(WESimIter.n_iter >= n_iter).delete()
            
            # Update sim manager
            prev_iter = dbsession.query(WESimIter).get([n_iter-1])
            sim_manager.we_driver.current_iteration = prev_iter.n_iter
            sim_manager.we_driver.bins = sim_manager.we_driver.make_bins()
            sim_manager.we_driver.n_particles = prev_iter.data.get('bin_n_particles')
            sim_manager.we_driver.bins_population = prev_iter.data.get('bin_populations')
            sim_manager.save_state()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
            self.exit(0)
    
    def cmd_console(self, args):
        parser = self.make_parser('open an interactive Python console with '
                                  +'the current simulation loaded')
        (opts, args) = parser.parse_args(args)
        
        import code
        try:
            import readline, rlcompleter
        except Exception, e:
            class HistoryConsole(code.InteractiveConsole):
                def __init__(self, locals, filename='<console>',
                             histfile = os.path.expanduser('~/.wemd_history')):
                    code.InteractiveConsole.__init__(self, locals, filename)
        else:
            class HistoryConsole(code.InteractiveConsole):
                def __init__(self, locals, filename='<console>',
                             histfile = os.path.expanduser('~/.wemd_history')):
                    code.InteractiveConsole.__init__(self, locals, filename)
                    self.init_history(histfile)
                
                def init_history(self, histfile):
                    import atexit
                    readline.parse_and_bind('tab: complete')
                    try:
                        readline.read_history_file(histfile)
                    except IOError:
                        pass
    
                    atexit.register(self.save_history, histfile)
                
                def save_history(self, histfile):
                    try:
                        readline.write_history_file(histfile)
                    except Exception, e:
                        self.error_stream.write('could not save history file: %s\n'
                                                % e)
        
        import wemd, sqlalchemy, numpy
        sim_manager = make_sim_manager(self.runtime_config)
        
        banner = 'WEMD Python Console'
        
        locals = {'os': os, 'sys': sys, 
                  'sqlalchemy': sqlalchemy, 'numpy': numpy,
                  'wemd': wemd,
                  'Segment': wemd.Segment, 'Particle': wemd.Particle,
                  'WESimIter': wemd.WESimIter, 'Trajectory': wemd.Trajectory,
                  'sim_manager': sim_manager}
        
        try:
            import mpi4py
            from mpi4py import MPI
        except ImportError:
            pass
        else:
            locals['mpi4py'] = mpi4py
            locals['MPI'] = MPI
        
        try:
            readline.set_completer(rlcompleter.Completer(locals).complete)
        except NameError:
            pass
        
        hc = HistoryConsole(locals = locals)
        hc.interact(banner)
        self.exit(0)        
        
    def cmd_status(self, args):
        parser = self.make_parser('report the status of a WEMD simulation')
        parser.add_option('-i', '--iteration', dest='we_iter', type='int',
                          help='report status of iteration number WE_ITER '
                              +'(default: current iteration)')
        (opts, args) = parser.parse_args(args)

        import numpy
        from sqlalchemy.sql import desc
        sim_manager = make_sim_manager(self.runtime_config)
        dbsession = sim_manager.data_manager.require_dbsession()
        current_iter = dbsession.query(WESimIter)\
                                .order_by(desc(WESimIter.n_iter))\
                                .first()
        if opts.we_iter is None:
            sim_iter = current_iter
            self.output_stream.write('Current WE iteration: %d\n' 
                                     % current_iter.n_iter)
        else:
            sim_iter = dbsession.query(WESimIter).get([opts.we_iter])
            
        segq = dbsession.query(Segment)\
                        .filter(Segment.n_iter == sim_iter.n_iter)
        
        n_total = segq.count()
        n_pending = segq.filter(Segment.status == Segment.SEG_STATUS_PREPARED).count()
        n_complete = segq.filter(Segment.status == Segment.SEG_STATUS_COMPLETE).count()
        pct_pending = n_pending / n_total * 100
        pct_complete = n_complete / n_total * 100
        
        self.output_stream.write('Data source: %s\n' % self.runtime_config['data.db.url'])
        self.output_stream.write('  Schema version: %s\n' % sim_manager.data_manager.get_schema_version())
        self.output_stream.write('Total of %d segments in iteration %d\n'
                                 % (sim_iter.n_particles, sim_iter.n_iter))
        self.output_stream.write('population norm = %.8f\n' % sim_iter.norm)
        self.output_stream.write('error in norm = %.8g (machine epsilon = %.16g)\n'
                                 % (sim_iter.norm - 1,
                                    numpy.finfo(numpy.double).eps))
        self.output_stream.write('  %d segments pending (%.0f%%)\n'
                                 % (n_pending, pct_pending))
        self.output_stream.write('  %d segments complete (%.0f%%)\n'
                                 % (n_complete, pct_complete))
        
        
        if sim_iter.data and sim_iter.data.get('bins_population') is not None:
            bpop = sim_iter.data['bins_population']
            bpart = sim_iter.data['bins_nparticles']
            n_bins = bpop.size
            n_populated = bpart[bpart != 0].size
            self.output_stream.write('%d bins in %d dimension(s)\n'
                                     % (bpop.size, len(bpop.shape)))
            self.output_stream.write('  %d populated (%.0f%%)\n'
                                     % (n_populated, n_populated/n_bins*100))
        
    def cmd_update_schema(self, args):
        parser = self.make_parser('update the database to the most recent schema')
        (opts, args) = parser.parse_args(args)
        
        from wemd.data_manager.versioning import update_schema
        sim_manager = make_sim_manager(self.runtime_config)
        update_schema(sim_manager.data_manager.dbengine)
        
        

        
        
                
if __name__ == '__main__':
    WEMDCtlTool().run()

