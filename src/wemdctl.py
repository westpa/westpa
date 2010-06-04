from __future__ import division
import wemd
from wemd import Segment, WESimIter
from wemd.sim_managers import get_sim_manager
from wemd.util.wetool import WECmdLineMultiTool
from wemd.rc import EX_USAGE_ERROR, EX_ENVIRONMENT_ERROR, EX_STATE_ERROR
from copy import copy

import logging
log = logging.getLogger(__name__)

import os, sys
from itertools import izip

class WEMDCtlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDCtlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('init', 'initialize a new WE simulation',
                        self.cmd_init, True)
        cop.add_command('reinit', 'reinitialize a WE simulation',
                        self.cmd_reinit, True)
        cop.add_command('status', 'show simulation status',
                        self.cmd_status, True)
        cop.add_command('reweight', 'apply new bin weights to a WE simulation',
                        self.cmd_reweight, True)
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
            sys.stderr.write('a simulation configuration file is required\n')
            parser.print_help(sys.stderr)
            sys.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        from wemd.util.config_dict import ConfigDict
        self.sim_config_src = ConfigDict()
        try:
            self.sim_config_src.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            sys.stderr.write('cannot open simulation config file: %s\n' % e)
            sys.exit(EX_ENVIRONMENT_ERROR)
        
        self.load_sim_manager(load_sim_config = False)
        self.sim_manager.sim_init(self.sim_config, self.sim_config_src)
        self.sim_manager.save_sim_state()
        sys.exit(0)

    def cmd_reinit(self, args):
        parser = self.make_parser('NEW_SIM_CONFIG', 
                                  'Start a new WEMD simulation from an old one according '
                                  +'to the directives in NEW_SIM_CONFIG')
        parser = self.make_parser('NEW_RUN_CONFIG', 
                                  'new_run.cfg that will be used for the reinit\'d simulation ')        
        (opts, args) = parser.parse_args(args)
        if len(args) != 2:
            sys.stderr.write('a simulation and new runtime configuration file is required\n')
            parser.print_help(sys.stderr)
            sys.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
            new_rc_file = args[1]
            
        from wemd.util.config_dict import ConfigDict
        self.sim_config_src = ConfigDict()
        try:
            self.sim_config_src.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            sys.stderr.write('cannot open simulation config file: %s\n' % e)
            sys.exit(EX_ENVIRONMENT_ERROR)

        try:
            new_runtime_config = wemd.rc.read_config(new_rc_file)
        except IOError, e:
            self.log.debug('cannot open new run config file', exc_info=True)
            sys.stderr.write('cannot open new run config file: %s\n' % e)
            sys.exit(EX_ENVIRONMENT_ERROR)            

        if self.runtime_config.get('data.state') == new_runtime_config.get('data.state')\
            or self.runtime_config.get('data.db.url') == new_runtime_config.get('data.db.url')\
            or self.runtime_config.get('data.sim_config') == new_runtime_config.get('data.sim_config'):
            self.log.debug('data parameters can\'t be the same \n see [data] in run config files\n', exc_info=True)
            sys.stderr.write('data parameters can\'t be the same \n see [data] in run config files\n')
            sys.exit(EX_ENVIRONMENT_ERROR)  

                
        sim_manager = self.load_sim_manager()
        sim_manager.load_sim_state()

        new_driver_name = new_runtime_config.get('sim_manager.driver', 'serial').lower()
        Manager = wemd.sim_managers.get_sim_manager(new_driver_name)   
        new_sim_manager = Manager()
        new_sim_manager.runtime_init(new_runtime_config, load_sim_config = False)
                    
        from wemd import Segment
        import numpy

        n_iter = sim_manager.we_driver.current_iteration
        
        new_segments = []
        old_segments = sim_manager.data_manager.get_segments(Segment.n_iter == n_iter, result_format='objects')
        for old_s in old_segments:
            assert old_s.status == Segment.SEG_STATUS_COMPLETE
            last_pcoord = numpy.array( [old_s.pcoord[-1][:]], dtype=numpy.float64)


            new_s = Segment( n_iter = 0,
                            status = Segment.SEG_STATUS_COMPLETE,
                            weight = old_s.weight,
                            pcoord = copy(last_pcoord),
                            data = {'old_we_iter': old_s.n_iter, 'old_seg_id': old_s.seg_id} )
            
            new_segments.append( new_s )

        new_sim_manager.sim_init(self.sim_config, self.sim_config_src, new_segments)
        new_sim_manager.save_sim_state() 
                
        sys.exit(0)
        
    def cmd_reweight(self, args):
        import cPickle as pickle
        import numpy
        parser = self.make_parser('SIM_CONFIG', 
                                  'apply new bin weights to a WE simulation '
                                  +' according to the directives in SIM_CONFIG')
        (opts, args) = parser.parse_args(args)

        if len(args) != 1:
            sys.stderr.write('a simulation configuration file is required\n')
            parser.print_help(sys.stderr)
            sys.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        from wemd.util.config_dict import ConfigDict
        sim_config_src = ConfigDict()
        try:
            sim_config_src.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            sys.stderr.write('cannot open simulation config file: %s\n' 
                                    % e)
            sys.exit(EX_ENVIRONMENT_ERROR)
            
        sim_config_src.require('bins.weights_from')
        
        new_weights = pickle.load(open(sim_config_src['bins.weights_from'], 'rb'))
                    
        from sqlalchemy.sql import select, delete
        from wemd.data_manager import schema
    
        sim_manager = self.load_sim_manager()
        sim_manager.load_sim_state()
        
        dbsession = sim_manager.data_manager.require_dbsession()
        n_iter = sim_manager.we_driver.current_iteration
        qsegs = dbsession.query(Segment).filter(Segment.n_iter == n_iter)
        
        dbsession.begin()
        n_total = qsegs.count()
        n_prepared = qsegs.filter(Segment.status == Segment.SEG_STATUS_PREPARED).count()
        
        try:
            if n_prepared and n_prepared != n_total:
                sys.stderr.write('%d segments are incomplete\n' 
                                        % n_prepared)
                sys.stderr.write('complete or truncate current iteration'
                                        +' before reweighting\n')
                sys.exit(EX_STATE_ERROR)
            
            dbsession.execute(delete(schema.segment_lineage_table,
                                     schema.segment_lineage_table.c.seg_id.in_(
                                         select([schema.segments_table.c.seg_id],
                                                schema.segments_table.c.n_iter == n_iter))))
            qsegs.delete()
            
            # set new bin parameters
            sim_manager.we_driver.sim_init(sim_manager.sim_config, sim_config_src)            
            sim_manager.we_driver.current_iteration = n_iter-1
            sim_manager.we_iter = sim_manager.data_manager.get_we_sim_iter(n_iter-1)
            dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).delete()
            dbsession.flush()
            sim_manager.run_we(reweight = new_weights)
            dbsession.flush()     
            sim_manager.save_sim_state()
            sim_manager.save_sim_config()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
            sys.exit(0)
        
        
        
    def cmd_rebin(self, args):
        parser = self.make_parser('SIM_CONFIG', 
                                  'apply new bin limits to a WE simulation '
                                  +' according to the directives in SIM_CONFIG')
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            sys.stderr.write('a simulation configuration file is required\n')
            parser.print_help(sys.stderr)
            sys.exit(EX_USAGE_ERROR)
        else:
            sim_config_file = args[0]
        
        from wemd.util.config_dict import ConfigDict
        sim_config_src = ConfigDict()
        try:
            sim_config_src.read_config_file(sim_config_file)
        except IOError, e:
            self.log.debug('cannot open simulation config file', exc_info=True)
            sys.stderr.write('cannot open simulation config file: %s\n' 
                                    % e)
            sys.exit(EX_ENVIRONMENT_ERROR)
            
        from sqlalchemy.sql import select, delete
        from wemd.data_manager import schema
        
        sim_manager = self.load_sim_manager()
        self.sim_manager.load_sim_state()
        
        dbsession = sim_manager.data_manager.require_dbsession()
        n_iter = sim_manager.we_driver.current_iteration
        qsegs = dbsession.query(Segment).filter(Segment.n_iter == n_iter)
        
        dbsession.begin()
        n_total = qsegs.count()
        n_prepared = qsegs.filter(Segment.status == Segment.SEG_STATUS_PREPARED).count()
        
        try:
            if n_prepared and n_prepared != n_total:
                sys.stderr.write('%d segments are incomplete\n' 
                                        % n_prepared)
                sys.stderr.write('complete or truncate current iteration'
                                        +' before rebinning\n')
                sys.exit(EX_STATE_ERROR)
            
            dbsession.execute(delete(schema.segment_lineage_table,
                                     schema.segment_lineage_table.c.seg_id.in_(
                                         select([schema.segments_table.c.seg_id],
                                                schema.segments_table.c.n_iter == n_iter))))
            qsegs.delete()
            # flush the deletes
            dbsession.flush()
            
            sim_manager.we_driver.sim_init(sim_manager.sim_config, sim_config_src)
            sim_manager.we_driver.current_iteration = n_iter-1
            sim_manager.we_iter = sim_manager.data_manager.get_we_sim_iter(n_iter-1)
            dbsession.query(WESimIter).filter(WESimIter.n_iter == n_iter).delete()
            # flush the delete of the WESimIter
            dbsession.flush()
            
            sim_manager.run_we()
            # flush the new WESimIter
            dbsession.flush()
            sim_manager.save_sim_state()
            sim_manager.save_sim_config()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
            sys.exit(0)
    
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
        
        sim_manager = self.load_sim_manager()
        sim_manager.load_sim_state()
        dbsession = sim_manager.data_manager.require_dbsession()
        
        from sqlalchemy.sql.functions import max as max_
        from sqlalchemy.sql import select, delete
        from wemd.data_manager.schema import segments_table, segment_lineage_table
        max_we_iter = dbsession.query(max_(WESimIter.n_iter)).scalar()
        
        try:
            n_iter = int(args[0])
        except IndexError:
            n_iter = sim_manager.we_driver.current_iteration
        except ValueError:
            sys.stderr.write('invalid iteration number %r\n' % args[0])
            sys.exit(EX_USAGE_ERROR)
        
        if opts.trunc_after:
            n_iter += 1
            
        if n_iter > max_we_iter:
            sys.stderr.write(('iteration number (%d) exceeds those found '
                              +'in the simulation (max = %d)')
                              % (n_iter, max_we_iter))
            sys.exit(EX_USAGE_ERROR)
            
        dbsession.begin()
        try:
            # Clean up lineage table
            dbsession.execute(delete(segment_lineage_table,
                                     segment_lineage_table.c.seg_id.in_(
                                         select([segments_table.c.seg_id],
                                                segments_table.c.n_iter >= n_iter))))
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
            sim_manager.save_sim_state()
        except:
            dbsession.rollback()
            raise
        else:
            dbsession.commit()
            sys.stdout.write('simulation truncated after iteration %d\n' % n_iter)
            self.cmd_status([])
            sys.exit(0)
    
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
                        sys.stderr.write('could not save history file: %s\n'
                                                % e)
        
        import wemd, sqlalchemy, numpy
        sim_manager = self.load_sim_manager()
        
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
        sys.exit(0)        
        
    def cmd_status(self, args):
        parser = self.make_parser('report the status of a WEMD simulation')
        parser.add_option('-i', '--iteration', dest='we_iter', type='int',
                          help='report status of iteration number WE_ITER '
                              +'(default: current iteration)')
        (opts, args) = parser.parse_args(args)

        import numpy
        from sqlalchemy.sql import desc
        sim_manager = self.load_sim_manager()
        dbsession = sim_manager.data_manager.require_dbsession()
        current_iter = dbsession.query(WESimIter)\
                                .order_by(desc(WESimIter.n_iter))\
                                .first()
        if opts.we_iter is None:
            sim_iter = current_iter
            sys.stdout.write('Current WE iteration: %d\n' 
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
        
        sys.stdout.write('Data source: %s\n' % self.runtime_config['data.db.url'])
        sys.stdout.write('  Schema version: %s\n' % sim_manager.data_manager.get_schema_version())
        sys.stdout.write('Total of %d segments in iteration %d\n'
                                 % (sim_iter.n_particles, sim_iter.n_iter))
        sys.stdout.write('population norm = %.8f\n' % sim_iter.norm)
        sys.stdout.write('error in norm = %.8g (machine epsilon = %.16g)\n'
                                 % (sim_iter.norm - 1,
                                    numpy.finfo(numpy.double).eps))
        sys.stdout.write('  %d segments pending (%.0f%%)\n'
                                 % (n_pending, pct_pending))
        sys.stdout.write('  %d segments complete (%.0f%%)\n'
                                 % (n_complete, pct_complete))
        
        
        if sim_iter.data and sim_iter.data.get('bins_population') is not None:
            bpop = sim_iter.data['bins_population']
            bpart = sim_iter.data['bins_nparticles']
            n_bins = bpop.size
            n_populated = bpart[bpart != 0].size
            sys.stdout.write('%d bins in %d dimension(s)\n'
                                     % (bpop.size, len(bpop.shape)))
            sys.stdout.write('  %d populated (%.0f%%)\n'
                                     % (n_populated, n_populated/n_bins*100))
        
    def cmd_update_schema(self, args):
        parser = self.make_parser('update the database to the most recent schema')
        (opts, args) = parser.parse_args(args)
        
        from wemd.data_manager.versioning import update_schema
        sim_manager = self.load_sim_manager()
        update_schema(sim_manager.data_manager.dbengine)
        
        

        
        
                
if __name__ == '__main__':
    WEMDCtlTool().run()

