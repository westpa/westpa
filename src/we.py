#!/usr/bin/python

from __future__ import division

import os, sys, re, logging, datetime
import we

from optparse import OptionParser
from we.util.wetool import WECmdLineMultiTool
from we.core import Segment
from we.core.we_sim import PropagationIncompleteError
from we.environment import *

class WEUberTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEUberTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('init', 'initialize a new WE simulation', 
                        self.cmd_init, True)
        cop.add_command('showsetup', 'show simulation setup parameters', 
                        self.cmd_showsetup, True)
        cop.add_command('status', 'show simulation status', 
                        self.cmd_status, True)
        cop.add_command('cstatus', 'command-line status utility',
                        self.cmd_cstatus, True)
        cop.add_command('nextseg', 'prepare for the next MD segment', 
                        self.cmd_nextseg, True)
        cop.add_command('updateseg', 'record results for a given segment', 
                        self.cmd_updateseg, True)
        cop.add_command('we', 'perform WE analysis', 
                        self.cmd_we, True)
        cop.add_command('scrub', 'prune data from crashed segments',
                        self.cmd_scrub, True)
        
    def cmd_init(self, args):
        parser = self.make_parser('[options] DB CONFIG')
        parser.add_option('-f', '--force', dest='force', action='store_true',
                          help='overwrite existing state')
        
        (opts, args) = parser.parse_args(args)
        if len(args) != 2:
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        else:
            (dbfile, configfile) = args
            
        from we.util.config_dict import ConfigDict
        from we.data_managers import make_data_manager
        from we.sim_drivers import make_sim_driver
        
        self.output_stream.write('reading simulation configuration file %s\n' 
                                 % configfile)
        config = ConfigDict()
        config.read_config_file(configfile)
        
        if os.path.exists(dbfile):
            if opts.force:
                os.unlink(dbfile)
            else:
                self.error_stream.write("will not overwrite existing state "
                                        + "file; use -f to force overwrite\n")
                self.exit(EX_ENVIRONMENT_ERROR)
        
        self.output_stream.write('creating state file %s\n' % dbfile)
        datamgr = self.datamgr = make_data_manager(dbfile)
        datamgr.configure(config)
        datamgr.initialize()
        datamgr.we_sim = we_sim = make_sim_driver(datamgr, config)
        we_sim.configure(config)
        we_sim.initialize()
        datamgr.save_state()        
        self.showsetup(opts)
        
    def showsetup(self, opts = None):
        try:
            verbose = opts.verbose
        except AttributeError:
            verbose = False

        
        self.output_stream.write('\nBin space:\n')
        for (idim, limits) in enumerate(self.datamgr.config['bins.boundaries']):
            self.output_stream.write('  Dimension %d: %d bins\n' % (idim, len(limits)-1))
            binwidth = len(str(len(limits)-1))
            for ibound in xrange(0, len(limits)-1):
                self.output_stream.write('    Bin %*d: [%g, %g)\n' 
                                         % (binwidth, ibound, \
                                            limits[ibound], limits[ibound+1]))
        
        if verbose:        
            self.output_stream.write('\nConfiguration dump:\n')
            config = self.datamgr.config
            maxlen = max(len(k) for k in config)
            for key in sorted(config):
                self.output_stream.write('  %-*s : %s\n' 
                                         % (maxlen, key, config[key]))
                        
    def cmd_showsetup(self, args):
        parser = self.make_parser()
        parser.add_option('-v', '--verbose', dest='verbose', action='store_true')
        self.add_state_option(parser)
                
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        self.showsetup(opts)
        
    def cmd_status(self, args):
        import numpy
        from we.util import datetime_from_iso
        parser = self.make_parser()
        self.add_state_option(parser)
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        datamgr = self.datamgr
        current_iteration = datamgr.we_sim.current_iteration
        
        curs = datamgr.conn.cursor()
        # Start a transaction so we get a coherent view of the entire DB
        curs.execute('''BEGIN IMMEDIATE TRANSACTION''')
        curs.execute('''SELECT status, COUNT(*) FROM segments
                        WHERE we_iter=? GROUP BY status''',
                     (current_iteration,))
        count_by_status = dict((row[0], row[1]) for row in curs)
        
        curs.execute('''SELECT MIN(starttime), MAX(endtime) FROM segments
                        WHERE we_iter=? AND status=?
                        AND starttime IS NOT NULL AND endtime IS NOT NULL''', 
                        (current_iteration, Segment.SEG_STATUS_COMPLETE))
        row = curs.fetchone()

        try:
            we_start, we_end = datetime_from_iso(row[0]), datetime_from_iso(row[1])
            iter_wallclock = we_end - we_start
        except TypeError:
            we_start = we_end = iter_wallclock = 0
        
        curs.execute('''SELECT cputime FROM segments 
                        WHERE we_iter=? AND cputime IS NOT NULL''',
                     (current_iteration,))
        cputimes = numpy.fromiter((row[0] for row in curs), numpy.single)
        
        curs.execute('''SELECT weight FROM segments WHERE we_iter=?''',
                     (current_iteration,))
        weights = numpy.fromiter((row[0] for row in curs), numpy.float_)
        nparticles = len(weights)
        
        curs.execute('''SELECT COUNT(*), SUM(cputime) FROM segments 
                        WHERE status=?''',
                     (Segment.SEG_STATUS_COMPLETE,))
        row=curs.fetchone()
        sim_segs, sim_cpu = row[0], row[1]
        
        curs.execute('''SELECT we_iter, 
                               STRFTIME('%s',MIN(starttime)),
                               STRFTIME('%s',MAX(endtime)), COUNT(*)
                        FROM segments WHERE status=?
                          AND starttime IS NOT NULL AND endtime IS NOT NULL
                        GROUP BY we_iter''',
                     (Segment.SEG_STATUS_COMPLETE,))
        throughputs = []
        for row in curs:
            throughputs.append(3600 / (long(row[2])-long(row[1])) * row[3])
        throughputs = numpy.array(throughputs)
                
        curs.execute('''COMMIT''')
                
        bins_population = datamgr.we_sim.bins_population
        bins_nparticles = datamgr.we_sim.bins_nparticles
        
        bins_populated = bins_nparticles[bins_nparticles>0].size
        nbins = bins_nparticles.size
                    
        self.output_stream.write("Current WE iteration: %d\n" 
                                 % current_iteration)
        self.output_stream.write('''\
  Number of segments: %d
    Prepared:         %d
    Running:          %d
    Completed:        %d
  Populated bins:     %d / %d (%.0f%%)
  Norm:               %.8f
  Cumulative CPU (s): %.1f
  Average CPU:        %.1f
  Min CPU:            %.1f
  Max CPU:            %.1f
  Wallclock time:     %s
''' % (len(weights), 
       count_by_status.get(Segment.SEG_STATUS_PREPARED, 0),
       count_by_status.get(Segment.SEG_STATUS_RUNNING, 0),
       count_by_status.get(Segment.SEG_STATUS_COMPLETE, 0),
       bins_populated, nbins, bins_populated/nbins*100,
       weights.sum(),
       cputimes.sum(),
       cputimes.mean(),
       min(cputimes) if cputimes.size else 0.0,
       max(cputimes) if cputimes.size else 0.0,
       iter_wallclock))
    
        self.output_stream.write('WE simulation:\n')
        self.output_stream.write('''\
  Completed segments:             %d
  Total CPU (s):                  %s
  Average throughput (segs/hour): %g
  Minimum throughput:             %g
  Maximum throughput:             %g        
''' % (sim_segs,
       sim_cpu,
       throughputs.mean(),
       min(throughputs),
       max(throughputs)))

        if datamgr.is_propagation_complete():
            self.output_stream.write('\nPropagation complete; ready for WE\n')

            
    def cmd_nextseg(self, args):
        parser = self.make_parser(description = \
    """Request information about the next simulation segment, updating the
    simulation state so that other calls to ``nextseg`` yield different
    segments.  Emits output suitable for ``eval`` in bash.""")
        parser.add_option('-p', '--pretend', dest='pretend', action='store_true',
                          help='Do not actually modify the simulation state. '
                              +'NOTE: subsequent calls may not yield the same '
                              +'information')
        parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
                          help='Print output on stderr as well as stdout')
        self.add_state_option(parser)
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
    
        nextseg = self.datamgr.get_next_segment(pretend = bool(opts.pretend))
        
        if nextseg is None:
            self.error_stream.write('Propagation complete for this iteration\n')
            self.exit(EX_ERROR)
        
        vars = {ENV_STATE: self.datamgr.source,
                ENV_CURRENT_ITER: self.datamgr.we_sim.current_iteration,
                ENV_CURRENT_SEG_ID: nextseg.seg_id,
                ENV_CURRENT_SEG_DATA_REF: nextseg.data_ref}
        
        if nextseg.p_parent is not None:
            vars[ENV_PARENT_ITER] = nextseg.p_parent.we_iter
            vars[ENV_PARENT_SEG_ID] = nextseg.p_parent.seg_id
            if nextseg.p_parent.data_ref:
                vars[ENV_PARENT_SEG_DATA_REF] = nextseg.p_parent.data_ref
            
        for (var, val) in ((key[19:], val) 
                           for key,val in self.datamgr.config.iteritems()
                           if key.startswith('additional_environ')):
            vars[var] = val
        
        if opts.verbose:
            for (var, val) in ((var, vars[var]) for var in sorted(vars)):
                self.error_stream.write("%s='%s'\n" % (var, val))
            self.error_stream.flush()
        
        for (var,val) in vars.iteritems():
            self.output_stream.write("%s='%s';" % (var, val))
            
        for var in vars:
            self.output_stream.write('export %s;' % var)
            
        self.output_stream.write('\n')
        self.output_stream.flush()
            
    def cmd_updateseg(self, args):
        from we.util import parse_elapsed_time
        from datetime import datetime
        
        parser = self.make_parser(description=\
    """Record segment results, including progress coordinate value.  Use the
    option -c/--complete to mark the segment as completed and its data fully
    stored.""")
        self.add_state_option(parser)
        self.add_seg_id_option(parser)
        parser.add_option('-v', '--verbose', dest='verbose',
                          action = 'store_true',
                          help = 'be descriptive')
        parser.add_option('-t', '--cputime', dest='cputime',
                          help='CPU time used for this segment')
        parser.add_option('-P', '--pcoord', dest='pcoord',
                          help='Progress coordinate. Separate components with '
                              +'whitespace or a comma')
        parser.add_option('-c', '--complete', dest='mark_complete',
                          action='store_true',
                          help='Mark this segment as complete.  After using '
                              +'this option, the segment data becomes '
                              +'immutable')
        
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        self.load_current_segment(opts)
        
        if self.current_segment.status == Segment.SEG_STATUS_COMPLETE:
            self.error_stream.write('segment already marked complete\n')
            self.exit(2)
        
        if opts.verbose:
            def vprint(msg):
                self.output_stream.write(msg + '\n')
        else:
            def vprint(msg):
                pass
        
        updates = {}
        
        vprint('updating segment %s' % self.current_segment.seg_id)
                
        if opts.cputime:
            updates['cputime'] = parse_elapsed_time(opts.cputime)
            vprint(' -> cputime = %f s' % updates['cputime'])
        
        if opts.pcoord:
            rePCSep = re.compile(r'\s+|,')
            pcoord_strings = rePCSep.split(opts.pcoord)
            if len(pcoord_strings) != len(self.datamgr.we_sim.bins.boundaries):
                self.error_stream.write(('%d progress coordinate elements required '
                                  +'for this simulation (%d provided)\n') 
                                  % (len(sim.bin_limits), len(pcoord_values)))
                self.exit(EX_USAGE_ERROR)
            else:            
                pcoord_values = [float(pcs) for pcs in pcoord_strings]
                updates['final_pcoord'] = tuple(pcoord_values)
                vprint(' -> progress coord: %r' % updates['final_pcoord'])
    
        if opts.mark_complete:
            updates['status'] = Segment.SEG_STATUS_COMPLETE
            updates['endtime'] = datetime.now()
            vprint(' -> marked as complete')
        
        self.datamgr.update_segment(self.current_segment, updates)
        
    def cmd_we(self, args):
        parser = self.make_parser(description=\
    """Perform weighted-ensemble bin/split/merge. Propagation must be complete
    for the current set of segments.""")
        parser.add_option('-p', '--pretend', dest='pretend', 
                          action='store_true',
                          help='Do not actually modify the simulation state '
                              +'(mainly useful for debugging)'
                          )
        self.add_state_option(parser)
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        try:        
            self.datamgr.we_sim.run_we(opts.pretend)
        except PropagationIncompleteError:
            self.error_stream.write('WE unavailable until propagation complete\n')
            self.exit(EX_STATE_ERROR)
        if not opts.pretend:
            self.datamgr.save_state()
    
    def cmd_cstatus(self, args):
        parser = self.make_parser(usage_tail='{-p VALUE|-t COND}',
                              description=\
    """Print (-p) or test (-t) a simulation condition.  Output of -p is
       suitable for assignment to environment variables.  Running with -t
       causes an exit code of 0 if the condition is true or 1 if the 
       condition is false.""")
        self.add_state_option(parser)
        self.add_seg_id_option(parser)
        parser.add_option('-p', '--print', dest='value',
                          help='print simulation value VALUE')
        parser.add_option('-t', '--test', dest='cond',
                          help='test simulation condition COND')
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        
        if not (opts.value or opts.cond):
            self.error_stream.write('either a value to print (-p) or a condition to '
                             +'test (-t) must be provided\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        elif (opts.value and opts.cond):
            self.error_stream.write('only one of -p and -t may be provided\n')
            parser.print_help(self.error_stream)
            self.exit(EX_USAGE_ERROR)
        elif opts.value:
            if opts.value == 'n_segs':
                segs = self.datamgr.get_segments(self.datamgr.we_sim.current_iteration)
                self.output_stream.write(str(len(segs)))
            elif opts.value == 'n_segs_remaining':
                segs = self.datamgr.get_segments(self.datamgr.we_sim.current_iteration)
                n_remaining = len([seg for seg in segs if seg.status != Segment.SEG_STATUS_COMPLETE])
                self.output_stream.write(str(n_remaining))
            elif opts.value == 'next_seg_id':
                seg = self.datamgr.get_next_segment(pretend = True)
                if seg:
                    self.output_stream.write(str(seg.seg_id))
            elif opts.value == 'we_iter':
                self.output_stream.write(str(self.datamgr.we_sim.current_iteration))
            else:
                self.error_stream.write('unrecognized VALUE\n')
                self.exit(EX_USAGE_ERROR)
            self.exit(0)
        elif opts.cond:
            if opts.cond == 'segs_remaining':
                segs = self.datamgr.get_segments(self.datamgr.we_sim.current_iteration)
                n_remaining = len([seg for seg in segs if seg.status != Segment.SEG_STATUS_COMPLETE])
                if n_remaining: 
                    self.exit(0)
                else: 
                    self.exit(1)
            else:
                self.error_stream.write('unrecognized COND\n')
                self.exit(EX_USAGE_ERROR)
                
    def cmd_scrub(self, args):
        parser = OptionParser(description=\
    """Remove database entries related to crashed segments.""")
        self.add_state_option(parser)
        (opts, args) = parser.parse_args(args)
        self.load_state(opts)
        self.datamgr.scrub_crashed_segments(self.datamgr.we_sim.current_iteration)
                    
if __name__ == '__main__':
    WEUberTool().run()
