import os, sys, subprocess, time, datetime, tempfile, string
import re
import cPickle as pickle
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
from wemd.core import Segment
from wemd.util.config_dict import ConfigError
from wemd.backend_drivers import BackendDriver

import math
from copy import copy

class ODLDBackendDriver(BackendDriver):
    re_pot_entry = re.compile(r'^\s*(\S+)\s*,\s*\(([^\)]+)\)\s*,\s*\(([^\)]+)\)')
    re_split_coords = re.compile('\s*,\s*')
    
    def __init__(self):
        super(ODLDBackendDriver,self).__init__()
        self.ld_variance = None
        self.ndim = None
        self.nsteps = None
        self.initial_pcoord = None
        self.potential_heights = None
        self.potential_centers = None
        self.potential_variances = None
        
    def runtime_init(self, runtime_config):
        try:
            self.ld_variance = self.sim_config['backend.odld.ld_variance']
        except (KeyError,TypeError):
            log.info('skipping backend configuration for now')
        else:
            self.ndim = self.sim_config['backend.odld.ndim']
            self.nsteps = self.sim_config['backend.odld.nsteps']
            self.initial_pcoord = self.sim_config['wemd.initial_pcoord']
            self.potential_heights = self.sim_config['backend.odld.potential_heights']
            self.potential_centers = self.sim_config['backend.odld.potential_centers']
            self.potential_variances = self.sim_config['backend.odld.potential_variances']
        
    def sim_init(self, sim_config, sim_config_src):
        self.sim_config = sim_config
        
        sim_config_src.require_all(('backend.odld.ld_stdev','backend.odld.ndim', 'backend.odld.nsteps'))
        
        sim_config['backend.odld.ld_variance'] = self.ld_variance = sim_config_src.get_float('backend.odld.ld_stdev')**2
        sim_config['backend.odld.ndim'] = self.ndim = sim_config_src.get_int('backend.odld.ndim')
        sim_config['backend.odld.nsteps'] = self.nsteps = sim_config_src.get_int('backend.odld.nsteps')
        self.initial_pcoord = sim_config['wemd.initial_pcoord']
        
        pot_heights = []
        pot_pos = []
        pot_vars = []
        for key in [key[:] for key in sim_config_src if key.startswith('backend.odld.potential')]:
            pot_line = sim_config_src[key]
            m = self.re_pot_entry.match(pot_line)
            if not m:
                raise ConfigError('malformed potential entry %r' % pot_line)
            height, pos_str, std_str = m.groups()
            pot_heights.append(float(height))
            pot_pos.append([float(q) for q in self.re_split_coords.split(pos_str)])
            pot_vars.append([float(s)**2 for s in self.re_split_coords.split(std_str)])
        sim_config['backend.odld.potential_heights'] = self.potential_heights = numpy.array(pot_heights)
        sim_config['backend.odld.potential_centers'] = self.potential_centers = numpy.array(pot_pos)
        sim_config['backend.odld.potential_variances'] = self.potential_variances = numpy.array(pot_vars)
        
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        
        for segment in segments:
            self.pre_segment(segment)
            segment.pcoord = numpy.zeros((self.nsteps+1, self.ndim), numpy.float64)
            
            if segment.n_iter == 1 or segment.p_parent is None:
                segment.pcoord[0] = copy(self.initial_pcoord)
            else:
                segment.pcoord[0] = copy(segment.p_parent.pcoord[-1])
                
            segment.starttime = datetime.datetime.now()
            init_walltime = time.time()
            init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            
            self.ld_propagate(segment)

            final_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            final_walltime = time.time()
            segment.endtime = datetime.datetime.now()
            segment.walltime = final_walltime - init_walltime
            segment.cputime = final_cputime - init_cputime
            
            segment.status = Segment.SEG_STATUS_COMPLETE
            
            self.post_segment(segment)
        
    def ld_propagate(self, segment):
        pcoord = segment.pcoord
        ld_variance = self.ld_variance
        potential_heights = self.potential_heights
        potential_centers = self.potential_centers
        potential_variances  = self.potential_variances
        gauss_random = (2*math.pi*ld_variance)**0.5 * numpy.random.normal(0.0, ld_variance**0.5, (self.nsteps, self.ndim))
        
        if potential_centers.size:
            for i in xrange(1, self.nsteps+1):
                pot_dist = pcoord[i-1] - potential_centers
                fric_term = ld_variance/2 * (potential_heights * pot_dist/potential_variances
                                             * numpy.exp(-pot_dist**2/(2*potential_variances))).sum(0)
                pcoord[i] = pcoord[i-1] - fric_term + gauss_random[i-1]
        else:
            for i in xrange(1, self.nsteps+1):
                pcoord[i] = pcoord[i-1] + gauss_random[i-1]
