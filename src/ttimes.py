from __future__ import division
import os, sys, re
from optparse import OptionParser
import numpy
import h5py
from wemd.analysis.probability import quantiles
from itertools import izip

class TransitionFinder(object):
    reSplitWithSpace = re.compile(r',\s*')
    def __init__(self):
        self.region_names = None
        self.region_edges = None
        self.region_name_index = None
        self.event_counts = None
        self.event_times = {}
        self.fpts = {}
        self.tbs = {}
        self.fpt_quantiles = {}
        self.tb_quantiles = {}
        self.transition_log = None
        self.dist = None
        self.dt = None
        self.ntraj = None
        self.ndim = None
        
        self.save_tb = {}
        self.save_fpt = {}
        self.save_tb_quantiles = {}
        self.save_fpt_quantiles = {}
        
        self.output_stream = sys.stdout
        self.error_stream = sys.stderr
    
    def load_pc_hdf5(self, filename, nodename, weight_nodename = None):
        self.h5file = h5py.File(filename, 'r')
        # This ellipsis is important, as it goes through the overhead of 
        # creating the HDF5 hyperslab; otherwise every single access
        # in the analysis loop will trigger an HDF5 selection creation,
        # slowing things down about 1000-fold
        self.dist = self.h5file[nodename][...]
        self.output_stream.write('file contains %d trajectories of %d frames\n'
                                 % (self.dist.shape[0], self.dist.shape[1]))
        if weight_nodename:
            self.weight = self.h5file[weight_nodename][...]
            assert self.dist.shape[0:2] == self.weight.shape[0:2]
            self.output_stream.write('using weights from node %s\n' % weight_nodename)
        else:
            self.weight = None
            self.output_stream.write('using uniform (1.0) weights\n')
            
        try:
            self.dt = self.h5file[nodename].attrs['timestep']
        except KeyError:
            self.dt = 1.0
        self.ntraj = self.dist.shape[0]
        self.ndim = self.dist.ndim - 2

    def read_config(self, config_filename):
        from ConfigParser import SafeConfigParser, NoSectionError, NoOptionError
        
        config_parser = SafeConfigParser()
        config_file = open(config_filename, 'rt')
        config_parser.readfp(config_file)
        
        # Configure the region names and boundaries
        region_name_string = config_parser.get('regions', 'region_names')
        self.region_names = self.reSplitWithSpace.split(region_name_string)
        region_edges_string = config_parser.get('regions', 'region_edges')
        self.region_edges = [float(v) for v in 
                             self.reSplitWithSpace.split(region_edges_string)]
        if len(self.region_edges) != len(self.region_names) + 1:
            raise ValueError('region names and boundaries do not match')
        self.region_name_index = dict((v,i) for (i,v) in enumerate(self.region_names))
        
        # Configure the input
        source_filename = config_parser.get('input', 'source')
        source_node = config_parser.get('input', 'source_node')
        try:
            source_weight_node = config_parser.get('input', 'source_weight_node')
        except NoOptionError:
            source_weight_node = None
        self.load_pc_hdf5(source_filename, source_node, source_weight_node)
        
        # Configure the output, as much as we can at this point
        try:
            if config_parser.getboolean('output', 'save_transitions'):
                self.transition_log = open('transitions.txt', 'wt')
        except (NoSectionError, NoOptionError):
            pass
        
        # Configure the calculations
        try:
            config_tb = config_parser.get('calc', 'calc_tb')
        except (NoOptionError,NoSectionError):
            pass
        else:
            self._config_calc(config_tb, self.tbs)

        try:
            config_fpts = config_parser.get('calc', 'calc_fpt')
        except (NoOptionError,NoSectionError):
            pass
        else:
            self._config_calc(config_tb, self.fpts)

        try:
            config_tb_quantiles = config_parser.get('calc', 'calc_tb_quantiles')
        except (NoOptionError,NoSectionError):
            pass
        else:
            self._config_calc(config_tb_quantiles, self.tb_quantiles)

        try:
            config_fpt_quantiles = config_parser.get('calc', 'calc_fpt_quantiles')
        except (NoOptionError,NoSectionError):
            pass
        else:
            self._config_calc(config_fpt_quantiles, self.fpt_quantiles)

            
        for name in ('save_tb', 'save_fpt', 'save_tb_quantiles',
                     'save_fpt_quantiles'):
            basename = name.split('_', 1)[-1]
            try:
                self._config_save(config_parser.get('output', name),
                                  getattr(self, name),
                                  basename,
                                  'txt')
            except (NoOptionError,NoSectionError):
                pass
        
        self.config_parser = config_parser

    def _config_calc(self, config_string, storage):
        config_string = config_string.strip()
        if config_string.lower() == 'all':
            for irr1 in xrange(0, len(self.region_names)):
                for irr2 in xrange(0, len(self.region_names)):
                    if abs(irr1-irr2) > 1:
                        storage[irr1, irr2] = []
        else:            
            calcs = self.reSplitWithSpace.split(config_string)
            for calc in calcs:
                (region1, region2) = calc.split('->')
                region1 = region1.strip()
                region2 = region2.strip()
                storage[self.region_name_index[region1],
                        self.region_name_index[region2]] = []
                        
    def _config_save(self, config_string, storage, basename, ext):
        config_string = config_string.strip()
        if config_string.lower() == 'all':
            for irr1 in xrange(0, len(self.region_names)):
                for irr2 in xrange(0, len(self.region_names)):
                    if abs(irr1-irr2) > 1:
                        storage[irr1, irr2] = '%s_%s-%s.%s' % (basename,
                                                               self.region_names[irr1],
                                                               self.region_names[irr2],
                                                               ext)
        else:
            saves = self.reSplitWithSpace.split(config_string)
            for save in saves:
                (region1, region2) = save.split('->')
                region1 = region1.strip()
                region2 = region2.strip()
                storage[self.region_name_index[region1],
                        self.region_name_index[region2]] = '%s_%s-%s.%s' \
                                                            % (basename,
                                                               region1,
                                                               region2,
                                                               ext)
       
    def find_transitions(self):
        nreg = len(self.region_names)
        
        event_count = numpy.zeros((nreg,nreg), numpy.int64)
        region_edges = self.region_edges
        region_names = self.region_names
        event_times = self.event_times
        transition_log = self.transition_log
        fpts = self.fpts
        tbs = self.tbs
        dist = self.dist
        dt = self.dt
        weight = self.weight
        
        for itraj in xrange(0, self.ntraj):
            completion_times = numpy.zeros((nreg,nreg), numpy.int64) - 1
            
            # Assign the particle's initial region
            for ir in xrange(0, nreg):
                q = dist[itraj,0,0]
                if region_edges[ir] <= q < region_edges[ir+1]:
                    last_iregion = ir
                    break
            else:
                raise ValueError('particle is outside specified coordinate space')
            
            for it in xrange(1, dist.shape[1]):
                if weight is None:
                    w = 1.0
                else:
                    w = weight[itraj, it]
                
                q = dist[itraj,it,0]
                for ir in xrange(0, nreg):
                    lb = region_edges[ir]
                    ub = region_edges[ir+1]
                    
                    if lb <= q < ub:
                        iregion = ir
                        break
                else:
                    raise ValueError('particle is outside specified coordinate space')
                
                if iregion != last_iregion:        
                    # Crossing event...do all the real work here
                
                    completion_times[last_iregion, iregion] = it
                    event_count[last_iregion, iregion] += 1
                    #try:
                    #    event_times[last_iregion, iregion].append(it)
                    #except KeyError:
                    #    pass
                    
                    if transition_log:
                        transition_log.write('%d    %g    %s -> %s\n' 
                                             % (it, w,
                                                region_names[last_iregion], 
                                                region_names[iregion]) )
                    
                    for irr in xrange(1, nreg):
                        if iregion > last_iregion:
                            trans_iregion = last_iregion - irr
                            tb_iregion = trans_iregion + 1
                        else:
                            trans_iregion = last_iregion + irr
                            tb_iregion = trans_iregion - 1
                            
                        if 0 <= trans_iregion < nreg and trans_iregion != iregion:
                            if completion_times[iregion, last_iregion] > completion_times[trans_iregion, last_iregion]:
                                # Fluctuation between regions
                                #print "  cannot be a %s->%s transition" % (region_names[trans_iregion], region_names[iregion])
                                pass
                            else:
                                # First passage time:  time since last transition ending in starting state
                                if completion_times[iregion, trans_iregion] >= 0:
                                    fpt = it-completion_times[iregion, trans_iregion]
                                    try:
                                        fpts[trans_iregion, iregion].append((fpt,w))
                                    except KeyError:
                                        pass
                            
                                # Event duration: time since last transition originating from starting state
                                # trans_iregion -> trans_iregion +/- 1
                                if completion_times[trans_iregion, tb_iregion] >= 0:
                                    # Since we record the completion time of the exit of the initial region
                                    # subtract one timestep to indicate the initial time of the 
                                    # transition
                                    tb = it-completion_times[trans_iregion, tb_iregion]-1
                                    try:
                                        tbs[trans_iregion, iregion].append((tb,w))
                                    except KeyError:
                                        pass
                                    
                                # Update completion times matrix
                                completion_times[trans_iregion, iregion] = it
                                event_count[trans_iregion, iregion] += 1
                    last_iregion = iregion
                
        self.event_counts = event_count
        for (k,v) in self.fpts.items():
            fpts = numpy.array(v, numpy.float64)
            if fpts.shape[0]:
                fpts[:,0] *= self.dt
                self.fpts[k] = fpts
            else:
                try:
                    del self.fpts[k]
                except KeyError:
                    pass
        for (k,v) in self.tbs.items():
            tbs = numpy.array(v, numpy.float64)
            if tbs.shape[0]:
                tbs[:,0] *= self.dt
                self.tbs[k] = tbs
            else:
                try:
                    del self.tbs[k]
                except KeyError:
                    pass
        #for (k,v) in self.tb_quantiles.iteritems():
        #    self.tb_quantiles[k] = self._calc_quantiles(self.tbs, k)
        #for (k,v) in self.fpt_quantiles.iteritems():
        #    self.fpt_quantiles[k] = self._calc_quantiles(self.fpts, k)
    
    def _calc_quantiles(self, bucket, key):
        data = bucket[key]
        if not data.any(): return numpy.zeros((0, 0), numpy.float64)
        #nq = min(len(data), int(len(data)/5), 100)
        #dq = 1/(nq or 1)
        #qs = numpy.arange(dq,1,dq)
        qs = numpy.arange(0.05, 1.0, 0.05)
        quants = quantiles(data, qs)
        return numpy.column_stack([qs, quants])
    
    def save_results(self):
        for key in self.save_tb:
            of = open(self.save_tb[key], 'wt')
            for t in self.tbs[key][:]:
                of.write('%20.16g %20.16g\n' % tuple(t))
            of.close()
                
        for key in self.save_fpt:
            of = open(self.save_fpt[key], 'wt')
            for t in self.fpts[key][:]:
                of.write('%20.16g %20.16g\n' % tuple(t))
            of.close()
            
        for key in self.save_tb_quantiles:
            of = open(self.save_tb_quantiles[key], 'wt')
            for t in self.tb_quantiles[key][:]:
                of.write('%8f    %20.16g\n' % tuple(t))
            of.close()

        for key in self.save_fpt_quantiles:
            of = open(self.save_fpt_quantiles[key], 'wt')
            for t in self.fpt_quantiles[key][:]:
                of.write('%8f    %20.16g\n' % tuple(t))
            of.close()
            
    def summarize_results(self):
        fpt_keys = set(self.fpts.iterkeys())
        tb_keys = set(self.tbs.iterkeys())
        all_keys = fpt_keys | tb_keys
        common_keys = fpt_keys & tb_keys
        fpts_only = fpt_keys - common_keys
        tbs_only = tb_keys - common_keys
        
        for key in all_keys:
            region1 = self.region_names[key[0]]
            region2 = self.region_names[key[1]]
            self.output_stream.write('%s->%s transitions:\n' 
                                     % (region1, region2))
            
            if key in tb_keys:
                tb = self.tbs[key]
                tb_wnorm = tb[:,1].sum()
                tb_mean = (tb[:,1] * tb[:,0]).sum()/tb_wnorm
                tb_std = ((tb[:,1] * (tb[:,0] - tb_mean)**2).sum()/tb_wnorm)**(0.5)
                #tb_mean = self.tbs[key][:].mean()
                #tb_std = self.tbs[key][:].std()
                tb_sem = tb_std / (tb.shape[0] ** 0.5)
                #tb_min = self.tbs[key][:].min()
                #tb_max = self.tbs[key][:].max()
                #tb_med = self.tbs[key][self.tbs[key].shape[0] // 2]
                self.output_stream.write(
'''    Average t_b: %g
    St.Dev. t_b: %g
    S.E.M.  t_b: %g
''' % (tb_mean, tb_std, tb_sem))
                
            if key in fpt_keys:
                fpt = self.fpts[key]
                fpt_wnorm = fpt[:,1].sum()
                fpt_mean = (fpt[:,1] * fpt[:,0]).sum() / fpt_wnorm
                fpt_std = ((fpt[:,1] * (fpt[:,0] - fpt_mean)**2).sum()/fpt_wnorm)**(0.5)
                #fpt_mean = self.fpts[key][:].mean()
                #fpt_std = self.fpts[key][:].std()
                fpt_sem = fpt_std / (fpt.shape[0] ** 0.5)
                rrate = 1/fpt_mean
                rrate_err = fpt_sem/(fpt_mean ** 2)
                self.output_stream.write(
'''    Average fpt: %g
    St.Dev. fpt: %g
    S.E.M.  fpt: %g
    Reaction rate (1/fpt):  %g
    Error in reaction rate: %g
''' % (fpt_mean, fpt_std, fpt_sem,
       rrate, rrate_err))
            self.output_stream.write('\n')
            
        self.output_stream.write('Event count (R->C):\n%s\n%s\n'
                                 % (self.region_names, self.event_counts))

################################################################################
parser = OptionParser('%prog [OPTIONS] CONFIG_FILE',
                      description = 'analyze trajectories for transition times')
(opts, args) = parser.parse_args()

if len(args) != 1:
    sys.stderr.write('exactly one non-option argument is required\n')
    parser.print_help(sys.stderr)
    sys.exit(2)
    
tf = TransitionFinder()
tf.read_config(args[0])
tf.find_transitions()
tf.summarize_results()
tf.save_results()
