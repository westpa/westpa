from __future__ import division, print_function

import os, sys

if sys.version_info[0] < 3 and sys.version_info[1] < 7:
    sys.stderr.write('wemd requires at least Python version 2.7\n')
    sys.exit(1)

import logging
log = logging.getLogger('wemd_cmd')
import argparse

import wemd
from wemd.util.config_dict import ConfigDict
from wemd.util import extloader

import numpy, operator, itertools
from wemd.types import Segment, Particle
from itertools import izip
# Runtime config file management
ENV_RUNTIME_CONFIG  = 'WEMDRC'
RC_DEFAULT_FILENAME = 'wemd.cfg'

def read_config(filename = None):
    if filename is None:
        filename = RC_DEFAULT_FILENAME
    
    cdict = ConfigDict()
    cdict.read_config_file(filename)
    
    return cdict

def load_sim_manager(runtime_config):
    drivername = runtime_config.get('drivers.sim_manager', 'default')
    if drivername.lower() == 'default':
        from wemd.sim_manager import WESimManager
        return WESimManager(runtime_config)
    else:
        pathinfo = runtime_config.get_pathlist('drivers.module_path')
        return extloader.get_object(drivername,pathinfo)(runtime_config)

class steady_state(object):
    """Calculate the Transition Probability Matrix for lag time [1]  
    """
    def __init__(self, sim_manager):
        self._sim_manager = sim_manager
        self._data_manager = sim_manager.data_manager
        #self._h5file = self._data_manager.h5file
        self._system_driver = sim_manager.system
        # Determine first and last iteration in the database
        self._miniter = 1
        self._maxiter = self._data_manager.current_iteration - 1 #the last complete iteration
        
    def get_new_weights(self):     
        Cij = self.calcCountMatrix()
        eps = numpy.finfo(type(Cij[0,0])).eps
        
        #Reduce Cij to remove columns/row which are empty
        #oldindex[Ri] == Ci (ie it maps the reduced i,j to the original i,j which correspond to the original bins)        
        Rij, oldindex = self.reduceCountMatrix(Cij)
        
        #prevents a problem due to a low # of transitions
        Rij += eps*10 
        Q = self.convCount2ProbMatrix(Rij)
            
        region_set = self._system_driver.region_set
        #Cij -> i or j is the index into bins
        bins = region_set.get_all_bins()
        
        target_regions = self._system_driver.target_states
        
        flat_target_regions = []
        for (name, target_region) in target_regions:
            for ibin in xrange(0,len(bins)):
                if target_region is bins[ibin]:
                    break

            new_ibin = oldindex.index(ibin)  #this is in terms of the Rij
            flat_target_regions.append(new_ibin)
            
        n_targets = len(flat_target_regions)
        
        bindim = self._nstates
        print("nstates:%r Q.shape:%r"%(self._nstates,Q.shape))
        print("flat_target_regions:%r"%flat_target_regions)
        SSP, ActiveBins = self.steady_state_approximation(Q, bindim, flat_target_regions)
        print("SSP:%r"%SSP)
        
        #new_weights = self.prepare_new_weights(SSP,(len(ActiveBins),))
        
        #now remap new_weights onto the bin indices
        MappedActiveBins = []
        for ibin in ActiveBins:
            MappedActiveBins.append(oldindex[ibin])
        print("oldindex:%r"%oldindex)
        mapped_new_weights = numpy.zeros(len(bins))
        for i in xrange(0,len(SSP)):            
            mapped_new_weights[MappedActiveBins[i]] = SSP[i]
            #else it will be 0
        
        for ibin in xrange(0,len(bins)):
            print("ibin:%r binweight:%r"%(ibin,bins[ibin].weight))
            
        binnorm = 0.0
        for ibin in MappedActiveBins:
            binnorm += bins[ibin].weight
        print("mapped_new_weights0:%r"%mapped_new_weights)
        print("MappedActiveBins:%r"%MappedActiveBins) 
        #scale weights so sum is 1.0 for active bins (and 0 for inactive bins)
        for ibin in xrange(0,len(bins)):
            if ibin in MappedActiveBins:
                mapped_new_weights[ibin] *= 1.0/binnorm
            else:
                mapped_new_weights[ibin] = 0.0

        print("mapped_new_weights:%r"%mapped_new_weights)
        return mapped_new_weights

    def calcCountMatrix(self):
        """ Calculate the count matrix 
            Returns
            Cij - the count matrix    
        """
        region_set = self._system_driver.region_set
        bins = region_set.get_all_bins()
        self._nstates = len(bins)
        self._Cij = numpy.zeros((self._nstates, self._nstates))

        #populate the bins
        segments = self._data_manager.get_segments(self._maxiter)
        #last_pcoords = [s.pcoord[-1] for s in last_segs] 
        #last_bins = region_set.map_to_bins(last_pcoords)        
        #for (bin, seg) in zip(last_bins, last_segs):                
        #    bin.add(seg)

        for iter in xrange(self._miniter,self._maxiter):
            print("iter:%r"%iter)
            
            parent_segs = self._data_manager.get_segments(iter)
            parent_pcoords = [s.pcoord[-1] for s in parent_segs] 
            parent_seg_ids = [s.seg_id for s in parent_segs]
            parent_weights = [s.weight for s in parent_segs]

            parent_bins = region_set.map_to_bins(parent_pcoords)        
                            
            del parent_segs

            child_segs = self._data_manager.get_segments(iter + 1)
            child_weights = [s.weight for s in child_segs]
            child_pcoords = [s.pcoord[-1] for s in child_segs] 
            child_seg_ids = [s.seg_id for s in child_segs]
            child_p_parent_ids = [s.p_parent_id for s in child_segs]
                        
            #bin child and corresponding p parent to determine transition matrix
            child_bins = region_set.map_to_bins(child_pcoords)

            del child_segs

            #print("child_bins:%r"%child_bins)
            #print("parent_bins:%r"%parent_bins)
            #print("bins:%r"%bins)
            #print("bins[0]:%r"%bins[0])
            #print("bins[1]:%r"%bins[1])
            #print("bins[0]==bins[1]:%r"%(bins[0] is bins[1]))
            
            child_linbin = []
            for child_bin in child_bins:
                for ibin in xrange(0,len(bins)):
                    if bins[ibin] is child_bin:
                        child_linbin.append(ibin)

            parents_linbin = []
            for parent_bin in parent_bins:
                for ibin in xrange(0,len(bins)):
                    if bins[ibin] is parent_bin:
                        parents_linbin.append(ibin)
                                
            #child_linbin = [bins.index(child_bin) for child_bin in child_bins]
            #parents_linbin = [bins.index(parent_bin) for parent_bin in parent_bins]
            
            #print("child_linbin:%r"%child_linbin)
            #print("parents_linbin:%r"%parents_linbin)
            for i in xrange(0,len(bins)):
                print("w:%r"%bins[i].weight)
                                        
            #match child with parent and update Cij
            for j in xrange(0,len(child_linbin)):
                Cj = child_linbin[j] #final state 
                jweight = child_weights[j]               
                pid = child_p_parent_ids[j]
                if pid < 0: #no parent -> no transition
                    continue
                i = parent_seg_ids.index(pid)
                iweight = parent_weights[i]
                Ci = parents_linbin[i]
                self._Cij[Ci,Cj] += jweight / iweight
                
        print("Cij:%r"%self._Cij)
        return self._Cij            

    def reduceCountMatrix(self, Cij):        
        #Determine which columns/rows are empty
        nonempty = range(0,Cij.shape[0])
        eps = numpy.finfo(type(Cij[0,0])).eps

        print("Cij:%r"%Cij)
                
        for i in xrange(0,Cij.shape[0]):
            if (Cij[i,:] < eps).all() and (Cij[:,i] < eps).all():
                print("nonempty:%r i:%r"%(nonempty,i))
                nonempty.pop(nonempty.index(i))

        print("nonempty:%r" % (nonempty)) 

        Rij = numpy.zeros((len(nonempty),len(nonempty)))
        
        for i in xrange(0,len(nonempty)):
            for j in xrange(0,len(nonempty)):
                Rij[i,j] = Cij[nonempty[i],nonempty[j]]
                
        return Rij, nonempty

    def convCount2ProbMatrix(self, Q):
        """Convert count matrix to a row transition probability matrix"""
        nstates = Q.shape[0]
        QT = Q/numpy.tile(numpy.sum(Q,axis=1).reshape(-1,1),(1,nstates))
        
        ii = numpy.isnan(QT)
        QT[ii] = 0.0
    
        return QT
    
    def steady_state_approximation(self, T, N, targetBinsIndex):
        ntarget = len(targetBinsIndex)          # Number of target states
        nstates = N                             # Number of total states

        # Number of active states
        nsolve = nstates - ntarget
       
        # list of active states
        sactive = sorted(list(set(xrange(nstates)) - set(targetBinsIndex)))
                                    
        W = numpy.zeros((nsolve,nsolve))
        S = numpy.zeros((nsolve,))
        
        #we have n equations, but only n-1 linearly independent equations
        #set one equation to constraint sum(Pi)=1
        S[0] = 1.0
        
        for iiw,iit in zip(xrange(0,W.shape[0]),sactive):
            for jjw,jjt in zip(xrange(0,W.shape[0]),sactive):
                W[iiw,jjw] = T[jjt,iit]
                
        for ii in xrange(0,W.shape[0]):
            W[ii,ii] = 0.0
            for jj in xrange(0,T.shape[0]): #nstates
                if jj != ii:       
                    W[ii,ii] -= T[ii,jj]
    
        #sum(Pi)=1
        for ii in xrange(W.shape[0]):
            W[0,ii] = 1.0
                    
        P = numpy.linalg.solve(W,S)
        
        return P, sactive
    
    def prepare_new_weights(self,Pbins,dim):
        
        new_weights = numpy.zeros(dim[0])
    
        IP = Pbins
    
        for i in xrange(dim[0]):
            new_weights[i] += IP[i]
            
        return new_weights
        
        
def cmd_steady_state(sim_manager):
    #sim_manager.load_work_manager()
    #aux_args = sim_manager.work_manager.parse_aux_args(aux_args)
    
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing()
        
    sim_manager.load_system_driver()
    sim_manager.load_we_driver()

    dm = sim_manager.data_manager
    dm.del_iter_group(dm.current_iteration)
    n_iter = dm.current_iteration - 1

    #populate the bins
    segments = sim_manager.data_manager.get_segments(n_iter)

    # Convert segments into particles representing their endpoints
    particles = [Particle(seg_id = segment.seg_id,
                          weight = segment.weight,
                          p_parent_id = None, # NOT the same as a segment's p_parent_id
                          parent_ids = set(),
                          pcoord = segment.pcoord[-1]) for segment in segments]

    next_iter_particles = sim_manager.we_driver.run_we(particles, sim_manager.system.region_set)

    # Bin particles after running we
    for (particle, bin) in izip(next_iter_particles,
                                sim_manager.system.region_set.map_to_bins(particle.pcoord for particle in next_iter_particles)):
        bin.add(particle)

    ss = steady_state(sim_manager)
    new_weights = ss.get_new_weights()
    #reweight bins
    region_set = sim_manager.system.region_set
    bins = region_set.get_all_bins()    
    
    #actually reweight the particles in the bins                
    for ibin in xrange(0,len(bins)):
        bins[ibin].reweight(new_weights[ibin]) 
    #update the particles in the data manager
    particles = []
    for bin in bins:
        particles.extend(bin)
    
   # Create a new iteration with new weights
    new_segments = []
    for particle in particles:
        if log.isEnabledFor(logging.DEBUG):
            log.debug('processing particle %r' % particle)
        segment = Segment(seg_id = None,
                          n_iter = n_iter+1,
                          weight = particle.weight,
                          status = Segment.SEG_STATUS_PREPARED,
                          )

        # Recall that each particle has only one progress coordinate point; we need
        # the correctly-shaped array
        segment.pcoord = sim_manager.system.new_pcoord_array()
        segment.pcoord[0,:] = particle.pcoord
        #print("particle.pcoord:%r"%particle.pcoord)
        
        if particle.p_parent_id is None:
            # Particle did not result from split or merge, but perhaps (if seg_id is negative) a recycle
            assert len(particle.parent_ids) == 0
            assert particle.seg_id is not None
            segment.p_parent_id = particle.seg_id
            segment.parent_ids = set((particle.seg_id,))
        else:
            # Particle did result from a split or a merge
            assert len(particle.parent_ids) > 0
            assert None not in particle.parent_ids
            assert particle.seg_id is None
            segment.p_parent_id = particle.p_parent_id
            segment.parent_ids = set(particle.parent_ids)

        segment.n_parents = len(segment.parent_ids)
        new_segments.append(segment)

    # Create new iteration group in HDF5
    sim_manager.data_manager.prepare_iteration(n_iter+1, new_segments,
                                        sim_manager.system.pcoord_ndim, sim_manager.system.pcoord_len, sim_manager.system.pcoord_dtype)

    sim_manager.data_manager.current_iteration = n_iter + 1
    dm.update_segments(n_iter + 1, new_segments)
    dm.flush_backing()
    dm.close_backing()
    
    sys.exit(0)        


parser = argparse.ArgumentParser()
parser.add_argument('-r', '--rcfile', metavar='RCFILE', dest='run_config_file',
                    help='use RCFILE as the WEMD run-time configuration file (default: %s)' 
                          % RC_DEFAULT_FILENAME)

# Parse command line arguments
args = parser.parse_args()

# Read runtime configuration file
runtime_config = read_config(args.run_config_file)

# Merge command line arguments into runtime config (for convenience)
runtime_config.update({'args.%s' % key : value for (key,value) in args.__dict__.viewitems() if not key.startswith('_')})

# Load SimManager
sim_manager = load_sim_manager(runtime_config)

try:
    cmd_steady_state(sim_manager)
except Exception as e:
    # The following won't show up if the log isn't set up properly
    log.error(str(e))
    sys.stderr.write('ERROR: {!s}\n'.format(e))
    import traceback
    traceback.print_exc()
    sys.exit(1)    
