from __future__ import division, print_function

import os, sys

if sys.version_info[0] < 3 and sys.version_info[1] < 7:
    sys.stderr.write('wemd requires at least Python version 2.7\n')
    sys.exit(1)

import logging
log = logging.getLogger('wemd_cli')
import argparse

import wemd
from wemd.util.config_dict import ConfigDict
from wemd.util import extloader

import numpy, operator, itertools

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
        self._system_driver = sim_manager.system_driver

    def get_new_weights(self):                            
        self.calcTransitionMatrix()        
        Cij = self.calcCountMatrix()
        eps = numpy.finfo(numpy.double).eps
        
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
            ibin = bins.index(target_region) #this is in terms of the original Cij
            new_ibin = oldindex.index(ibin)  #this is in terms of the Rij
            flat_target_regions.append(new_ibin)
            
        n_targets = len(flat_target_regions)
        
        bindim = self._nstates
        print("nstates:%r Q.shape:%r"%(self._nstates,Q.shape))
        
        SSP, ActiveBins = self.steady_state_approximation(Q, bindim, flat_target_regions)
        new_weights = self.prepare_new_weights(SSP,(len(ActiveBins),))
        
        #now remap new_weights onto the bin indices
        MappedActiveBins = []
        for ibin in ActiveBins:
            MappedActiveBins.append(oldindex[ibin])
        
        mapped_new_weights = numpy.zeros(len(bins))
        for i in xrange(0,len(new_weights)):
            mapped_new_weights[oldindex[i]] = new_weights[i]
        
        binnorm = 0.0
        for ibin in MappedActiveBins:
            binnorm += bins[ibin].weight
        
        #scale weights so sum is 1.0 for active bins (and 0 for inactive bins)
        for ibin in xrange(0,len(bins)):
            if ibin in MappedActiveBins:
                mapped_new_weights[ibin] *= 1.0/binnorm
            else:
                mapped_new_weights[ibin] = 0.0

        return mapped_new_weights
    
    def calcTransitionMatrix(self):
        # Determine first and last iteration in the database
        self._miniter = 1
        self._maxiter = self._data_manager.current_iteration - 1 #the last complete iteration
        self._nstates = numpy.max(self._assignments) + 1

    def calcCountMatrix(self):
        """ Calculate the count matrix 
            Returns
            Cij - the count matrix    
        """
        region_set = self._system_driver.region_set
        bins = region_set.get_all_bins()
        self._nstates = len(bins)
        self._Cij = numpy.zeros((self._nstates, self._nstates))
        
        for iter in xrange(self._miniter,self._maxiter):
            parent_segs = self._data_manager.get_segments(iter)
            parent_pcoords = [s.pcoord[-1] for s in parent_segs] 
            parent_seg_ids = [s.seg_id for s in parent_segs]
            parent_weights = [s.weight for s in parent_weights]
            del parent_segs
            child_segs = self._data_manager.get_segments(iter + 1)
            child_weights = [s.weight for s in child_segs]
            child_pcoords = [s.pcoord[-1] for s in child_segs] 
            child_seg_ids = [s.seg_id for s in child_segs]
            child_p_parent_ids = [s.p_parent_id for s in child_segs]
            del child_segs
            
            #bin child and corresponding p parent to determine transition matrix
            child_bins = region_set.map_to_bins(child_pcoords)
            parent_bins = region_set.map_to_bins(parent_pcoords)
            
            child_linbin = [bins.index(child_bin) for child_bin in child_bins]
            parents_linbin = [bins.index(parent_bin) for parent_bin in parent_bins]
            
            #match child with parent and update Cij
            for j in xrange(0,len(child_linbin)):
                Cj = child_linbin[j] #final state 
                jweight = child_weight[j]               
                pid = child_p_parent_ids[j]
                i = parent_seg_ids.index(pid)
                iweight = parent_weights[i]
                Ci = parents_linbin[i]
                self._Cij[Ci,Cj] += jweight / iweight
                 
        return self._Cij            

    def reduceCountMatrix(self, Cij):        
        #Determine which columns/rows are empty
        nonempty = range(0,Cij.shape[0])

        for i in xrange(0,Cij.shape[0]):
            if (Cij[i,:] < eps).all() and (Cij[:,i] < eps).all():
                nonempty.pop(i)

        print("Cij:%r"%Cij)
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
        # nsolve = nstates - ntarget
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
        
        new_weights = numpy.zeros(dim)
    
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

    ss = steady_state(sim_manager)
    new_weights = ss.get_new_weights()
    #reweight bins
    dm = sim_manager.data_manager
    n_iter = dm.current_iteration

    region_set = sim_manager.system_driver.region_set
    bins = region_set.get_all_bins()    
    
    #actually reweight the particles in the bins                
    for ibin in xrange(0,len(bins)):
        bins[ibin].reweight(new_weights[ibin]) 
    
    #update the particles in the data manager
    particles = []
    for bin in bins:
        particles.extend(bin)
    
    dm.update_segments(particles)
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