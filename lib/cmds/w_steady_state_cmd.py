from __future__ import division, print_function

'''Enhanced steady state sampling (Bhatt, Zhang, Zuckerman, 2010) for WE,
originally implemented by Josh Adelman (2010), re-implemented for this version
of WEMD by Joe Kaus (February 2011).'''

import os, sys

import logging
log = logging.getLogger('w_steady_state')
import argparse

import wemd
from wemd import Segment
import numpy, operator, itertools
from itertools import izip

class steady_state(object):
    """Calculate the Transition Probability Matrix for lag time [1]  
    """
    def __init__(self, sim_manager, count_adjust = 0, debug = None):
        self._sim_manager = sim_manager
        self._data_manager = sim_manager.data_manager

        self._system_driver = sim_manager.system
        # Determine first and last iteration in the database
        self.miniter = 1
        self.maxiter = self._data_manager.current_iteration - 1 #the last complete iteration
        self.count_adjust = count_adjust
        self.debug = None
        if debug is not None:
            import h5py
            self.debug = h5py.File(debug,'w')

    def get_new_weights(self):                    
        Cij = self.calcCountMatrix()
        if self.debug is not None:        
            self.debug.create_dataset('Cij',data=Cij)
        
        eps = numpy.finfo(type(Cij[0,0])).eps

        #Reduce Cij to remove columns/row which are empty
        #oldindex[Ri] == Ci (ie it maps the reduced i,j to the original i,j which correspond to the original bins)        
        Rij, oldindex = self.reduceCountMatrix(Cij)
        #prevents a problem due to a low # of transitions
        Rij += self.count_adjust

        if self.debug is not None:          
            self.debug.create_dataset('Rij',data=Rij)
            
        Q = self.convCount2ProbMatrix(Rij)
        
        if self.debug is not None:             
            self.debug.create_dataset('Q',data=Q)
                
        region_set = self._system_driver.region_set
        #Cij -> i or j is the index into bins
        bins = region_set.get_all_bins()
        
        target_regions = self._system_driver.target_states
        
        flat_target_regions = []
        for (name, target_region) in target_regions:
            for ibin in xrange(0,len(bins)):
                if target_region is bins[ibin]:
                    break

            if ibin in oldindex: #it is possible that the target region was removed (ie if no recycling has occurred)
                new_ibin = oldindex.index(ibin)  #this is in terms of the Rij
                flat_target_regions.append(new_ibin)
            
        n_targets = len(flat_target_regions)
        
        bindim = Rij.shape[0]
        
        SSP, ActiveBins = self.steady_state_approximation(Q, bindim, flat_target_regions)
        if self.debug is not None:     
            self.debug.create_dataset('SSP',data=SSP)
        
        #now remap new_weights onto the bin indices
        MappedActiveBins = []
        for ibin in ActiveBins:
            MappedActiveBins.append(oldindex[ibin])

        mapped_new_weights = numpy.zeros(len(bins))
        for i in xrange(0,len(SSP)):
            mapped_new_weights[MappedActiveBins[i]] = SSP[i]
            #else it will be 0
  
        binnorm = 0.0
        for ibin in MappedActiveBins:
            binnorm += bins[ibin].weight
            
        if self.debug is not None:     
            self.debug.create_dataset('binnorm',data=binnorm)
            if MappedActiveBins:
                self.debug.create_dataset('MappedActiveBins',data=MappedActiveBins)
            self.debug.create_dataset('mapped_new_weights_not_norm',data=mapped_new_weights)

        #scale weights so sum is 1.0 for active bins (and 0 for inactive bins)
        for ibin in xrange(0,len(bins)):
            if ibin in MappedActiveBins:
                mapped_new_weights[ibin] *= 1.0/binnorm
            else:
                mapped_new_weights[ibin] = 0.0

        if self.debug is not None:     
            self.debug.create_dataset('final_weights',data=mapped_new_weights)
            self.debug.close()
                
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

        for iter in xrange(self.miniter,self.maxiter):

            segs = self._data_manager.get_segments(iter)
            weights = [s.weight for s in segs]
            initial_pcoords = [s.pcoord[0] for s in segs] 
            final_pcoords = [s.pcoord[-1] for s in segs] 
            endpoint_types = [s.endpoint_type for s in segs]
            
            #bin child and corresponding p parent to determine transition matrix
            initial_bins = region_set.map_to_bins(initial_pcoords)
            final_bins = region_set.map_to_bins(final_pcoords)
            
            del segs
            
            initial_linbin = []
            for initial_bin in initial_bins:
                for ibin in xrange(0,len(bins)):
                    if bins[ibin] is initial_bin:
                        initial_linbin.append(ibin)

            final_linbin = []
            for final_bin in final_bins:
                for ibin in xrange(0,len(bins)):
                    if bins[ibin] is final_bin:
                        final_linbin.append(ibin)                                        
                                  
            #match child with parent and update Cij
            for j in xrange(0,len(initial_linbin)):
                
                if endpoint_types[j] != Segment.SEG_ENDPOINT_TYPE_CONTINUES and endpoint_types[j] != Segment.SEG_ENDPOINT_TYPE_RECYCLED:    
                    continue
                
                Cj = final_linbin[j] #final state 
                Ci = initial_linbin[j] #initial state
                self._Cij[Ci,Cj] += 1.0
      
        return self._Cij            

    def reduceCountMatrix(self, Cij):        
        #Determine which columns/rows are empty
        nonempty = range(0,Cij.shape[0])
        eps = numpy.finfo(type(Cij[0,0])).eps
                
        for i in xrange(0,Cij.shape[0]):
            if (Cij[i,:] < eps).all() and (Cij[:,i] < eps).all():
                nonempty.pop(nonempty.index(i))

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
  
def cmd_steady_state(sim_manager, count_adjust, ss_h5):
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing()
        
    sim_manager.load_system_driver()
    sim_manager.load_we_driver()
    sim_manager.load_propagator()
    
    dm = sim_manager.data_manager    
    n_iter = dm.current_iteration

    segments = sim_manager.data_manager.get_segments(n_iter)
    
    #If the segments are all incomplete, delete this iteration and re-run we in case the binning has changed
    seg_status = numpy.array([s.status for s in segments],dtype=numpy.uint)
    
    if (seg_status != Segment.SEG_STATUS_COMPLETE).any():
        dm.del_iter_group(n_iter)
        dm.del_iter_summary(n_iter)
        sim_manager.data_manager.current_iteration = n_iter - 1
        n_iter -= 1     
        segments = sim_manager.data_manager.get_segments(n_iter)
               
    new_segments = sim_manager.we_driver.run_we(segments, sim_manager.system.region_set)
    
    sim_manager.data_manager.update_segments(n_iter, segments)
        
    ss = steady_state(sim_manager, count_adjust, ss_h5)
    new_weights = ss.get_new_weights()

    #reweight bins
    region_set = sim_manager.system.region_set
    
    region_set.clear()
    bins = region_set.get_all_bins()    
   
    for (seg, bin) in izip(new_segments, sim_manager.system.region_set.map_to_bins(s.pcoord[0] for s in new_segments)):
        bin.add(seg)
    
    #actually reweight the particles in the bins                
    for ibin in xrange(0,len(bins)):
        bins[ibin].reweight(new_weights[ibin]) 

    # Create new iteration group in HDF5
    for new_segment in new_segments:
        assert new_segment.weight is not None
        assert new_segment.p_parent_id is not None
        assert new_segment.parent_ids is not None and new_segment.p_parent_id in new_segment.parent_ids
        new_segment.n_iter = n_iter+1
        new_segment.status = Segment.SEG_STATUS_PREPARED
        new_segment.n_parents = len(new_segment.parent_ids)
        
    dm.prepare_iteration(n_iter+1, new_segments, 
                                        sim_manager.system.pcoord_ndim, sim_manager.system.pcoord_len, sim_manager.system.pcoord_dtype)

    sim_manager.data_manager.current_iteration = n_iter + 1
    dm.flush_backing()
    dm.close_backing()
    
    sys.exit(0)        

#to allow others to use the steady_state class if desired 
if __name__ == "__main__":
    parser = wemd.rc.common_arg_parser('w_steady_state', description='''Reweight a simulation to help achieve steady state''')

    eps = numpy.finfo(numpy.float64).eps
    
    parser.add_argument('--count-adjust', dest='count_adjust', type=numpy.float64, default=10.0*eps,
                            help='adjust the count matrix to prevent singular matrix or negative probabilities (default: 10*machine epsilon)\
                                 If there are problems, try starting with this set to one.')

    parser.add_argument('--ss-h5', dest='ss_h5', default=None, help='write debugging info to this h5 file')
        
    # Parse command line arguments
    args = parser.parse_args()
    wemd.rc.config_logging(args)
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)

    try:
        cmd_steady_state(sim_manager, args.count_adjust, args.ss_h5)
    except Exception as e:
        # The following won't show up if the log isn't set up properly
        log.error(str(e))
        sys.stderr.write('ERROR: {!s}\n'.format(e))
        import traceback
        traceback.print_exc()
        sys.exit(1)    
