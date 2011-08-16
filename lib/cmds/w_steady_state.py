from __future__ import division, print_function

'''Enhanced steady state sampling (Bhatt, Zhang, Zuckerman, 2010) for WE,
originally implemented by Joshua L. Adelman (2010), re-implemented for this version
of WEMD by Joe Kaus (February 2011) and updated by Joshua L. Adelman (May 2011).'''

import os, sys

import logging
log = logging.getLogger('w_steady_state')
import argparse

import wemd, wemdtools
from wemd import Segment
import numpy, operator, itertools
from itertools import izip
from math import ceil, floor, log10

class steady_state(object):
    """Calculate the Transition Probability Matrix for lag time [1]  
    """
    def __init__(self, sim_manager, args):
        self._sim_manager = sim_manager
        self._data_manager = sim_manager.data_manager

        self._system_driver = sim_manager.system
        self.miniter = args.start_iter
        self.maxiter = args.stop_iter
        self.count_adjust = args.count_adjust
        self.debug = None
        if args.ss_h5 is not None:
            import h5py
            self.debug = h5py.File(args.ss_h5,'w')
            
        self.quiet_mode = args.quiet_mode
        self.symmetrize = args.symmetrize

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
        for target_region in target_regions:
            for ibin in xrange(0,len(bins)):
                if target_region.bin is bins[ibin]:
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
            
        #assert not (mapped_new_weights < 0).any()
                
        return mapped_new_weights

    def calcCountMatrix(self):
        """ Calculate the count matrix 
            Returns
            Cij - the count matrix    
        """
        if sys.stdout.isatty() and not self.quiet_mode:
            sys.stdout.write('Calculating Count Matrix\n')
        
        
        region_set = self._system_driver.region_set
        bins = region_set.get_all_bins()
        self._nstates = len(bins)
        self._Cij = numpy.zeros((self._nstates, self._nstates))

        for n_iter in xrange(self.miniter,self.maxiter):
            if sys.stdout.isatty() and not self.quiet_mode:
                sys.stdout.write('\rProcessing iteration {:d} of {:d}'.format(n_iter,self.maxiter))
                sys.stdout.flush()

            segs = self._data_manager.get_segments(n_iter)
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
        
        if self.symmetrize:
            self._Cij = 0.5*(self._Cij + self._Cij.T)
        
        if sys.stdout.isatty() and not self.quiet_mode:
            sys.stdout.write('\n')
            
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
  
def cmd_steady_state(sim_manager, args):
    sim_manager.load_data_manager()
    sim_manager.data_manager.open_backing()
        
    sim_manager.load_system_driver()
    sim_manager.load_we_driver()
    sim_manager.load_propagator()
    
    dm = sim_manager.data_manager    
    n_iter = dm.current_iteration
    
    # Adjust iteration start/stop to ensure that they are in storage
    args.start_iter = max(1,args.start_iter)
    maxiter = sim_manager.data_manager.current_iteration - 1
    if args.stop_iter:
        args.stop_iter = min(maxiter,args.stop_iter)
    else:
        args.stop_iter = maxiter

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
        
    ss = steady_state(sim_manager, args)
    new_weights = ss.get_new_weights()
    
    # Calculate pre-reweighting weights
    region_set = sim_manager.system.region_set
    n_bins = new_weights.size
    binprobs = numpy.empty((n_bins,), numpy.float64)
    binprobs[:] = 0.0
    try:
        iter_group = sim_manager.data_manager.get_iter_group(dm.current_iteration)
    except KeyError:
        # current_iteration is set but segments don't exist
        raise
       
    weights = iter_group['seg_index'][:, 'weight']

    # Read only the initial point
    pcoords = iter_group['pcoord'][:,0,:]
    bin_indices = region_set.map_to_all_indices(pcoords)
    for (seg_id, ibin) in enumerate(bin_indices):
        binprobs[ibin] += weights[seg_id]

    # Write new weights to storage
    if args.no_reweight is False:
        if sys.stdout.isatty() and not args.quiet_mode:
            sys.stdout.write('\n')
            sys.stdout.write('Reweighting segments in storage\n')
        #reweight bins
        region_set = sim_manager.system.region_set
    
        region_set.clear()
        bins = region_set.get_all_bins()    
    
        new_seg_coords = [segment.pcoord[0] for segment in new_segments]
        for (seg, bin) in izip(new_segments, sim_manager.system.region_set.map_to_bins(new_seg_coords)):
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
        
        dm.prepare_iteration(n_iter+1, new_segments)

        sim_manager.data_manager.current_iteration = n_iter + 1
        dm.flush_backing()
        dm.close_backing()
        
    # Write original and new bin weights to output file
    if not args.suppress_headers:
        args.output_file.write('''\
# WE re-weighting analysis
# ----
''')
    
        args.output_file.write('''\
# column 0:  bin index
# column 1:  Weight before re-weighting
# column 2:  Weight calculated by analysis
''')
    
        max_bin_width = int(ceil(log10(n_bins)))
        prob_width = min(12, args.precision+5)
        prob_fmt = '{:d}.{:d}e'.format(prob_width, args.precision)
        for ibin in xrange(0, n_bins):
            prew = binprobs[ibin]
            postw = new_weights[ibin]
            args.output_file.write(( '{ibin:<{max_bin_width}d}    {prew:{prob_fmt}}    {postw:{prob_fmt}}\n')
                                    .format(ibin=ibin, prew=prew, postw=postw, max_bin_width=max_bin_width, prob_fmt=prob_fmt))
    
    sys.exit(0)        

#to allow others to use the steady_state class if desired 
if __name__ == "__main__":
    parser = wemd.rc.common_arg_parser('w_steady_state', description='''Reweight a simulation to help achieve steady state''')

    eps = numpy.finfo(numpy.float64).eps
    # Subset options
    parser.add_argument('-b', '--begin', '--start', dest='start_iter', type=int, default=1,
                        help='Begin collecting statistics at iteration START_ITER (default: first iteration)')
    parser.add_argument('-e', '--end', '--stop', dest='stop_iter', type=int,
                        help='Stop collecting statistics at iteration STOP_ITER (default: last completed iteration)')
    # Count adjust options
    parser.add_argument('--count-adjust', dest='count_adjust', type=numpy.float64, default=10.0*eps,
                            help='adjust the count matrix to prevent singular matrix or negative probabilities (default: 10*machine epsilon)\
                            If there are problems, try starting with this set to one.')
    parser.add_argument('--symmetrize', dest='symmetrize', action='store_true',
                            help='''Symmetrize count matrix (Equilibrium only) 
                            (default: Count matrix is not symmetrized)''')
    # Reweighting option
    parser.add_argument('--noreweight', dest='no_reweight', action='store_true',
                        help='''Do not reweight simulation; only output calculated values to OUTPUT_FILE 
                        (default: reweight segment weights in data storage)''')
    # Output
    parser.add_argument('-o', '--output', dest='output_file', type=argparse.FileType('wt'), default=sys.stdout,
                        help='Store bin weights before and after reweighting in OUTPUT_FILE (default: write to standard output).')
    parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                        help='Do not write headers to output files (default: write headers).')
    parser.add_argument('-p', '--precision', dest='precision', type=int, 
                        help='Number of significant figures for probability display (default: 6)',
                        default=6)
    parser.add_argument('--ss-h5', dest='ss_h5', default=None, help='write debugging info to this h5 file')
    parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                        help='''Do not emit periodic status messages (default: emit status messages if standard output
                        is a terminal)''')
        
    # Parse command line arguments
    args = parser.parse_args()
    wemd.rc.config_logging(args)
    runtime_config = wemd.rc.read_config(args.run_config_file)
    runtime_config.update_from_object(args)
    sim_manager = wemd.rc.load_sim_manager(runtime_config)
    
    if sys.stdout.isatty() and not args.quiet_mode:
        sys.stdout.write('Pushing new weights to storage: {}\n'.format(not args.no_reweight))
        sys.stdout.write('Symmetrizing count matrix: {}\n'.format(args.symmetrize))

    try:
        cmd_steady_state(sim_manager, args)
    except Exception as e:
        # The following won't show up if the log isn't set up properly
        log.error(str(e))
        sys.stderr.write('ERROR: {!s}\n'.format(e))
        import traceback
        traceback.print_exc()
        sys.exit(1)    
