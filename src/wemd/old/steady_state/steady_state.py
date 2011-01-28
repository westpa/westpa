import networkx as nx
import wemd
from wemd import Segment, WESimIter
from wemd.util.wetool import WECmdLineTool
import numpy
from numpy import digitize

import logging
log = logging.getLogger(__name__)

class steady_state(object):
    """Calculate the Transition Probability Matrix for varying lag times
        Arguments
        assignments - list of state assignments for all segments in simulation. 
                    If no list is given, the progress coordinate bins will be used
                    as default states
        runtime_config - runtime configuration file name. Default: run.cfg
        lag - a list of lag times in units of tau to calculate the transition 
                probability matrix. Default: [1]        
    """
    def __init__(self, runtime_config_file='run.cfg', sim_manager = None, data_manager = None, max_iters = 60):
        super(steady_state, self).__init__()
        if not sim_manager or not data_manager:
            raise "Need to pass sim and data managers"

        self._sim_manager = sim_manager
        self._data_manager = data_manager
        self.runtime_config = wemd.rc.read_config(runtime_config_file)   
        self.MAX_ITERS = max_iters # Maximum number of iterations used to assemble a graph
        
    def get_new_weights(self, assignments = None, weights = None, pcoord = None, lag = [1], bb2 = None):
        
        sim_manager = self._sim_manager        
        boundaries = sim_manager.we_driver.bins.boundaries
        target_pcoords = sim_manager.we_driver.target_pcoords
        ndim = sim_manager.we_driver.bins.ndim
        bin_shape = sim_manager.we_driver.bins.shape
        pops = sim_manager.we_driver.bins.population_array()

        bindim = 1
        for idim in xrange(0,ndim):
            bindim *= bin_shape[idim]
            
        print "bindim:%r"%bindim        
                
        self.calcTransitionMatrix(assignments = assignments, weights = weights, pcoord = pcoord,lag = lag, bb2 = bb2)
        
        Cij = self.calcCountMatrix()
        print "Cij:%r" % Cij[1]
        epsilon = numpy.finfo(numpy.double).eps*10
        
        Q = {}
        for k in Cij.keys():
            Cij[k] += epsilon #-> prevents a problem due to a low # of transitions
            Q[k] = self.convCount2ProbMatrix(Cij[k])
            
        print "Q:%r"% Q[1]

        n_targets = len(target_pcoords)

        print target_pcoords
        #assign target_pcoords to bins
        pcoords = numpy.empty((len(target_pcoords)*2, ndim), numpy.float64)
        for ipcoord in xrange(0,len(target_pcoords)):
            pcoord = target_pcoords[ipcoord]
            
            for idim in xrange(0,ndim):                                
                pcoords[ipcoord*2+0][idim] = pcoord[idim][0] + epsilon                             
                pcoords[ipcoord*2+1][idim] = pcoord[idim][1] - epsilon                                               
       
        indices = numpy.empty((pcoords.shape[0], ndim), numpy.uintp)    
        
        for idim in xrange(0, ndim):
            indices[:, idim] = numpy.digitize(pcoords[:, idim], boundaries[idim])
            if (indices[indices[:, idim] == 0].size
               or indices[indices[:, idim] == len(boundaries[idim])].size):
                indices[:, idim] = 1 #this will become bin 0
                log.warning('target_pcoord not in bins, setting to bin 0 for dim:%r'%idim)
            
        # Rebase so that we have a valid index into the array of bins
        indices -= 1

        # gives the ranges of an occupied region of the bins (as an index to the bin boundaries)
        target_regions = [[] for x in xrange(0,int(len(indices)/2))]
        
        for i in xrange(0,len(indices),2):
            low = indices[i]
            upper = indices[i+1]
            for idim in xrange(0,ndim):
                target_regions[int(i/2)].append([low[idim],upper[idim]])
                
        print "target_regions:%r"%target_regions

        flattened_target_regions = []     
        for target_pcoord in target_regions:
            
            print "target_pcoord:%r" % target_pcoord
            target_pcoord_expanded = [[] for x in xrange(0,ndim)]

            for idim in xrange(0,ndim):
                print "target_pcoord:%r"%target_pcoord
                target_pcoord_expanded[idim].extend(range(int(target_pcoord[idim][0]),int(target_pcoord[idim][1])+1))
            print "target_pcoord_expanded:%r"%target_pcoord_expanded
            enum_target_regions = [[]]  
                          
            for x in target_pcoord_expanded:
                t = []
                for y in x:
                    for i in enum_target_regions:
                        t.append(i+[y])
                enum_target_regions = t
            print "enum_target_regions:%r" % enum_target_regions
            for target in enum_target_regions:
                flat = target[-1]
                
                for idim in xrange(1,ndim):
                    dim_size = 1
                    for jdim in xrange(0,idim):
                        dim_size *= bin_shape[-1 - jdim]
               
                    flat += target[-1 - idim] * dim_size
                    
                if flat not in flattened_target_regions:
                    flattened_target_regions.append(flat)                                   

        #add empty bins to flattened_target_regions so no probability is assigned to that location 
        #*** WARNING *** This seems to cause problems (Pi < 0, poor distribution of probability, sometimes unphysical values result)    
        #flat_pops = pops.reshape(bindim,)
        #for i in xrange(0,len(flat_pops)):
        #    print "flat_pops[%r]:%r" % (i,flat_pops[i])
        #    if flat_pops[i] < epsilon and i not in flattened_target_regions:
        #        flattened_target_regions.append(i)

        print "unoccupied regions: \ntarget:%r" % (flattened_target_regions)
                          
        SSP = {}
        ActiveBins = {}
        for k in Q.keys():                                  
            SSP[k],ActiveBins = self.steady_state_approximation(Q[k],bindim,flattened_target_regions) #flattened_target_regions

        #index pops according to sactive
        #flat_pops = pops.reshape(bindim,)
        #sa_pops = [flat_pops[i] for i in ActiveBins]
        
        new_weights = {}
        for k in SSP.keys():
            new_weights[k] = self.prepare_new_weights(SSP[k],(len(ActiveBins),))

        #index new_weights according to original bins -> new_weights_by_bin[i] is weight for bin given by bins[i],bins[i+1]
        empty_bins = []
        flat_pops = pops.reshape(bindim,)
        for i in xrange(0,len(flat_pops)):
            if flat_pops[i] < epsilon:
                empty_bins.append(i)
                        
        new_weights_by_bin = {}
        for k in new_weights.keys():
            new_weights_by_bin[k] = numpy.zeros((bindim,))
            binnorm = 0.0
            for i in xrange(0,bindim):
                if i in ActiveBins and i not in empty_bins:
                    new_weights_by_bin[k][i] = new_weights[k][ActiveBins.index(i)]
                    binnorm += new_weights[k][ActiveBins.index(i)]
            
            #scale weight according to binnorm 
            for i in xrange(0,bindim):
                new_weights_by_bin[k][i] *= 1.0/binnorm
                
            new_weights_by_bin[k] = new_weights_by_bin[k].reshape(bin_shape)
                
        print new_weights_by_bin[1]
        return new_weights_by_bin
                         
    def calcTransitionMatrix(self,assignments=None, weights=None, pcoord=None, lag=[1], bb2=None):

        sim_manager = self._sim_manager
        data_manager = self._data_manager
        runtime_config = self.runtime_config

        # Determine first and last iteration in the database
        self._miniter = self._data_manager.get_first_iter().n_iter # Iteration zero does not contain any segments
        self._maxiter = self._data_manager.get_last_complete_iter().n_iter  # Do not consider last iteration in case incomplete
                                
        last_segs = self._data_manager.get_segments(self._maxiter)
        last_segids = sorted([s.seg_id for s in last_segs])
        self._maxsegid = last_segids[-1]
        
        if assignments == None:
            if bb2 == None:
                print 'Generating assignment list from progress coordinate bins'
                self._assignments = self._generate_assignments_from_pcoord(runtime_config,pcoord)
            
                # save assignments to pickle for next time
                import cPickle
                f = open(runtime_config['data.SIM_NAME'] + '.assignments.pkl','wb')
                cPickle.dump(self._assignments,f,cPickle.HIGHEST_PROTOCOL)
                f.close()
            else:
                print 'Generating assignment list from progress coordinate bins and in-bin drmsd'
                self._assignments = self._generate_assignments_from_pcoord_2d(runtime_config,bb2)
            
                # save assignments to pickle for next time
                import cPickle
                f = open(runtime_config['data.SIM_NAME'] + '.assignments_rmsd2d.pkl','wb')
                cPickle.dump(self._assignments,f,cPickle.HIGHEST_PROTOCOL)
                f.close()
        else:
            self._assignments = numpy.array(assignments)
        
        if weights == None:
            print 'Retrieving weights from database'
            self._weights = self._retrieve_weights()
            
            # save assignments to pickle for next time
            import cPickle
            f = open(runtime_config['data.SIM_NAME'] + '.weights.pkl','wb')
            cPickle.dump(self._weights,f,cPickle.HIGHEST_PROTOCOL)
            f.close()
        else:
            self._weights = numpy.array(weights)
        
        self._lag = lag
        
        #self._nstates = numpy.unique(self._assignments).size
        self._nstates = numpy.max(self._assignments) + 1
        print self._assignments
        print numpy.max(self._assignments) + 1
        # Initialize count matrix 
        self._Cij = {}
        for k in lag:
            self._Cij[k] = numpy.zeros((self._nstates,self._nstates))        
    
    def calcCountMatrix(self,iters=None,reinit=True):
        """ Calculate the count matrix 
            Arguments
            iters - the list of iterations over which to calculate the count matrix
            reinit - boolean flag to reset the count matrix
            
            Returns
            Cij - Dictionary containing the count matrix for each lag time
        
        """
        
        if reinit: # Reinitialize Count Matrix
            for k in self._lag:
                self._Cij[k] = numpy.zeros((self._nstates,self._nstates)) 
        
        maxlag = max(self._lag)
         
        if iters == None:
            # Compute Transition Matrix for the whole graph
            super_iterlist = map(None, *(iter(xrange(self._miniter,self._maxiter)),) * self.MAX_ITERS)
            
            # Remove None padding from last list
            super_iterlist[-1] = tuple(x for x in super_iterlist[-1] if x is not None)
            niters = self._maxiter - self._miniter
            
        else:
            imin = max(min(iters),self._miniter)
            imax = min(max(iters),self._maxiter)    
        
            super_iterlist = map(None, *(iter(xrange(imin,imax)),) * self.MAX_ITERS)
            
            # Remove None padding from last list
            super_iterlist[-1] = tuple(x for x in super_iterlist[-1] if x is not None)
            niters = imax - imin
        
        
        for iterlist in super_iterlist:
            # Construct graph
            G = self._construct_graph(iterlist)
            
            segments = self._data_manager.get_segments(min(iterlist),n_iter_upper = max(iterlist))
        
            print 'Computing contribution for iterations %d to %d' % (min(iterlist),max(iterlist))
        
            for curseg in segments:
                if G.has_node(curseg.seg_id) is False:
                    #print curseg.seg_id
                    continue
            
                if curseg.endpoint_type == Segment.SEG_ENDPOINT_TYPE_RECYCLED:
                    # particle entered target state
                    continue
                elif curseg.endpoint_type == Segment.SEG_ENDPOINT_TYPE_MERGED:
                    # particle termini of branch
                    continue
                elif curseg.endpoint_type == Segment.SEG_ENDPOINT_TYPE_UNKNOWN:
                    # should possibly catch this, but for now just continue
                    continue
                else: #curseg.endpoint_type == Segment.SEG_ENDPOINT_TYPE_CONTINUATION:

                    # Source state and weight
                    statei = self._assignments[curseg.seg_id]
                    sourceweight = curseg.weight
  
                    # tsegs gives a list of all nodes/segs (value) at the lag time (key) away from the source
                    tsegs = {}
                    tsegs[0] = [curseg.seg_id]
                    c = [curseg.seg_id]
            
                    for step in xrange(1,maxlag+1):
                        tt = nx.predecessor(G,curseg.seg_id,cutoff=step).keys()
                        tsegs[step] = list(set(tt).difference(set(c)))
                        c.extend(tt)
                
                    ## If current segment is the termini of a branch, go to next segment    
                    #if tsegs[1] == curseg.seg_id:
                    #    continue
 
                    for step in self._lag:
                        # Traverse subgraph and assign probabilities
                        for n in tsegs[step]:
                            statej = self._assignments[n]
                            #statej_weight = self._data_manager.get_segments(Segment.seg_id.between(n,n),
                            #                                                result_format='rows')[0].weight      
                            #self._Cij[step][statei,statej] += statej_weight/sourceweight
                            self._Cij[step][statei,statej] += self._weights[n]/sourceweight
        
        return self._Cij
        
    def getNumSegments(self):
        """docstring for getNumSegments"""
        return self._G.order()
        
    def _construct_graph(self,iterlist):
        """docstring for _construct_graph"""
        
        G = nx.DiGraph()
        print 'Building Graph'
        imin = int(min(iterlist) - 2)
        if imin < 1:
            imin = 1 #minimum is 1st iteration
            
        imax = int(min(max(iterlist) + max(self._lag) + 2,self._maxiter-1))
        
        # This returns a list of (parent_id, child_id) pairs; merges are not included
        conn_pairs = self._data_manager.get_connectivity(imin,n_iter_upper = imax)
        
        G.add_edges_from(conn_pairs)
        # Clean graph of any extra segments that do not have assignments
        maxnode = max(G.nodes()) + 1
        G.remove_nodes_from(xrange(self._maxsegid+1,maxnode))
        
        print 'Graph Completed'    
        return G
        
    def _retrieve_weights(self):
        """docstring for _retrieve_weights"""      
        
        # Extract weights from segments. Do it in pieces to avoid bringing
        # all of the segments objects into memory simultaneously
        weights = numpy.zeros((self._maxsegid+1,))
        
        for k in xrange(self._miniter,self._maxiter+1):
            segs = self._data_manager.get_segments(k)
            wk = numpy.array([x.weight for x in segs])
            seg_ids = sorted([x.seg_id for x in segs])
            imin = seg_ids[0]
            imax = seg_ids[-1]
            weights[imin:imax+1] = wk
            
        return weights
        
    def _generate_assignments_from_pcoord(self,runtime_config,pcoord=None):
        """docstring for _generate_assignments_from_pcoord"""

        import cPickle 
        
        # Determine Bin Boundaries from state pickle
        f = open(runtime_config['data.state'],'rb')
        stdict = cPickle.load(f)
        f.close()
        
        bin_bounds = stdict['we_driver'].bin_boundaries
        ndim = len(bin_bounds)
        print 'ndim: ',ndim
        print 'bin_bounds: ',bin_bounds
        nbins = []
        
        for d in range(ndim):
            nb = bin_bounds[d].size
            nbins.append(nb)
                   
        # Extract pcoord for all segments. Do it in pieces to avoid bringing lots of Segment objects
        # into memory simultaneously
        if pcoord == None:
            pcoord = numpy.zeros((self._maxsegid+1,ndim))
            for k in xrange(self._miniter,self._maxiter+1):
                if k % 20 == 0:
                    print k
                segs = self._data_manager.get_segments(k)
                pk = numpy.array([x.pcoord[-1] for x in segs]).reshape(-1,ndim) #only use last pcoord
                seg_ids = sorted([x.seg_id for x in segs])

                imin = seg_ids[0]
                imax = seg_ids[-1]

                print "imin %r imax %r pcoord[] %r pk %r" %(imin,imax,pcoord[imin:imax+1,:].shape,pk.shape)
                pcoord[imin:imax+1,:] = pk
                                
            import cPickle
            f = open(runtime_config['data.SIM_NAME'] + '.pcoords.pkl','wb')
            cPickle.dump(pcoord,f,cPickle.HIGHEST_PROTOCOL)
            f.close() 
        
        # Assign segments to bins
        bassign = ndim*[None]
        for idim in xrange(ndim):
            bassign[idim] = numpy.digitize(pcoord[:,idim],bin_bounds[idim]) - 1
            
        #if ndim == 1:
        #    linbassign = bassign[0]
        #elif ndim == 2:
        #    linbassign = bassign[1] + bassign[0]*(nbins[1]-1)   
        
        print "nbins:%r"%nbins
        
        linbassign = bassign[ndim - 1]            
        for idim in xrange(1,ndim):
            binsize = 1
            for jdim in xrange(0,idim):
                binsize *= (nbins[ndim - 1 - jdim] - 1)
                
            linbassign += bassign[ndim - 1 - idim] * binsize

        return linbassign 
        
    def _generate_assignments_from_pcoord_2d(self,runtime_config,bb2):
        """docstring for _generate_assignments_from_pcoord_2d"""    
        import cPickle    
        
        # Determine Bin Boundaries from state pickle
        f = open(runtime_config['data.state'],'rb')
        stdict = cPickle.load(f)
        f.close()
        
        bin_bounds = stdict['we_driver'].bin_boundaries
        ndim = len(bin_bounds)
        nbins = []
        
        for d in range(ndim):
            nb = bin_bounds[d].size
            nbins.append(nb)
        nbins.append(bb2.size)
        
        # Extract pcoord for all segments. Do it in pieces to avoid bringing lots of Segment objects
        # into memory simultaneously
        
        pcoord = numpy.zeros((self._maxsegid+1,ndim+1))
        for k in xrange(self._miniter,self._maxiter+1):
            segs = self._data_manager.get_segments(k)
            pk = numpy.array([x.pcoord[-1] for x in segs]).reshape(-1,)
            seg_ids = sorted([x.seg_id for x in segs])
            imin = seg_ids[0]
            imax = seg_ids[-1]
            pcoord[imin:imax+1,0] = pk
                                          
        # Assign segments to bins
        bassign = (ndim+1)*[None]
        for idim in xrange(ndim):
            bassign[idim] = numpy.digitize(pcoord[:,idim],bin_bounds[idim]) - 1
        
        # Find reference structure in each bin
        maxbin = numpy.max(bassign[0])
        first_segid = {}
        for k in xrange(maxbin+1):
            #print k, numpy.where(assignments == k)[0].shape
            dd = int(0.75*numpy.where(bassign[0] == k)[0].shape[0])
            #print dd
            first_segid[k] = numpy.where(bassign[0] == k)[0][dd]
          
        # Get coordinates of reference structure in each bin
        print 'Retreiving reference structure coordinates'
        import netCDF4 as netcdf
        from scipy.spatial.distance import pdist
        nc = netcdf.Dataset('segfile.nc','r')
        print 'Finished opening netcdf file'
        x = nc.variables['coords'] 
        binRefCoord = {}
        binRefPdist = {}
        for k in xrange(maxbin+1):
            binRefCoord[k] = x[first_segid[k],:,:]
            binRefPdist[k] = pdist(binRefCoord[k],'euclidean')
         
        # Loop through all coordinates and measure the dRMSD between coords and binRefCoord
        print 'Calculating dRMSD'
        prefac = 1.0/binRefPdist[0].shape[0]

        #drmsd = {}
        segids = {}
        for k in xrange(maxbin+1):
            print 'Bin %d' % (k)
            #drmsd[k] = []
            ii = numpy.where(bassign[0] == k)[0]
            segids[k] = ii
            #print ii
            for jj in ii:
                #print jj
                curpd = pdist(x[jj,:,:],'euclidean')
                pcoord[jj,1] = -10.0*numpy.sqrt(prefac*numpy.sum((curpd - binRefPdist[k])**2))
        
        nc.close()
 
        
        bassign[1] = numpy.digitize(pcoord[:,1],bb2) - 1
        
        f = open('assignlist.pkl','wb')
        cPickle.dump(bassign,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
            
        linbassign = bassign[0]*nbins[1] + bassign[1]   
        ii = numpy.where(bassign[0] == maxbin)[0]
        linbassign[ii] = bassign[0][ii]*nbins[1]
                
        return linbassign
        
    
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
        nsolve = nstates - ntarget              # Number of active states
    
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
    
    