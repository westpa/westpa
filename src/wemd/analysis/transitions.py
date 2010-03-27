__metaclass__ = type

from copy import copy
import numpy

class RegionSet:
    def __init__(self, names, boundaries, adjacency):
        self.names = names
        self.boundaries = boundaries
        self.adjacency = numpy.array(adjacency, numpy.bool_)
        if self.adjacency.shape != (len(self.names), len(self.names)):
            raise TypeError('invalid adjacency matrix supplied')
        
        self.ids_by_name = dict((name, i) for (i, name) in enumerate(names))
        
    def localize(self, q):
        raise NotImplementedError
    
    def localize_all(self, qs):
        return [self.localize(q) for q in qs]
    
    def nonadjacent_pairs(self):
        return numpy.argwhere(~self.adjacency)
        
class OneDimRegionSet(RegionSet):
    def __init__(self, names, boundaries, adjacency = None):
        if adjacency is None:
            adjacency = numpy.ones((len(names),len(names)), numpy.bool_)
            for ix in xrange(0,len(names)):
                for iy in xrange(0, len(names)):
                    if abs(ix-iy) > 1:
                        adjacency[ix,iy] = 0
        super(OneDimRegionSet,self).__init__(names, boundaries, adjacency)
                
    def localize(self, q):
        for (i, (lb,ub)) in enumerate(self.boundaries):
            if lb <= q < ub:
                return i
        else:
            raise ValueError('coordinate is not in any defined region')
    
    def localize_all(self, qs):
        regions = numpy.zeros((len(qs),), numpy.uintp)
        for (ireg, (lb, ub)) in enumerate(self.boundaries):
            regions[(qs[:,0] >= lb) & (qs[:,0] < ub)] = ireg
        return regions
    
class TransitionEventAccumulator:
    def __init__(self, regions,
                 timestep = 1.0,
                 data_overlaps = False, 
                 accumulate_eds = True,
                 accumulate_fpts = True,
                 transition_log = None):
        
        self.regions = regions
        self.timestep = timestep
        self.transition_log = transition_log
        self.data_overlaps = data_overlaps
        
        if accumulate_eds is True:
            self.accumulate_eds = set(tuple(pair) for pair in regions.nonadjacent_pairs())
        elif not accumulate_eds:
            self.accumulate_eds = set()
        else:
            self.accumulate_eds = set(tuple(pair) for pair in accumulate_eds)
        
        if accumulate_fpts is True:
            self.accumulate_fpts = set(tuple(pair) for pair in regions.nonadjacent_pairs())
        elif not accumulate_fpts:
            self.accumulate_fpts = set()
        else:
            self.accumulate_fpts = set(tuple(pair) for pair in accumulate_fpts)

        self.eds = dict((key, []) for key in self.accumulate_eds)
        self.fpts = dict((key, []) for key in self.accumulate_fpts)

        self.time_index = 0
        self.last_iregion = None
        
        self.completion_indices = numpy.zeros(self.regions.adjacency.shape, numpy.uint64)
        self.event_counts = numpy.zeros(self.regions.adjacency.shape, numpy.uint64)

    def get_state(self):
        return {'time_index': self.time_index,
                'last_iregion': self.last_iregion,
                'completion_indices': copy(self.completion_indices)}
    
    def set_state(self, state):
        self.time_index = state['time_index']
        self.last_iregion = state['last_iregion']
        self.completion_indices = state['completion_indices']
    
    def identify_transitions(self, pcoords, weights):
        pcoord_regions = self.regions.localize_all(pcoords)
        
        if len(pcoords) == 1:
            if self.last_iregion is None:
                self.last_iregion = pcoord_regions[0]
                self.time_index += 1
                return
            else:
                ti_start = 0
        else:
            if self.data_overlaps or self.last_iregion is None:
                self.last_iregion = pcoord_regions[0]
                ti_start = 1
            else:
                ti_start = 0

        for ti in xrange(ti_start, len(pcoords)):
            self.time_index += 1
            iregion = pcoord_regions[ti]
            q = pcoords[ti]
            w = weights[ti]
            
            if iregion != self.last_iregion:
                # A crossing has occurred
                self.completion_indices[self.last_iregion, iregion] = self.time_index
                self.event_counts[self.last_iregion, iregion] += 1
                
                if self.transition_log:
                    self.transition_log.write('%-12s    %21.16g    %21.16g    %s->%s\n'
                                              % ('', self.time_index * self.timestep, 
                                                 w, 
                                                 self.regions.names[self.last_iregion],
                                                 self.regions.names[iregion]))
                
                sources_nonadjacent = numpy.argwhere(~self.regions.adjacency[:,iregion])
                for isrc in sources_nonadjacent:
                    isrc = int(isrc)
                    if self.completion_indices[isrc,iregion] < self.completion_indices[isrc, self.last_iregion]:
                        #print "found %s->%s transition at time %s (last %s)" \
                        #      % (self.regions.names[isrc], self.regions.names[iregion],
                        #         self.time_index * self.timestep,
                        #         self.completion_indices[isrc,iregion]*self.timestep)
                        #print "adjacency matrix:"
                        #print self.regions.adjacency
                        #print "completion indices (not yet updated):"
                        #print self.completion_indices
                        
                        if (isrc, iregion) in self.accumulate_eds and self.completion_indices[isrc,self.last_iregion]>0:
                            self.eds[isrc,iregion].append((self.time_index - self.completion_indices[isrc, self.last_iregion] - 1, w))
                            #print "found event duration %f" % ((self.time_index - self.completion_indices[isrc, self.last_iregion] - 1) * self.timestep)
                        if (isrc, iregion) in self.accumulate_fpts and self.completion_indices[iregion, isrc] > 0:
                            # completed first passage
    
                            self.fpts[isrc,iregion].append((self.time_index - self.completion_indices[iregion,isrc], w))

                        self.completion_indices[isrc,iregion] = self.time_index
                        self.event_counts[isrc,iregion] += 1
                    else:
                        # fluctuation across boundary without transition
                        pass
            
            self.last_iregion = iregion
