import numpy
from itertools import izip

class TransitionEventFinder(object):
    def __init__(self):
        self.indices_by_region = {}
        
        self.crossings = None
        self.ordered_crossings = None
        self.transition_indices = None
        self.passage_indices = None
        
        self.npoints = None
        
    def assign_regions(self, Q):
        raise NotImplementedError
        
    def find_crossings(self):
        ordered_crossings = []
        crossing_indices = {}
        
        regions = (list(self.regions) 
                   + list(reversed(self.regions[:-1]))) 
        boundaries = list(izip(regions[:-1], regions[1:])) 
        
        for i in xrange(0, self.npoints-1):
            for (region1, region2) in boundaries:
                sr1 = self.indices_by_region[region1]
                sr2 = self.indices_by_region[region2]            
                if i in sr1 and i+1 in sr2: 
                    ordered_crossings.append((i, region1, region2))
        self.ordered_crossings = ordered_crossings

    
    def find_transitions(self):
        if not self.ordered_crossings:
            self.find_crossings()
        
        transitions = zip(self.regions[:-2], self.regions[2:])
        transitions += list(tuple(reversed(item)) for item in transitions) 
            
        wait_starts = {}
        transition_indices = {}
        passage_indices = {}
        
        for ((idx1, region1, junk1), (idx2, junk2, region2)) \
        in izip(self.ordered_crossings[:-1], self.ordered_crossings[1:]):
            if (region1, region2) in transitions:
                # Record the transition from region1 to region2
                try:
                    transition_indices[region1, region2].append((idx1, idx2))
                except KeyError:
                    transition_indices[region1, region2] = [(idx1, idx2)]
                
                try:
                    wait_start = wait_starts[region1, region2]
                except KeyError:
                    wait_starts[region1, region2] = idx2
                else:
                    del wait_starts[region1, region2]
                    try:
                        passage_indices[region1, region2].append((wait_start, idx2))
                    except KeyError:
                        passage_indices[region1, region2] = [(wait_start, idx2)]
        
        self.transition_indices = dict((k, numpy.array(v, numpy.uintp))
                                       for (k,v) in transition_indices.iteritems())
        self.passage_indices = dict((k, numpy.array(v, numpy.uintp))
                                    for (k,v) in passage_indices.iteritems())
        
class PC1DCrossingFinder(TransitionEventFinder):
    regions = ('A', 'T', 'B')
    def __init__(self, Q_trans_lb, Q_trans_ub):
        super(PC1DCrossingFinder,self).__init__()
        self.Q_trans_lb = Q_trans_lb
        self.Q_trans_ub = Q_trans_ub
        
    def assign_regions(self, Q):
        self.npoints = Q.shape[0]
        
        indices = numpy.fromiter(xrange(0, self.npoints), numpy.uintp)
        self.indices_by_region['A'] = frozenset(indices[Q < self.Q_trans_lb])
        self.indices_by_region['T'] = frozenset(indices[ (Q >= self.Q_trans_lb)
                                                        &(Q <= self.Q_trans_ub)])
        self.indices_by_region['B'] = frozenset(indices[Q > self.Q_trans_ub])
