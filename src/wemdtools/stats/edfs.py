from __future__ import division, print_function; __metaclass__ = type

import numpy

class EDF:
    '''A class for creating and manipulating empirical distribution functions (cumulative 
    distribution functions derived from sample data).
    
    Note that, like ``numpy.histogram()``, all intervals are half-open except for the
    right-most.  That is, the maximum abcissa has an ordinate of 1.0.
    '''
    
    def __init__(self, values, weights = None):
        '''Construct a new EDF from the given values and (optionally) weights.'''
        
        if values = None:
            self.abcissae = None
            self.edf = None
            return
        
        if weights is None:
            weights = numpy.ones((len(values)), numpy.float64)
        elif numpy.isscalar(weights):
            tweights = numpy.empty((len(values)), numpy.float64)
            tweights[:] = weights
            weights = tweights
        else:
            if len(weights) != len(values):
                raise TypeError('values and weights have different lengths')
            
        # Sort values
        sort_indices = numpy.argsort(values)
        values = values[sort_indices]
        weights = weights[sort_indices]
        
        # Determine unique abcissae; this is essentially stolen from numpy.lib.arraysetops.unique()         
        abcissae = values[numpy.concatenate(([True], values[1:] != values[:-1]))]        
        (hist,binbounds) = numpy.histogram(values, bins=abcissae, weights=weights)
        edf = numpy.add.accumulate(hist)
        edf /= edf[-1]
        
        self.abcissae = abcissae
        self.edf = edf
                
    def __call__(self, abcissae):
        '''Evaluate this EDF at the given abcissae.'''
        results = numpy.empty((len(abcissae),), numpy.float64)
        indices = numpy.digitize(abcissae, self.abcissae)
        results[indices >= len(self.abcissae)] = 1.0
        results[indices < len(self.abcissae)] = self.edf[indices]
        return results
    
    def as_array(self):
        '''Return this EDF as a (N,2) array, where N is the number of unique values passed to
        the constructor.  Numpy type casting rules are applied (so, for instance, integral abcissae
        are converted to floating-point values).'''
        
        result = numpy.empty((len(self.edf),), dtype=numpy.result_type(self.abcissae, self.edf))
        result[:,0] = self.abcissae
        result[:,1] = self.edf
        return result
    
    
    
        
        