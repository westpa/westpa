# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.



import numpy

class EDF:
    '''A class for creating and manipulating empirical distribution functions (cumulative 
    distribution functions derived from sample data).    
    '''
    
    @staticmethod
    def from_array(array):
        edf = EDF(None,None)
        edf.x = array[:,0]
        edf.F = array[:,1]
        edf.dF = numpy.diff(edf.F)
        return edf
        
    @staticmethod
    def from_arrays(x, F):
        edf = EDF(None,None)
        edf.x = x
        edf.F = F
        edf.dF = numpy.diff(edf.F)
        return edf
    
    def __init__(self, values, weights = None):
        '''Construct a new EDF from the given values and (optionally) weights.'''
        
        if values is None:
            self.x = None
            self.F = None
            self.dF = None
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
        x = values[numpy.concatenate(([True], values[1:] != values[:-1]))]
        F = numpy.empty((len(x),), numpy.float64)
        
        # ``values`` is arranged in increasing order, so we can walk along it and add up weights
        # as we go 
        ival_last = 0
        ival = 0
        for ibin in range(0, len(x)):
            while ival < len(values) and values[ival] <= x[ibin]:
                ival+=1
            F[ibin] = weights[ival_last:ival].sum()
            ival_last = ival
        F = numpy.add.accumulate(F)
        F /= F[-1]
        
        self.x = x
        self.F = F
        self.dF = numpy.diff(F)
        
    def __len__(self):
        return len(self.x)
                
    def __call__(self, x):
        '''Evaluate this EDF at the given abcissae.'''
        indices = numpy.digitize(x, self.x)
        indices[indices >= len(self.x)] = len(self.x) - 1
        return self.F[indices]

    
    def as_array(self):
        '''Return this EDF as a (N,2) array, where N is the number of unique values passed to
        the constructor.  Numpy type casting rules are applied (so, for instance, integral abcissae
        are converted to floating-point values).'''
        
        result = numpy.empty((len(self.F),2), dtype=numpy.result_type(self.x, self.F))
        result[:,0] = self.x
        result[:,1] = self.F
        return result
    
    def quantiles(self, p):
        '''Treating the EDF as a quantile function, return the values of the (statistical) variable whose
        probabilities are at least p.  That is, Q(p) = inf {x: p <= F(x) }.'''
        
        indices = numpy.searchsorted(self.F, p)
        indices[indices >= len(self.x)] = len(self.x) - 1
        return self.x[indices]
    
    def quantile(self, p):
        return self.quantiles([p])[0]
    
    def median(self):
        return self.quantiles([0.5])[0]
    
    def moment(self, n):
        '''Calculate the nth moment of this probability distribution
        
        <x^n> = int_{-inf}^{inf} x^n dF(x)
        '''
        
        if n == 1:
            return (self.x[:-1] * self.dF).sum()
        else:
            return (self.x[:-1]**n * self.dF).sum()
        
    def cmoment(self, n):
        '''Calculate the nth central moment of this probability distribution'''
        
        if n < 2:
            return 0
        return ((self.x[:-1]-self.moment(1))**n * self.dF).sum()
    
    def mean(self): 
        return self.moment(1)
    
    def var(self):
        '''Return the second central moment of this probability distribution.'''
        return self.cmoment(2)
    
    def std(self):
        '''Return the standard deviation (root of the variance) of this probability distribution.'''
        return self.cmoment(2)**0.5

    