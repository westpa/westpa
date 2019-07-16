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

import logging
log = logging.getLogger(__name__)

EPS = numpy.finfo(numpy.float64).eps

class Bin(set):
    def __init__(self, iterable=None, label=None):
        super(Bin,self).__init__(iterable or [])
        self.label = label
    
    def __repr__(self):
        return ('<{classname} at 0x{id:x}, label={label!r}, count={count:d}, weight={weight:g}>'
                .format(classname=self.__class__.__name__, id=id(self), label=self.label,
                        count=len(self), weight=self.weight))
    
    @property
    def weight(self):
        'Total weight of all walkers in this bin'

        weight = 0.0
        for particle in self:
            weight += particle.weight
        return weight

    def reweight(self, new_weight):
        """Reweight all walkers in this bin so that the total weight is new_weight"""

        if len(self) == 0 and new_weight == 0:
            return
        
        if len(self) == 0 and new_weight != 0:
            raise ValueError('cannot reweight empty ParticleCollection')         
        
        current_weight = self.weight
        log.debug('reweighting collection of {:d} particles from {:g} to {:g}'.format(len(self), current_weight, new_weight))
        assert (new_weight == 0 and current_weight == 0) or new_weight > 0
                
        wrat = new_weight / current_weight
        for p in self:
            p.weight *= wrat
            
        log.debug('new weight: {:g}'.format(self.weight))
        assert abs(new_weight-self.weight) <= EPS*len(self)
    
    


    
    