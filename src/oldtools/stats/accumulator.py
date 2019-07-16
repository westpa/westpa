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
NAN=float('nan')
class RunningStatsAccumulator:
    def __init__(self, shape, dtype=numpy.float64, count_dtype=numpy.uint, weight_dtype=numpy.float64, mask_value=NAN):
        self.sum = numpy.zeros(shape, dtype)
        self.sqsum = numpy.zeros(shape, dtype)
        self.weight = numpy.zeros(shape, weight_dtype)
        self.count  = numpy.zeros(shape, count_dtype)
        self.mask_value = mask_value
        
    def incorporate(self, index, value, weight):
        self.count[index] += 1
        self.weight[index] += weight
        self.sum[index] += weight*value
        self.sqsum[index] += weight*value*value
        
    def average(self):
        valid = (self.count > 0)
        avg = numpy.empty_like(self.sum)
        avg[valid] = self.sum[valid] / self.weight[valid]
        avg[~valid] = self.mask_value
        return avg
    mean = average
    
    def std(self):
        valid = (self.count > 0)
        vavg = self.average()[valid]
        std = numpy.empty_like(self.sqsum)
        std[valid] = (self.sqsum[valid] / self.weight[valid] - vavg*vavg)**0.5
        std[~valid] = self.mask_value
        return std
