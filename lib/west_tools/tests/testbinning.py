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

from __future__ import division, print_function
from westpa.binning.assign import (RectilinearBinMapper, PiecewiseBinMapper, FuncBinMapper, VectorizingFuncBinMapper, 
                                 VoronoiBinMapper, RecursiveBinMapper)
from westpa.binning.assign import index_dtype, coord_dtype
from westpa.binning._assign import testfunc #@UnresolvedImport

import numpy
from scipy.spatial.distance import cdist
import nose
import nose.tools

import logging
log = logging.getLogger(__name__)

class TestRectilinearBinMapper:
    def test1dAssign(self):
        bounds = [0.0, 1.0, 2.0, 3.0]
        coords = numpy.array([0, 0.5, 1.5, 1.6, 2.0, 2.0, 2.9])[:,None]
        
        assigner = RectilinearBinMapper([bounds])
        assert (assigner.assign(coords) == [0, 0, 1, 1, 2, 2, 2]).all()
        
    def test2dAssign(self):
        boundaries = [(-1,-0.5,0,0.5,1), (-1,-0.5,0,0.5,1)]
        coords = numpy.array([ (-0.75, -0.75), (-0.25,-0.25), (0,0), (0.25,0.25), (0.75,0.75), (-0.25, 0.75), (0.25,-0.75)])
        assigner = RectilinearBinMapper(boundaries)
        
        """bin structure: [(a,b), (c,d)] => x in [a,b), y in [c, d)
        0:[(-1, -0.5), (-1, -0.5)]
        1:[(-1, -0.5), (-0.5, 0)]
        2:[(-1, -0.5), (0, 0.5)]
        3:[(-1, -0.5), (0.5, 1)]
        4:[(-0.5, 0), (-1, -0.5)]
        5:[(-0.5, 0), (-0.5, 0)]
        6:[(-0.5, 0), (0, 0.5)]
        7:[(-0.5, 0), (0.5, 1)]
        8:[(0, 0.5), (-1, -0.5)]
        9:[(0, 0.5), (-0.5, 0)]
        10:[(0, 0.5), (0, 0.5)]
        11:[(0, 0.5), (0.5, 1)]
        12:[(0.5, 1), (-1, -0.5)]
        13:[(0.5, 1), (-0.5, 0)]
        14:[(0.5, 1), (0, 0.5)]
        15:[(0.5, 1), (0.5, 1)]"""
        
        assert (assigner.assign(coords) == [0, 5, 10, 10, 15, 7, 8]).all()

class TestPiecewiseBinMapper:
    def test_bin_mapping(self):
        coords = numpy.array([[-0.5], [0.0], [0.5]], dtype=numpy.float32)
        fr1 = (lambda x: x<0)
        fr2 = (lambda x: x>=0)
        pm = PiecewiseBinMapper([fr1,fr2])
        assert list(pm.assign(coords)) == [0, 1, 1]
        
class TestFuncBinMapper:
    @staticmethod
    def fn(coords, mask, output):
        output[mask & (coords[:,0] <= 0.5)] = 0
        output[mask & (coords[:,0] > 0.5)] = 1
        
    def test_fmapper(self):
        mapper = FuncBinMapper(self.fn, 2)
        coords = numpy.array([0.0, 0.1, 0.5, 0.7])
        coords.shape = (coords.shape[0], 1)
        print(repr(coords))
        output = mapper.assign(coords)
        print(repr(output))
        assert list(output) == [0,0,0,1]
    

class TestVectorizingFuncBinMapper:
    @staticmethod
    def fn(coord):
        if coord[0] <= 0.5:
            return 0
        else:
            return 1
        
    def test_vmapper(self):
        mapper = VectorizingFuncBinMapper(self.fn, 2)
        coords = numpy.array([0.0, 0.1, 0.5, 0.7])
        coords.shape = (coords.shape[0], 1)
        output = mapper.assign(coords)
        assert list(output) == [0,0,0,1]
        
class TestVoronoiBinMapper:
    @staticmethod
    def distfunc(coordvec, centers):
        if coordvec.ndim < 2:
            new_coordvec = numpy.empty((1,coordvec.shape[0]), dtype=coord_dtype)
            new_coordvec[0,:] = coordvec[:]
            coordvec = new_coordvec 
        distmat = numpy.require(cdist(coordvec, centers), dtype=coord_dtype)
        return distmat[0,:]
        
    def test_vmapper(self):
        centers = numpy.array([[0,0], [2,2]], dtype=coord_dtype)
        coords = numpy.array([[0,0], [2,2], [0.9,0.9], [1.1, 1.1]], dtype=coord_dtype)
        
        mapper = VoronoiBinMapper(self.distfunc, centers)
        output = mapper.assign(coords)
        assert list(output) == [0,1,0,1]
        
 
class TestNestingBinMapper:
    #pass
    '''         
         0                            1                      2
         +----------------------------+----------------------+
         |            0.5             |         1.5          |
         | +-----------+------------+ | +--------+---------+ |
         | |    0.25   |            | | |        |         | |
         | | +---+---+ |            | | |        |         | |
         | | |   |   | |            | | |        |         | |
         | | +---+---+ |            | | |        |         | |
         | +-----------+------------+ | +--------+---------+ |
         +---------------------------------------------------+    
    '''
    
    @staticmethod
    def fn1(coords, mask, output):
        test = coords[:,0] < 1
        output[mask & test] = 0
        output[mask & ~test] = 1
    
    @staticmethod
    def fn2(coords, mask, output):
        test = coords[:,0] < 0.5
        output[mask & test] = 0
        output[mask & ~test] = 1

    @staticmethod        
    def fn3(coords, mask, output):
        test = coords[:,0] < 0.25
        output[mask & test] = 0
        output[mask & ~test] = 1
    
    @staticmethod
    def fn4(coords, mask, output):
        test = coords[:,0] < 1.5
        output[mask & test] = 0
        output[mask & ~test] = 1
        
    def testOuterMapper(self):
        #pass
        '''         
             0                            1                      2
             +----------------------------+----------------------+
             |              0             |          1           |
             +---------------------------------------------------+    
        '''
        
        mapper = FuncBinMapper(self.fn1,2)
        rmapper = RecursiveBinMapper(mapper)
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [0,0,0,0,0,1]
        
    def testOuterMapperWithOffset(self):
        #pass
        '''         
             0                            1                      2
             +----------------------------+----------------------+
             |              1             |          2           |
             +---------------------------------------------------+    
        '''

        mapper = FuncBinMapper(self.fn1,2)
        rmapper = RecursiveBinMapper(mapper,1)
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [1,1,1,1,1,2]
        
    def testSingleRecursion(self):
        #pass
    
        '''         
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |                      |
             | +-----------+------------+ |                      |
             | |           |            | |                      |
             | |     1     |     2      | |          0           |
             | |           |            | |                      |
             | |           |            | |                      |
             | +-----------+------------+ |                      |
             +---------------------------------------------------+    
        '''

        outer_mapper = FuncBinMapper(self.fn1,2)
        inner_mapper = FuncBinMapper(self.fn2,2)
        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(inner_mapper, [0.5])
        assert rmapper.nbins == 3
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [1, 1, 1, 1, 2, 0]
        
    def testDeepRecursion(self):
        #pass
    
        '''
         0                            1                      2
         +----------------------------+----------------------+
         |            0.5             |                      |
         | +-----------+------------+ |                      |
         | |    0.25   |            | |                      |
         | | +---+---+ |            | |           0          |
         | | | 2 | 3 | |     1      | |                      |
         | | +---+---+ |            | |                      |
         | +-----------+------------+ |                      |
         +---------------------------------------------------+    
        '''
        
        outer_mapper = FuncBinMapper(self.fn1,2)
        middle_mapper = FuncBinMapper(self.fn2,2)
        inner_mapper  = FuncBinMapper(self.fn3,2)
        
        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper, [0.5])
        rmapper.add_mapper(inner_mapper, [0.25])
        
        assert rmapper.nbins == 4
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1]])
        output = rmapper.assign(coords)
        assert list(output) == [2, 2, 3, 3, 1, 0]
        
    def testSideBySideRecursion(self):
        #pass
        '''         
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |         1.5          |
             | +-----------+------------+ | +--------+---------+ |
             | |           |            | | |        |         | |
             | |     0     |     1      | | |   2    |    3    | |
             | |           |            | | |        |         | |
             | |           |            | | |        |         | |
             | +-----------+------------+ | +--------+---------+ |
             +---------------------------------------------------+    
        '''
        outer_mapper = FuncBinMapper(self.fn1,2)
        middle_mapper2 = FuncBinMapper(self.fn2,2)
        middle_mapper1 = FuncBinMapper(self.fn4,2)

        
        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper1, [1.5])
        rmapper.add_mapper(middle_mapper2, [0.5])
        
        assert rmapper.nbins == 4
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        output = rmapper.assign(coords)
        assert list(output) == [0,0,0,0,1,2,3]
        

    def testMegaComplexRecursion(self):
        #pass
        '''         
             0                            1                      2
             +----------------------------+----------------------+
             |            0.5             |         1.5          |
             | +-----------+------------+ | +--------+---------+ |
             | |    0.25   |            | | |        |         | |
             | | +---+---+ |            | | |        |         | |
             | | | 1 | 2 | |     0      | | |   3    |    4    | |
             | | +---+---+ |            | | |        |         | |
             | +-----------+------------+ | +--------+---------+ |
             +---------------------------------------------------+    
        '''
        outer_mapper = FuncBinMapper(self.fn1,2)
        middle_mapper1 = FuncBinMapper(self.fn2,2)
        middle_mapper2 = FuncBinMapper(self.fn4,2)
        inner_mapper = FuncBinMapper(self.fn3,2)
        
        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(middle_mapper1, [0.5])
        rmapper.add_mapper(inner_mapper, [0.25])        
        rmapper.add_mapper(middle_mapper2, [1.5])
        
        
        assert rmapper.nbins == 5
        coords = numpy.array([[0.1], [0.2], [0.3], [0.4], [0.6], [1.1], [1.6]])
        output = rmapper.assign(coords)
        assert list(output) == [1,1,2,2,0,3,4]

    def test2dRectilinearRecursion(self):
        '''
             0                            1                      2
             +----------------------------+----------------------+
             |                            |         1.5          |
             |                            | +--------+---------+ |
             |                            | |        |         | |
             |             0              | |   4    |   5     | |
             |                            | |        |         | |
             |                            | |        |         | |
             |                            | +--------+---------+ |
            1+---------------------------------------------------+
             |            0.5             |                      |
             | +-----------+------------+ |                      |
             | |           |            | |                      |
             | |    2      |     3      | |           1          |
             | |           |            | |                      |
             | |           |            | |                      |
             | +-----------+------------+ |                      |
            2+---------------------------------------------------+

        '''

        outer_mapper = RectilinearBinMapper([[0,1,2],[0,1,2]])

        upper_right_mapper = RectilinearBinMapper([[1,1.5,2],[0,2]])
        lower_left_mapper = RectilinearBinMapper([[0,0.5,1], [0,2]])

        rmapper = RecursiveBinMapper(outer_mapper)
        rmapper.add_mapper(upper_right_mapper, [1.5,0.5])
        rmapper.add_mapper(lower_left_mapper, [0.5, 1.5])


        pairs = [(0.5, 0.5), (1.25, 0.5), (1.75, 0.5),
                 (0.25, 1.5), (0.75, 1.5), (1.5, 1.5)]

        assert rmapper.nbins == 6
        assert (rmapper.assign(pairs) == [0,4,5,2,3,1]).all()


