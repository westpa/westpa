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
from west.we_driver import WEDriver
from west.systems import WESTSystem
from westpa.binning import RectilinearBinMapper
from west.states import TargetState, InitialState
from west import Segment
import numpy

EPS = numpy.finfo(numpy.float64).eps

import nose
import nose.tools

class TestWEDriver:    
    def setup(self):
        system = WESTSystem()
        system.bin_mapper = RectilinearBinMapper([[0.0, 1.0, 2.0]])
        system.bin_target_counts = numpy.array([4,4])
        system.pcoord_len = 2
        self.we_driver = WEDriver(system=system)
        self.system = system
        self._seg_id = 0

    def segment(self, init_pcoord, final_pcoord, weight=1.0):
        segment= Segment(n_iter=1, seg_id=self._seg_id,
                         pcoord=self.system.new_pcoord_array(),
                         weight=weight)
        segment.pcoord[0] = init_pcoord
        segment.pcoord[1] = final_pcoord
        self._seg_id += 1
        return segment
        
    def teardown(self):
        self.we_driver.clear()
        del self.we_driver
        del self.system
        del self._seg_id
                
    def test_assign(self):
        segments = [self.segment(0.0, 1.5, weight=0.5),
                    self.segment(1.5, 0.5, weight=0.5)]
        self.we_driver.new_iteration()
        n_recycled = self.we_driver.assign(segments)
        assert n_recycled == 0
        assert len(self.we_driver.initial_binning[0]) == 1
        assert len(self.we_driver.initial_binning[1]) == 1
        assert len(self.we_driver.final_binning[0]) == 1
        assert len(self.we_driver.final_binning[1]) == 1
        assert (self.we_driver.flux_matrix == numpy.array([[0.0, 0.5], [0.5,0.0]])).all()
        
    def test_passthrough(self):
        segments = ([self.segment(0.0, 1.5, weight=0.125) for _i in xrange(4)]
                   +[self.segment(1.5, 0.5, weight=0.125) for _i in xrange(4)])
        segs_by_id = {segment.seg_id: segment for segment in segments}
        self.we_driver.new_iteration()
        self.we_driver.assign(segments)
        self.we_driver.construct_next()
        out_segs = set()
        for bin in self.we_driver.next_iter_binning:
            out_segs.update(bin)
        assert out_segs == set(self.we_driver.next_iter_segments)
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_segments) - 1.0) < 8*EPS
        
        for segment in self.we_driver.next_iter_segments:
            # has n_iter been advanced?
            assert segment.n_iter == 2
            
            # is weight correct?
            assert segment.weight == 0.125 #rigorous floating point comparison okay here because no math was done
            
            # is parent set correctly?
            assert segment.parent_id is not None
            
            # was pcoord set correctly?
            assert (segs_by_id[segment.parent_id].pcoord[-1] == segment.pcoord[0]).all()
            
            # was status set correctly?
            assert segment.status == Segment.SEG_STATUS_PREPARED

        # were parent endpoint types set correctly            
        for segment in segments:
            assert segment.endpoint_type == Segment.SEG_ENDPOINT_CONTINUES
            
        
    def test_split_by_weight(self):
        segments = [self.segment(1.5, 0.5, weight=0.25),
                    self.segment(0.0, 1.5, weight=0.75)]
        self.we_driver.new_iteration()
        self.we_driver.assign(segments)
        self.we_driver.construct_next()
        assert len(self.we_driver.next_iter_binning[0]) == 4
        assert len(self.we_driver.next_iter_binning[1]) == 4
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_binning[0]) - 0.25) < 4*EPS
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_binning[1]) - 0.75) < 4*EPS 
        assert numpy.allclose([seg.weight for seg in self.we_driver.next_iter_binning[0]],
                              [0.25/4.0 for _i in xrange(4)])
        assert numpy.allclose([seg.weight for seg in self.we_driver.next_iter_binning[1]],
                              [0.75/4.0 for _i in xrange(4)])
        
        for ibin in xrange(2):
            for segment in self.we_driver.next_iter_binning[ibin]:
                print(segment)
                assert segment.n_iter == 2
                assert segment.parent_id is not None 
                assert segment.parent_id == segments[ibin].seg_id
                assert segment.status == Segment.SEG_STATUS_PREPARED

    # this test will fail up to alpha of the time
    def test_merge_by_weight(self):
        selected_counts = {0: 0, 1: 0}
        alpha = 0.01
        nrounds = 1000
        from scipy.stats import binom
        # lower and upper bounds of 95% CI for selecting the segment with weight 1/3
        lb = binom.ppf(alpha/2.0, nrounds, 1.0/3.0)
        ub = binom.ppf(1.0-alpha/2.0, nrounds, 1.0/3.0)
        
        system = WESTSystem()
        system.bin_mapper = RectilinearBinMapper([[0.0, 1.0]])
        system.bin_target_counts = numpy.array([1])
        system.pcoord_len = 2
        self.we_driver = WEDriver(system=system)
        self.system = system
        self._seg_id = 0
        
        segments = [Segment(n_iter=1, seg_id=0, pcoord=numpy.array([[0],[0.25]], dtype=numpy.float32),weight=1.0/3.0),
                    Segment(n_iter=1, seg_id=1, pcoord=numpy.array([[0],[0.75]], dtype=numpy.float32),weight=2.0/3.0)]
        
        for _iround in xrange(nrounds):
            for segment in segments:
                segment.endpoint_type = Segment.SEG_ENDPOINT_UNSET
                
            self.we_driver.new_iteration()
            self.we_driver.assign(segments)
            self.we_driver.construct_next()
            
            assert len(self.we_driver.next_iter_binning[0]) == 1
            newseg = self.we_driver.next_iter_binning[0].pop()
            
            assert segments[newseg.parent_id].endpoint_type == Segment.SEG_ENDPOINT_CONTINUES
            assert segments[~newseg.parent_id].endpoint_type == Segment.SEG_ENDPOINT_MERGED
            
            selected_counts[newseg.parent_id] += 1
            
        print(selected_counts)
        assert lb <= selected_counts[0] <= ub, ('Incorrect proportion of histories selected.'
                                                'this is expected about {:%} of the time; retry test.'.format(alpha))
        
    def test_split_with_adjust(self):
        # this is a split followed by merge
        self.system.bin_target_counts = numpy.array([5,5])
        segments = [self.segment(1.5, 0.5, weight=0.125), self.segment(1.5, 0.5, weight=0.125),
                    self.segment(0.0, 1.5, weight=0.375), self.segment(0.0, 1.5, weight=0.375)]
        self.we_driver.new_iteration()
        self.we_driver.assign(segments)
        self.we_driver.construct_next()
        assert len(self.we_driver.next_iter_binning[0]) == 5
        assert len(self.we_driver.next_iter_binning[1]) == 5
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_binning[0]) - 0.25) < 5*EPS
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_binning[1]) - 0.75) < 5*EPS 
        
        for ibin in xrange(2):
            for segment in self.we_driver.next_iter_binning[ibin]:
                print(segment)
                assert segment.n_iter == 2
                assert segment.parent_id is not None 
                assert segment.status == Segment.SEG_STATUS_PREPARED

    def test_split_with_adjust_istates(self):
        # this is a split followed by merge, for segments which are initial states
        self.system.bin_target_counts = numpy.array([5,5])
        segments = [self.segment(1.5, 0.5, weight=0.125), self.segment(1.5, 0.5, weight=0.125),
                    self.segment(0.0, 1.5, weight=0.375), self.segment(0.0, 1.5, weight=0.375)]
        self.we_driver.new_iteration()
        self.we_driver._prep_we()
        self.we_driver.used_initial_states[-1] = None
        self.we_driver.used_initial_states[-2] = None
        
        for ibin,bin in enumerate(self.we_driver.next_iter_binning):
            pc = numpy.array([[0.5+ibin],[0.0]])
            for iseg in xrange(6):
                segment = Segment(n_iter=1, seg_id=None, weight=1.0/12.0,
                                  parent_id=-(ibin+1), pcoord=pc)
                bin.add(segment)
        
                    
        for ibin in xrange(len(self.we_driver.next_iter_binning)):
            # This will raise KeyError if initial state tracking is done improperly
            self.we_driver._adjust_count(ibin)

        assert len(self.we_driver.next_iter_binning[0]) == 5
        assert len(self.we_driver.next_iter_binning[1]) == 5

                
    def test_recycle(self):
        segments = [self.segment(0.0, 1.5, weight=0.5),
                    self.segment(0.0, 0.5, weight=0.5)]
        tstate = TargetState('recycle', [1.5], 0)
        istate = InitialState(0, 0, 0, pcoord=[0.0])
        
        self.we_driver.new_iteration(initial_states=[istate], target_states=[tstate])
        n_needed = self.we_driver.assign(segments)
        assert n_needed == 0
        self.we_driver.construct_next()
        
        n_recycled = len(list(self.we_driver.recycling_segments))
        assert n_recycled == 1
                
        assert len(self.we_driver.next_iter_binning[0]) == 4
        assert len(self.we_driver.next_iter_binning[1]) == 0
        assert abs(sum(seg.weight for seg in self.we_driver.next_iter_binning[0]) - 1.0) < 4*EPS 
        assert numpy.allclose([seg.weight for seg in self.we_driver.next_iter_binning[0]],
                              [0.25 for _i in xrange(4)])
        assert segments[0].endpoint_type == Segment.SEG_ENDPOINT_RECYCLED

        
    def test_multiple_merge(self):
        
        # This weight and count combination is known to trigger a split to 51
        # followed by a count adjustment to 50 (thanks to Josh Adelman)
        segment = self.segment(0.0, 0.5, weight=0.9999999999970001)
        
        # Initial state ID 0
        segment.parent_id = -1
        segment.wtg_parent_ids = set([-1])
        assert segment.initpoint_type == segment.SEG_INITPOINT_NEWTRAJ
        
        self.system.bin_target_counts = numpy.array([50,50])
        self.we_driver.new_iteration()
        self.we_driver.assign([segment])
        self.we_driver.construct_next()
        
        assert len(self.we_driver.next_iter_binning[0]) == 50

    def check_populate_initial(self, prob, target_counts):
        istate = InitialState(0, 0, 0, pcoord=[0.0])
        self.system.bin_target_counts = numpy.array([target_counts, target_counts])

        self.we_driver.populate_initial([istate], [prob], system=self.system)
        assert len(self.we_driver.next_iter_binning[0]) == target_counts

    @nose.SkipTest
    def test_populate_initial(self):
        for prob in [0.1, 1.0 / 3.0, 0.9999999999970001, 1.0]:
            for tcount in xrange(30, 60):
                yield self.check_populate_initial, prob, tcount

    # TODO: add test for seeding the flux matrix based on recycling 
    # TODO: add test for split after merge in adjust count
