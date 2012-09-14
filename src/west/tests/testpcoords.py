from __future__ import division; __metaclass__ = type

import numpy
import west
import west.pcoords
import nose.tools

inf = float('inf')

class TestParticleSet:    
    def test_get_all_bins(self):
        pset = west.pcoords.ParticleSet()
        pset.add(west.Segment(weight=0.1))
        pset.add(west.Segment(weight=0.1))        
        assert list(pset.get_all_bins()) == [pset]
        
    def test_weight(self):
        pset = west.pcoords.ParticleSet()
        pset.add(west.Segment(weight=0.1))
        pset.add(west.Segment(weight=0.1))
        nose.tools.assert_almost_equal(pset.weight, 0.2)
        
    def test_reweight(self):
        pset = west.pcoords.ParticleSet()
        pset.add(west.Segment(weight=0.1))
        pset.add(west.Segment(weight=0.1))
        pset.reweight(1.0)
        nose.tools.assert_almost_equal(pset.weight, 1.0)
        
class TestParticleList:
    def test_clear(self):
        plist = west.pcoords.ParticleList()
        plist.append(west.Segment(weight=0.1))
        plist.append(west.Segment(weight=0.1))
        plist.clear()
        assert len(plist) == 0
        
class TestRectilinearRegionSet:
    def test_1d_index_mapping(self):
        boundaries = [[-inf] + [x for x in numpy.arange(-1, 0, 0.2)] + [x for x in numpy.arange(0, 1.2, 0.2)] + [inf]]
        coords = [[x] for x in (-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)]
        rrs = west.pcoords.RectilinearRegionSet(boundaries)
    
        """
        bin structure: (a,b) => x in [a,b)
        0: [(-inf, -1.0)
        1: [(-1.0, -0.80000000000000004)]
        2: [(-0.80000000000000004, -0.60000000000000009)]
        3: [(-0.60000000000000009, -0.40000000000000013)]
        4: [(-0.40000000000000013, -0.20000000000000018)]
        5: [(-0.20000000000000018, 0.0)]
        6: [(0.0, 0.20000000000000001)]
        7: [(0.20000000000000001, 0.40000000000000002)]
        8: [(0.40000000000000002, 0.60000000000000009)]
        9: [(0.60000000000000009, 0.80000000000000004)]
        10:[(0.80000000000000004, 1.0)]
        11:[(1.0, inf)]
        """
        indices = rrs._map_to_indices(numpy.array(coords))
        assert list(indices) == [1,2,3,4,6,7,8,9,11]
        
    def test_1d_bin_mapping(self):
        boundaries = [[-inf] + [x for x in numpy.arange(-1, 0, 0.2)] + [x for x in numpy.arange(0, 1.2, 0.2)] + [inf]]
        coords = [[x] for x in (-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)]
        rrs = west.pcoords.RectilinearRegionSet(boundaries)
        bins = rrs.get_all_bins()
        mapped_bins = rrs.map_to_bins(coords)
        assert list(mapped_bins) == [bins[1], bins[2], bins[3], bins[4], bins[6], bins[7], bins[8], bins[9], bins[11]]

    @nose.tools.raises(ValueError)        
    def test_1d_bin_out_of_space(self):
        boundaries = [[x for x in numpy.arange(-1, 0, 0.2)] + [x for x in numpy.arange(0, 1.2, 0.2)] + [inf]]
        rrs = west.pcoords.RectilinearRegionSet(boundaries)
        rrs._map_to_indices(numpy.array([[-2]]))

    def test_2d_index_mapping(self):
        boundaries = [(-1,-0.5,0,0.5,1), (-1,-0.5,0,0.5,1)]
        coords = [ (-0.75, -0.75), (-0.25,-0.25), (0,0), (0.25,0.25), (0.75,0.75), (-0.25, 0.75), (0.25,-0.75)]
        rrs = west.pcoords.RectilinearRegionSet(boundaries)
        
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
        
        indices = rrs._map_to_indices(numpy.array(coords))
        assert list(indices) == [0, 5, 10, 10, 15, 7, 8]
        
    def test_2d_bin_mapping(self):
        boundaries = [(-1,-0.5,0,0.5,1), (-1,-0.5,0,0.5,1)]
        coords = [ (-0.75, -0.75), (-0.25,-0.25), (0,0), (0.25,0.25), (0.75,0.75), (-0.25, 0.75), (0.25,-0.75)]
        rrs = west.pcoords.RectilinearRegionSet(boundaries)
        bins = rrs.get_all_bins()
        mapped_bins = rrs.map_to_bins(coords)
        assert list(mapped_bins) == [bins[0], bins[5], bins[10], bins[10], bins[15], bins[7], bins[8]]

class TestPiecewiseRegionSet:
    def test_bin_mapping(self):
        coords = [[-0.5], [0.0], [0.5]]
        fr1 = (lambda x: x<0)
        fr2 = (lambda x: x>=0)
        prs = west.pcoords.PiecewiseRegionSet([fr1,fr2], [float('-inf'),float('inf')], 1)
        bins = prs.get_all_bins()
        mapped_bins = prs.map_to_bins(coords)
        assert list(mapped_bins) == [bins[0], bins[1], bins[1]]
        
class TestRegionSet:
    def test_replace_region(self):
        boundaries_outer = [[-1,0,1]]
        boundaries_inner = [[-1, -0.5, 0]]
        rrso = west.pcoords.RectilinearRegionSet(boundaries_outer)
        rrsi = west.pcoords.RectilinearRegionSet(boundaries_inner)
        rrso.replace_region_containing([-0.75], rrsi)
        assert rrsi in rrso.regions
        assert rrso.get_bin_containing([-0.75]) is rrsi.get_bin_containing([-0.75])
        