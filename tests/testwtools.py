'''

Created 3/19/2013

author: nbr

Test suite for west_tools

'''

from __future__ import division, print_function; __metaclass__ = type
import nose.tools
from nose.tools import raises

import numpy

from w_assign import WAssign
from westpa.binning._assign import assign_and_label
from westpa.binning.assign import RectilinearBinMapper

class Test_W_Assign_Single:
    '''
    Tests w_assign over a single iteration, using random pcoords for 4 segments over 10 timepoints'''

    ##Random pcoords for the 4 segments, over 10 time points
    pcoords = numpy.array([[  8.6250176 ,   0.94418644,   3.45604335,   8.06916037,
                          5.42291169,   1.75481853,   1.67204737,   5.18241533,
                          3.99005367,   7.1339341 ],
                       [  4.89593824,   3.65634435,   2.95876746,   2.1180446 ,
                          3.13918444,   1.17418433,  12.54785261,   4.0411646 ,
                         14.90786241,  13.57244566],
                       [  1.42639175,   0.84988666,  14.47342808,   9.85372988,
                         11.53806799,   9.48601102,  12.49882247,  13.23377607,
                          8.48030094,   2.63759144],
                       [ 14.57179029,   6.31455601,   5.99651334,  10.49378797,
                          5.03735065,  13.99786996,   1.62226425,  11.01858827,
                         14.39726662,  13.57312496]])

    pcoords = pcoords.reshape(4,10,1)
    parent_ids = numpy.array([-1, -1, -1, -1])
    nsegs, npts = pcoords.shape[0:2]

    #Expected bin assignments for each segment
    seg1_bins = numpy.array([ 8,  0,  3,  8,  5,  1,  1,  5,  3,  7], dtype=numpy.uint16)
    seg2_bins = numpy.array([ 4,  3,  2,  2,  3,  1, 12,  4, 14, 13], dtype=numpy.uint16)
    seg3_bins = numpy.array([ 1,  0, 14,  9, 11,  9, 12, 13,  8,  2], dtype=numpy.uint16)
    seg4_bins = numpy.array([14,  6,  5, 10,  5, 13,  1, 11, 14, 13], dtype=numpy.uint16)

    expected_bins = numpy.vstack((seg1_bins, seg2_bins, seg3_bins, seg4_bins))

    def setUp(self):

        ##15 bins
        self.bins = [0.0] + [1.0+i for i in xrange(0,13)] + [14.0,float('inf')]

        self.bm = RectilinearBinMapper([self.bins])

        self.assign = self.bm.assign
    

    def test_assign_and_label1d_no_state(self):
        '''WAssign.assign_and_label : 1d binning assignments with no specified states successful'''

        assignments, trajlabels = WAssign.assign_and_label(self.nsegs, self.npts, self.parent_ids,
                                                           self.assign, None, None, self.pcoords)

        assignments = numpy.array(assignments)
        assert assignments.shape == (4, 10)
        
        numpy.testing.assert_array_equal(assignments, self.expected_bins)




