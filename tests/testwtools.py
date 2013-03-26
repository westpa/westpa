'''

Created 3/19/2013

author: nbr

Test suite for west_tools

'''

from __future__ import division, print_function; __metaclass__ = type
import nose.tools
from nose.tools import raises
from nose.plugins.skip import SkipTest

import numpy
import sys, os

from w_assign import WAssign, _assign_and_label
from westtools.h5io import WESTPAH5File
from westpa.binning._assign import assign_and_label
from westpa.binning.assign import RectilinearBinMapper
import h5py

class DummyDataReader:
    '''Simple 'data_reader' facade to test full functionality of tools with sample data'''

    def __init__(self, max_iter=2, iter_groups=None):
        self.current_iteration = max_iter

        #Itergroups is a list of n 'iter groups'
        self.iter_groups = iter_groups

    def get_iter_group(self, iter_index = 0):

        return self.iter_groups[0]

class WToolBase:
    '''
    Base class containing the test case'''
    

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

    ##15 bins
    bins = [0.0] + [1.0+i for i in xrange(0,13)] + [14.0,float('inf')]

    bm = RectilinearBinMapper([bins])

    assign = bm.assign

class Test_W_Assign_assign_and_label(WToolBase):
    '''
    Tests w_assign over a single iteration, using random pcoords for 4 segments over 10 timepoints'''
    

    def test_assign_and_label1d_all_segs(self):
        '''WAssign.assign_and_label : 1d binning assignments (all segs at once) successful'''

        assignments, trajlabels, lb, ub = _assign_and_label(0, self.nsegs, self.npts, self.parent_ids,
                                                           self.bm, None, None, self.pcoords)

        assignments = numpy.array(assignments)
        assert assignments.shape == (4, 10)
        
        numpy.testing.assert_array_equal(assignments, self.expected_bins)

    def test_assign_and_label_slice(self):
        '''WAssign.assign_and_label : Binning assignments successful over a slice of segments'''

        #First two segments
        assignments, trajlabels, lb, ub = _assign_and_label(0, self.nsegs//2, self.npts, self.parent_ids,
                                                           self.bm, None, None, self.pcoords)       

        assignments = numpy.array(assignments)
        assert assignments.shape == (2, 10)

        numpy.testing.assert_array_equal(assignments, self.expected_bins[0:2, :])

        #Second two segments
        assignments, trajlabels, lb, ub = _assign_and_label(self.nsegs//2, self.nsegs, self.npts, self.parent_ids,
                                                           self.bm, None, None, self.pcoords)       

        assignments = numpy.array(assignments)
        assert assignments.shape == (2, 10)

        numpy.testing.assert_array_equal(assignments, self.expected_bins[2:4, :])


class Test_W_Assign(WToolBase):
    '''Tests the full WAssign functionality'''

    def setUp(self):
        assert 'WM_WORK_MANAGER' not in os.environ
        self.w = WAssign()
        self.w.binning.mapper = self.bm
        self.w.data_reader = DummyDataReader(iter_groups = [{'pcoord': self.pcoords, 'seg_index':{'parent_id':self.parent_ids}}])
        self.w.output_file = WESTPAH5File('test_w_assign.h5', 'w')

    def tearDown(self):
        try:
            self.w.output_file.close()
            os.unlink('test_w_assign.h5')
            self.w.work_manager.shutdown()
        finally:
            del self.w.binning.mapper
            del self.w.work_manager
            del self.w.data_reader
            del self.w
            os.environ.pop('WM_WORK_MANAGER', None)

    def test_go_serial(self):
        '''WAssign: works as expected using 'serial' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'serial'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1

        self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)

    def test_go_threads(self):
        '''WAssign: works as expected using 'threads' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'threads'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1
        self.w.work_manager.startup()

        self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)

    #@SkipTest
    def test_go_processes(self):
        '''WAssign: works as expected using 'processes' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'processes'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1
        self.w.work_manager.startup()

        self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)

    @SkipTest
    def test_go_simple_zmq(self):
        '''WAssign: works as expected using a simple 'ZMQ' server with a 1-worker internal client'''

        os.environ['WM_WORK_MANAGER'] = 'zmq' 
        os.environ['WM_N_WORKERS'] = '1'
        self.w.work_manager = self.w.wm_env.make_work_manager() ##Defaults should give a non-dedicated server
        assert self.w.work_manager.n_local_workers == 1
        self.w.work_manager.startup()

        self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)





