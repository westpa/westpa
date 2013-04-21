# Copyright (C) 2013 Nick Rego
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
import sys, os, time

from w_assign import WAssign, _assign_and_label
from w_kinetics import WKinetics
from westtools.h5io import WESTPAH5File
from westpa.binning._assign import assign_and_label
from westpa.binning.assign import RectilinearBinMapper, index_dtype, UNKNOWN_INDEX
import h5py
from work_managers import *

class DummyDataReader:
    '''Simple 'data_reader' facade to test full functionality of tools with sample data'''

    def __init__(self, max_iter=2, iter_groups=None):
        self.current_iteration = max_iter

        #Itergroups is a list of n 'iter groups'
        self.iter_groups = iter_groups

    def get_iter_group(self, iter_index = 1):

        return self.iter_groups[iter_index-1]

class DummyFile:
    '''mimics an hdf5 file that is the input for w_kinetics or w_kinavg'''

    def __init__(self, info_dict=dict(), attrs=dict()):
        self.info = info_dict
        self.attrs = attrs

    def __setitem__(self, key, value):
        self.info[key] = value

    def __getitem__(self, key):
        return self.info[key]

class WToolBase:
    '''
    Base class containing the test case'''

    
    ##Random, reverse ordered pcoords for the 4 segments, over 10 time points
    pcoords = numpy.array([[  8.6250176 ,   8.06916037,   7.1339341 ,   5.42291169, 5.18241533,
                              3.99005367,   3.45604335,   1.75481853,   1.67204737,   0.94418644],
                           [ 14.90786241,  13.57244566,  12.54785261,   4.89593824, 4.0411646 , 
                             3.65634435,   3.13918444,   2.95876746, 2.1180446 ,   1.17418433],
                           [ 14.47342808,  13.23377607,  12.49882247,  11.53806799, 9.85372988, 
                             9.48601102,   8.48030094,   2.63759144, 1.42639175,   0.84988666],
                           [ 14.57179029,  14.39726662,  13.99786996,  13.57312496, 11.01858827, 
                             10.49378797,   6.31455601,   5.99651334, 5.03735065,   1.62226425]])

    pcoords = pcoords.reshape(4,10,1)
    parent_ids = numpy.array([-1, -1, -1, -1])
    weights = numpy.array([0.25, 0.25, 0.25, 0.25])
    nsegs, npts = pcoords.shape[0:2]

    data_reader = DummyDataReader(iter_groups = [{'pcoord': pcoords, 'seg_index':{'parent_id':parent_ids,
                                                                                  'weight':weights}}])

    #2 states
    states = [{'coords': numpy.array([[ 0.],[ 1.],[ 2.],[ 3.]]),'label': 'bound'},
              {'coords': numpy.array([[ 4.],[ 5.],[ 6.],[ 7.],[ 8.],[ 9.],[ 10.],[ 11.],[ 12.],[ 13.],[ 14.],[ 15.]]),
               'label': 'unbound'}]

    state_labels = ['bound', 'unbound']
    state_map = numpy.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.uint16)

    last_labels = numpy.empty((4,), index_dtype)
    last_labels[:] = UNKNOWN_INDEX

    #Expected bin assignments for each segment
    expected_bins = numpy.array([[ 8,  8,  7,  5,  5,  3,  3,  1,  1,  0],
                                 [14, 13, 12,  4,  4,  3,  3,  2,  2,  1],
                                 [14, 13, 12, 11,  9,  9,  8,  2,  1,  0],
                                 [14, 14, 13, 13, 11, 10,  6,  5,  5,  1]], dtype=numpy.uint16)

    #Expected macrostate assignments for each segment
    expected_states = numpy.array([[1, 1,  1,  1,  1,  0,  0,  0,  0,  0],
                                   [1, 1,  1,  1,  1,  0,  0,  0,  0,  0],
                                   [1, 1,  1,  1,  1,  1,  1,  0,  0,  0],
                                   [1, 1,  1,  1,  1,  1,  1,  1,  1,  0]], dtype=numpy.uint16)

    #Flux from state 1 ('unbound') to state 0 ('bound') is 1 - all segments start unbound, finish bound
    expected_macro_fluxes = numpy.array([[0., 0.],
                                         [1., 0.]])

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
                                                           self.bm, self.state_map, self.last_labels, self.pcoords)

        assignments = numpy.array(assignments)
        trajlabels = numpy.array(trajlabels)
        assert assignments.shape == (4, 10)
        assert trajlabels.shape == (4, 10)
        
        numpy.testing.assert_array_equal(assignments, self.expected_bins)
        numpy.testing.assert_array_equal(trajlabels, self.expected_states)

    def test_assign_and_label_slice(self):
        '''WAssign.assign_and_label : Binning assignments successful over a slice of segments'''

        #First two segments
        assignments, trajlabels, lb, ub = _assign_and_label(0, self.nsegs//2, self.npts, self.parent_ids,
                                                           self.bm, self.state_map, self.last_labels, self.pcoords)       

        assignments = numpy.array(assignments)
        trajlabels = numpy.array(trajlabels)
        assert assignments.shape == (2, 10)
        assert trajlabels.shape == (2, 10)

        numpy.testing.assert_array_equal(assignments, self.expected_bins[0:2, :])
        numpy.testing.assert_array_equal(trajlabels, self.expected_states[0:2, :])

        #Second two segments
        assignments, trajlabels, lb, ub = _assign_and_label(self.nsegs//2, self.nsegs, self.npts, self.parent_ids,
                                                           self.bm, self.state_map, self.last_labels, self.pcoords)       

        assignments = numpy.array(assignments)
        trajlabels = numpy.array(trajlabels)
        assert assignments.shape == (2, 10)
        assert trajlabels.shape == (2, 10)

        numpy.testing.assert_array_equal(assignments, self.expected_bins[2:4, :])
        numpy.testing.assert_array_equal(trajlabels, self.expected_states[2:4, :])

class Test_W_Assign(WToolBase):
    '''Tests the full WAssign functionality'''

    def setUp(self):
        assert 'WM_WORK_MANAGER' not in os.environ
        assert 'WM_N_WORKERS' not in os.environ
        assert 'WM_ZMQ_INFO' not in os.environ
        self.w = WAssign()
        self.w.binning.mapper = self.bm
        self.w.states = self.states
        self.w.data_reader = self.data_reader
        self.w.output_file = WESTPAH5File('test_w_assign.h5', 'w')

    def tearDown(self):
        try:
            self.w.output_file.close()
            os.unlink('test_w_assign.h5')
            self.w.work_manager.shutdown()
        finally:
            del self.w.binning.mapper
            del self.w.work_manager
            del self.w
            os.environ.pop('WM_WORK_MANAGER', None)
            os.environ.pop('WM_N_WORKERS', None)
            os.environ.pop('WM_ZMQ_INFO', None)

    def test_go_serial(self):
        '''WAssign: works as expected using 'serial' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'serial'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1

        with self.w.work_manager:
            self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]
            trajlabels = data['trajlabels'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)
            assert trajlabels.shape == (self.nsegs, self.npts), 'Shape of trajlabels ({!r}) does not match expected ({},{})'.format(trajlabels.shape, self.nseg, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)
            numpy.testing.assert_array_equal(trajlabels, self.expected_states)

    def test_go_threads(self):
        '''WAssign: works as expected using 'threads' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'threads'
        os.environ['WM_N_WORKERS'] = '1'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1

        with self.w.work_manager:
            self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]
            trajlabels = data['trajlabels'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)
            assert trajlabels.shape == (self.nsegs, self.npts), 'Shape of trajlabels ({!r}) does not match expected ({},{})'.format(trajlabels.shape, self.nseg, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)
            numpy.testing.assert_array_equal(trajlabels, self.expected_states)

    def test_go_processes(self):
        '''WAssign: works as expected using 'processes' work manager'''

        os.environ['WM_WORK_MANAGER'] = 'processes'
        os.environ['WM_N_WORKERS'] = '1'
        self.w.work_manager = self.w.wm_env.make_work_manager()
        assert self.w.work_manager.n_workers == 1
        
        with self.w.work_manager:
            self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]
            trajlabels = data['trajlabels'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)
            assert trajlabels.shape == (self.nsegs, self.npts), 'Shape of trajlabels ({!r}) does not match expected ({},{})'.format(trajlabels.shape, self.nseg, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)
            numpy.testing.assert_array_equal(trajlabels, self.expected_states)

    def test_go_simple_zmq(self):
        '''WAssign: works as expected using a simple 'ZMQ' server with a 1-worker internal client'''

        os.environ['WM_WORK_MANAGER'] = 'zmq' 
        os.environ['WM_N_WORKERS'] = '1'
        self.w.work_manager = self.w.wm_env.make_work_manager() ##Defaults should give a non-dedicated server
        assert self.w.work_manager.n_local_workers == 1
        
        with self.w.work_manager:
            self.w.go()

        try: 
            data = h5py.File('test_w_assign.h5')
        except IOError:
            raise IOError('Error opening hdf5 file')
        else:
            assignments = data['assignments'][0]
            trajlabels = data['trajlabels'][0]

            assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)
            assert trajlabels.shape == (self.nsegs, self.npts), 'Shape of trajlabels ({!r}) does not match expected ({},{})'.format(trajlabels.shape, self.nseg, self.npts)

            numpy.testing.assert_array_equal(assignments, self.expected_bins)
            numpy.testing.assert_array_equal(trajlabels, self.expected_states)

    #@SkipTest
    def test_go_complex_zmq(self):
        '''WAssign: works as expected using a ZMQWorkManager with an external client'''

        os.environ['WM_WORK_MANAGER'] = 'zmq'
        os.environ['WM_N_WORKERS'] = '0' #Dedicated server
        os.environ['WM_ZMQ_INFO'] = 'server_info.json' #Server endpoint info

        server = self.w.work_manager = self.w.wm_env.make_work_manager() 
        assert server.n_local_workers == 0
        
        with server:
            os.environ['WM_N_WORKERS'] = '4' #for a client with 4 workers

            client = ZMQClient.from_environ()
            client.startup()

            time.sleep(0.11) #Some time for server and client to sync up

            self.w.go()

            try:
                assert self.w.work_manager.n_workers == 4
                try:
                    data = h5py.File('test_w_assign.h5')
                except IOError:
                    raise IOError('Error opening hdf5 file')
                else:
                    assignments = data['assignments'][0]
                    trajlabels = data['trajlabels'][0]

                    assert assignments.shape == (self.nsegs, self.npts), 'Shape of assignments ({!r}) does not match expected ({},{})'.format(assignments.shape, self.nsegs, self.npts)
                    assert trajlabels.shape == (self.nsegs, self.npts), 'Shape of trajlabels ({!r}) does not match expected ({},{})'.format(trajlabels.shape, self.nseg, self.npts)

                    numpy.testing.assert_array_equal(assignments, self.expected_bins)
                    numpy.testing.assert_array_equal(trajlabels, self.expected_states)

            finally:
                client.shutdown()


class Test_W_Kinetics(WToolBase):

    assignments_file = DummyFile()
    assignments_file.attrs['nbins'] = WToolBase.bm.nbins
    assignments_file.attrs['iter_start'] = 1
    assignments_file.attrs['iter_stop'] = 2
    assignments_file['state_labels'] = numpy.array(WToolBase.state_labels)
    assignments_file['state_map'] = WToolBase.state_map
    assignments_file['assignments'] = WToolBase.expected_bins.reshape(1,4,10)
    assignments_file['trajlabels'] = WToolBase.expected_states.reshape(1,4,10)

    def setUp(self):
        self.w = WKinetics()
        self.w.assignments_file = self.assignments_file
        self.w.data_reader = self.data_reader
        self.w.iter_range.iter_start, self.w.iter_range.iter_stop = 1,2
        self.w.window_size = 1
        self.w.output_file = WESTPAH5File('test_w_kinetics.h5', 'w')

    def tearDown(self):
        self.w.output_file.close()
        os.unlink('test_w_kinetics.h5')
        del self.w

    def test_go(self):

        self.w.go()

        try:
            data = h5py.File('test_w_kinetics.h5')

        except IOError:
            raise IOError('Error opening hdf5 file')

        else:
            labeled_bin_fluxes = data['labeled_bin_fluxes'][0]

            ##Bin fluxes is 4 dimensional - not going to reconstruct entire expected array
            # we know that all segments start 'unbound' (state 1) and end up 'bound' (state 0)
            # and all segments start at either bin 9 or bin 15, and end up in the first or second bins,
            # so the only bin fluxes we worry about are labeled_bin_fluxes[1][0][8 | 14][0 | 1]
            assert labeled_bin_fluxes[1][0][8][0] == 0.25
            assert labeled_bin_fluxes[1][0][14][0] == 0.25
            assert labeled_bin_fluxes[1][0][14][1] == 0.50

            #Macrostate fluxes is easier to compare directly to expected array - only 2-d (2 states by 2 states)
            macrostate_fluxes = data['trace_macro_fluxes'][0]

            numpy.testing.assert_array_equal(macrostate_fluxes, self.expected_macro_fluxes)
