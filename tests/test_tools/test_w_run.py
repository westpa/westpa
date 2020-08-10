import unittest
import argparse
import os

import h5py

import westpa

from westpa.cli.core.w_run import entry_point
from unittest import mock

from shutil import copy


class Test_W_Run(unittest.TestCase):
    test_name = 'W_RUN'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_run')

    def setUp(self):

        self.starting_path = os.getcwd()
        # self.initial_PWD = os.environ['PWD']

        self.odld_path = os.path.dirname(__file__) + '/ref'
        os.chdir(self.odld_path)
        os.environ['WEST_SIM_ROOT'] = self.odld_path
        westpa.rc.config = westpa.core.yamlcfg.YAMLConfig()

        copy(self.odld_path + '/west_ref.h5', self.odld_path + '/west.h5')

    def test_run_w_run(self):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                n_workers=None,
                only_one_segment=False,
                rcfile=self.odld_path + '/west.cfg',
                verbosity='0',
                work_manager=None,
                zmq_comm_mode=None,
                zmq_downstream_ann_endpoint=None,
                zmq_downstream_rr_endpoint=None,
                zmq_master_heartbeat=None,
                zmq_mode=None,
                zmq_read_host_info=None,
                zmq_shutdown_timeout=None,
                zmq_startup_timeout=None,
                zmq_timeout_factor=None,
                zmq_upstream_ann_endpoint=None,
                zmq_upstream_rr_endpoint=None,
                zmq_worker_heartbeat=None,
                zmq_write_host_info=None,
            ),
        ):
            entry_point()

        # Since the dynamics produce slightly different output in each iteration, just check that
        #    the third iteration exists in the h5 file.
        # TODO: Modify the ODLD system to run identically across runs, and then remove the excluded
        #   datasets.

        _hfile = h5py.File(self.odld_path + '/west.h5')
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()

    def tearDown(self):

        westpa.rc._sim_manager = None
        westpa.rc._system = None
        westpa.rc._data_manager = None
        westpa.rc._we_driver = None
        westpa.rc._propagator = None
        os.environ['WEST_SIM_ROOT'] = ''

        os.remove(self.odld_path + '/west.h5')
        os.chdir(self.starting_path)
