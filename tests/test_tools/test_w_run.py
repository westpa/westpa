import unittest
import argparse
import os

from .hdiff import H5Diff

from westpa.cli.core.w_run import entry_point
from unittest import mock

from shutil import copy


class Test_W_Run(unittest.TestCase):
    test_name = 'W_RUN'

    def __init__(self, methodName):

        super().__init__(methodName='test_run_w_run')

    def setUp(self):

        self.starting_path = os.getcwd()
        self.initial_PWD = os.environ['PWD']

        self.odld_path = os.path.dirname(__file__) + '/ref'

        copy(self.odld_path + '/west_ref.h5', self.starting_path + '/west.h5')
        copy(self.odld_path + '/west.cfg', self.starting_path + '/west.cfg')

    def test_run_w_run(self):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                n_workers=None,
                only_one_segment=False,
                rcfile=self.starting_path + '/west.cfg',
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

        # Verify that the generated h5 file is similar to the 3-iteration reference
        #   I say similar and not identical because I haven't found the right set of flags
        #   for running identically reproducible WESTPA runs.
        # TODO: Modify the ODLD system to run identically across runs, and then remove the excluded
        #   datasets.
        diff = H5Diff(self.odld_path + '/west_3iter.h5', self.starting_path + '/west.h5')
        diff.check(
            float_thresh=0.1,
            excluded_datasets=[
                'iterations/iter_00000002/wtgraph',
                'iterations/iter_00000001/seg_index',
                'iterations/iter_00000003/seg_index',
                'iterations/iter_00000003/wtgraph',
                'summary',
            ],
        )

    def tearDown(self):

        pass
        os.remove(self.starting_path + '/west.h5')
        os.remove(self.starting_path + '/west.cfg')
