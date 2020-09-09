import argparse

import h5py

from westpa.cli.core.w_run import entry_point
from unittest import mock


class Test_W_Run:
    def test_run_w_run(self, ref_initialized):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        with mock.patch(
            target='argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(
                n_workers=None,
                only_one_segment=False,
                rcfile=self.cfg_filepath,
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

        _hfile = h5py.File(self.h5_filepath)
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()
