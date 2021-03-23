import h5py

import westpa.cli.core.w_run
import westpa.work_managers as work_managers


class Test_W_Run:
    def test_run_w_run(self, ref_initialized):
        '''Tests running an initialized WESTPA system for 3 iterations'''

        class TestArgs:
            def __init__(
                self,
                n_workers,
                only_one_segment,
                rcfile,
                verbosity,
                work_manager,
                zmq_comm_mode,
                zmq_downstream_ann_endpoint,
                zmq_downstream_rr_endpoint,
                zmq_master_heartbeat,
                zmq_mode,
                zmq_read_host_info,
                zmq_shutdown_timeout,
                zmq_startup_timeout,
                zmq_timeout_factor,
                zmq_upstream_ann_endpoint,
                zmq_upstream_rr_endpoint,
                zmq_worker_heartbeat,
                zmq_write_host_info,
            ):

                self.n_workers = n_workers
                self.only_one_segment = only_one_segment
                self.rcfile = rcfile
                self.verbosity = verbosity
                self.work_manager = work_manager
                self.zmq_comm_mode = zmq_comm_mode
                self.zmq_downstream_ann_endpoint = zmq_downstream_ann_endpoint
                self.zmq_downstream_rr_endpoint = zmq_downstream_rr_endpoint
                self.zmq_master_heartbeat = zmq_master_heartbeat
                self.zmq_mode = zmq_mode
                self.zmq_read_host_info = zmq_read_host_info
                self.zmq_shutdown_timeout = zmq_shutdown_timeout
                self.zmq_startup_timeout = zmq_startup_timeout
                self.zmq_timeout_factor = zmq_timeout_factor
                self.zmq_upstream_ann_endpoint = zmq_upstream_ann_endpoint
                self.zmq_upstream_rr_endpoint = zmq_upstream_rr_endpoint
                self.zmq_worker_heartbeat = zmq_worker_heartbeat
                self.zmq_write_host_info = zmq_write_host_info

        args = TestArgs(
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
        )

        westpa.rc.process_args(args)
        work_managers.environment.process_wm_args(args)

        westpa.cli.core.w_run.run_simulation()

        # Since the dynamics produce slightly different output in each iteration, just check that
        #    the third iteration exists in the h5 file.
        # TODO: Modify the ODLD system to run identically across runs, and then remove the excluded
        #   datasets.

        _hfile = h5py.File(self.h5_filepath)
        assert 'iter_00000003' in _hfile['/iterations'].keys()
        _hfile.close()
