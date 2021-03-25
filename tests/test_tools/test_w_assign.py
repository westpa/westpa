import shutil

from h5diff import H5Diff

from westpa.cli.tools import w_assign
from common import MockArgs


class Test_W_Assign:
    def test_run_w_assign(self, ref_50iter):

        args = MockArgs(
            verbosity='debug',
            rcfile=self.cfg_filepath,
            max_queue_length=None,
            we_h5filename=self.h5_filepath,
            construct_dataset=None,
            dsspecs=None,
            output='assign.h5',
            subsample=None,
            config_from_file=True,
            scheme='TEST',
        )

        # This basically some logic that's wrapped up in WESTTool.main() for convenience.
        # It needs to be explicitly called like this because the args are captured and set in make_parser_and_process()
        #   which we don't want to call, because we don't want to make a parser.
        #   We just want to process the args that "would've" been captured if called from CLI.
        tool = w_assign.WAssign()

        # Prepare and instantiate work manager
        tool.wm_env.process_wm_args(args)
        tool.work_manager = tool.wm_env.make_work_manager()

        tool.process_all_args(args)
        with tool.work_manager:
            if tool.work_manager.is_master:
                tool.go()
            else:
                tool.work_manager.run()

        diff = H5Diff('./assign_ref.h5', './ANALYSIS/TEST/assign.h5')
        diff.check()

        # clean up
        shutil.rmtree('ANALYSIS')
