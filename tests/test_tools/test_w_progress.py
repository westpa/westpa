from datetime import datetime
from unittest import mock
import argparse
import subprocess
import sys

import pytest

from westpa.cli.tools.w_progress import entry_point, format_progress, read_progress


class Test_W_Progress:
    """Test class for w_progress tool."""

    def args(self, we_h5filename=None, refresh=1.0):
        return argparse.Namespace(
            verbosity=None,
            rcfile=None,
            we_h5filename=we_h5filename or self.h5_filepath,
            refresh=refresh,
        )

    def test_completed_run(self, ref_50iter):
        progress = read_progress(self.h5_filepath, max_total_iterations=50)

        assert progress['n_iter'] == 51
        assert progress['n_completed'] == 50
        assert progress['n_particles'] == 9985
        assert len(progress['recent_walltimes']) == 5
        assert progress['avg_walltime'] > 0
        assert progress['eta'] == 0

    def test_eta(self, ref_50iter):
        progress = read_progress(self.h5_filepath, max_total_iterations=60)

        assert progress['eta'] == pytest.approx(10 * progress['avg_walltime'])

    def test_initialized_run(self, ref_initialized):
        progress = read_progress(self.h5_filepath, max_total_iterations=2)

        assert progress['n_iter'] == 1
        assert progress['n_completed'] == 0
        assert progress['n_segs'] == 10
        assert progress['n_complete'] == 0
        assert progress['n_particles'] == 0
        assert progress['recent_walltimes'] == []
        assert progress['eta'] is None

    def test_file_open_for_writing(self, ref_50iter):
        # Lock west.h5 from another process like w_run does
        code = f"import h5py, time; f = h5py.File({self.h5_filepath!r}, 'r+'); print('ready', flush=True); time.sleep(60)"
        writer = subprocess.Popen([sys.executable, '-c', code], stdout=subprocess.PIPE, text=True)
        try:
            assert writer.stdout.readline().strip() == 'ready'
            progress = read_progress(self.h5_filepath)
        finally:
            writer.kill()
            writer.wait()

        assert progress['n_iter'] == 51

    def test_format_progress(self):
        progress = {
            'mtime': datetime(2026, 1, 2, 3, 4, 5).timestamp(),
            'n_iter': 37,
            'n_completed': 36,
            'max_total_iterations': 50,
            'n_segs': 100,
            'n_complete': 82,
            'n_failed': 2,
            'n_particles': 9985,
            'walltime': 642.0,
            'recent_walltimes': [12.4, 11.8, 13.1],
            'avg_walltime': 12.43,
            'eta': 174.07,
        }
        expected = """\
Last written:              2026-01-02 03:04:05
Current iteration:         37
Progress:                  36 / 50 iterations (72.0%)
Segments complete:         82 / 100
Segments failed:           2
Recent walltimes:          0:00:12, 0:00:12, 0:00:13
Avg iteration time:        0:00:12
Completed walltime:        0:10:42
Completed segments:        9985
ETA:                       0:02:54
"""
        assert format_progress(progress) == expected

    def test_format_progress_unknown(self):
        progress = {
            'mtime': 0,
            'n_iter': 1,
            'n_completed': 0,
            'max_total_iterations': None,
            'n_segs': 10,
            'n_complete': 0,
            'n_failed': 0,
            'n_particles': 0,
            'walltime': 0.0,
            'recent_walltimes': [],
            'avg_walltime': None,
            'eta': None,
        }
        output = format_progress(progress)

        assert 'Progress:                  0 completed\n' in output
        assert 'Recent walltimes:          unknown\n' in output
        assert 'ETA:                       unknown\n' in output

    def test_default(self, ref_50iter, capsys):
        with mock.patch('argparse.ArgumentParser.parse_args', return_value=self.args()):
            with mock.patch('time.sleep', side_effect=KeyboardInterrupt):
                entry_point()
        output = capsys.readouterr().out

        assert output.startswith(f'WESTPA progress for {self.h5_filepath} (updated ')
        assert 'Current iteration:         51\n' in output
        assert 'Progress:                  50 / 50 iterations (100.0%)\n' in output

    def test_missing_file(self, ref_50iter, capsys):
        with mock.patch('argparse.ArgumentParser.parse_args', return_value=self.args(we_h5filename='missing.h5')):
            with mock.patch('time.sleep', side_effect=KeyboardInterrupt):
                entry_point()
        output = capsys.readouterr().out

        assert 'Could not read missing.h5' in output

    def test_refresh_must_be_positive(self, ref_50iter):
        with mock.patch('argparse.ArgumentParser.parse_args', return_value=self.args(refresh=0)):
            with pytest.raises(SystemExit):
                entry_point()
