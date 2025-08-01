from unittest import mock
import argparse

from westpa.cli.tools.w_timings import entry_point


class Test_W_Timings:
    """Test class for w_timings tool."""

    def args(self, first_iter=None, last_iter=None, tau=None):
        return argparse.Namespace(
            verbosity=None,
            rcfile=None,
            we_h5filename=self.h5_filepath,
            first_iter=first_iter,
            last_iter=last_iter,
            tau=tau,
        )

    def test_default(self, ref_50iter, capsys):
        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=self.args(),
        ):
            entry_point()
        captured = capsys.readouterr()
        output = captured.out
        expected = """\
Iterations:                50
Total segments:            9985
Wall-clock time:           0:00:01.452625
"""
        assert output == expected

    def test_tau(self, ref_50iter, capsys):
        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=self.args(tau='100_ps'),
        ):
            entry_point()
        captured = capsys.readouterr()
        output = captured.out
        expected = """\
Iterations:                50
Total segments:            9985
Wall-clock time:           0:00:01.452625
Maximum trajectory length: 5.0 ns
Aggregate simulation time: 998.5 ns
"""
        assert output == expected

    def test_tau_invalid_unit(self, ref_50iter, capsys):
        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=self.args(tau='100_u'),
        ):
            try:
                entry_point()
            except SystemExit as e:
                assert e.code == 2
        captured = capsys.readouterr()
        assert "'u' is not a recognized time unit" in captured.err

    def test_tau_invalid_format(self, ref_50iter, capsys):
        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=self.args(tau='100'),
        ):
            try:
                entry_point()
            except SystemExit as e:
                assert e.code == 2
        captured = capsys.readouterr()
        assert 'must be formatted as <value>_<unit>' in captured.err

    def test_iter_range(self, ref_50iter, capsys):
        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=self.args(first_iter=11, last_iter=40),
        ):
            entry_point()
        captured = capsys.readouterr()
        output = captured.out
        expected = """\
Iterations:                30
Total segments:            6320
Wall-clock time:           0:00:00.868997
"""
        assert output == expected
