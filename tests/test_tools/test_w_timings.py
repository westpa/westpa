import os
from unittest import mock
import argparse
import pytest
from westpa.cli.tools.w_timings import entry_point

class Test_W_Timings:
    """Test class for w_timings tool."""

    @pytest.fixture
    def ref_west_file(self):
        return "tests/refs/west_3iter.h5"

    def test_run_w_timings(self, ref_west_file, capsys):
        """
        Currently only checks if sections headers are printing correctly
        """

        with mock.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=argparse.Namespace(
                verbosity="debug",
                rcfile=None,
                we_h5filename=ref_west_file,
                first_iter=None,
                last_iter=None,
                tau=100,
                count_events=True,
            ),
        ):
            entry_point()

        captured = capsys.readouterr()
        output = captured.out

        # Basic checks that output sections are present
        #TODO: Replace with actual known values from reference
        assert "===== WALLCLOCK  =====" in output
        assert "Total Wallclock Time:" in output
        assert "===== SIMULATION  =====" in output
        assert "Simulation time:" in output
        assert "===== RECYCLING =====" in output
        assert "Recycled walkers:" in output

