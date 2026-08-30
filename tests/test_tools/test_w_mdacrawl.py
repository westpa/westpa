import argparse
import pytest
import os
from unittest import mock

from westpa.cli.tools.w_mdacrawl import entry_point

hdf5_file_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
hdf5_file = os.path.join(hdf5_file_parent, 'refs', 'west_mdacrawl.h5')


class Test_W_MDACrawl:
    """Test class for w_mdacrawl tool"""

    def test_cli_execution(self, capsys):
        mock_args = argparse.Namespace(
            verbosity=None,
            rcfile=None,
            work_manager=None,
            west_h5=hdf5_file,
            analysis="MDAnalysis.analysis.rms.RMSD",
            column=None,
            save=None,
            select="all",  # Prevents NaN math errors in CI
            n_workers=1,
            overwrite=False,
        )

        with mock.patch("argparse.ArgumentParser.parse_args", return_value=mock_args):
            entry_point()

        output = capsys.readouterr().out

        assert "Analysis : RMSD" in output
        assert "Results were NOT saved" in output

    def test_invalid_module_import(self):
        mock_args = argparse.Namespace(
            verbosity=None,
            rcfile=None,
            work_manager=None,
            west_h5=hdf5_file,
            analysis="MDAnalysis.analysis.rms.FakeTool",
            column=None,
            save=None,
            select="all",
            n_workers=1,
            overwrite=False,
        )

        with mock.patch("argparse.ArgumentParser.parse_args", return_value=mock_args):
            with pytest.raises(AttributeError):
                entry_point()
