import argparse
import os
import sys

import h5py
import pytest
from mdtraj import Trajectory
from numpy.testing import assert_allclose

from westpa.core._rc import WESTRC
from westpa.core.data_manager import WESTDataManager
from westpa.core.segment import Segment
from westpa.core.trajectory import WESTTrajectory, load_mdanalysis, load_mdtraj, load_netcdf
from westpa.core.propagators.loaders import mdtraj_trajectory_loader, mdanalysis_trajectory_loader, netcdf_trajectory_loader


class TestHDF5Framework:
    """Class to test HDF5 Framework"""

    @pytest.mark.filterwarnings("ignore:Element information is missing")
    def test_load_mdanalysis(self, traj_setup):
        pytest.importorskip('MDAnalysis')

        test_traj = load_mdanalysis(self.current_path)

        assert isinstance(test_traj, WESTTrajectory)

        # load_mdanalysis automatically converts to nm
        assert_allclose(test_traj.xyz / 10, self.ref_coords)
        assert_allclose(test_traj.time, self.ref_time)

    def test_load_netcdf(self, traj_setup):
        test_traj = load_netcdf(self.current_path)

        assert isinstance(test_traj, WESTTrajectory)

        assert_allclose(test_traj.xyz, self.ref_coords)
        assert_allclose(test_traj.time, self.ref_time)

    def test_load_mdtraj(self, traj_setup):
        test_traj = load_mdtraj(self.current_path)

        assert isinstance(test_traj, Trajectory)

        assert_allclose(test_traj.xyz, self.ref_coords)
        assert_allclose(test_traj.time, self.ref_time)

    @pytest.mark.filterwarnings("ignore:Element information is missing")
    def test_mdanalysis_trajectory_loader(self, traj_setup, monkeypatch):
        pytest.importorskip('MDAnalysis')

        dummy_segment = Segment()

        mdanalysis_trajectory_loader('dummy', self.current_path, dummy_segment, False)

        assert_allclose(dummy_segment.data['iterh5/trajectory'].xyz / 10, self.ref_coords)
        assert_allclose(dummy_segment.data['iterh5/trajectory'].time, self.ref_time)

    def test_mdanalysis_trajectory_loader_fail_import(self, traj_setup, monkeypatch):
        """Test fallback to MDTraj with `mdanalysis_trajectory_loader`"""

        dummy_segment = Segment()

        with monkeypatch.context() as m:
            m.setitem(sys.modules, 'MDAnalysis', None)
            mdanalysis_trajectory_loader('dummy', self.current_path, dummy_segment, False)

        assert_allclose(dummy_segment.data['iterh5/trajectory'].xyz, self.ref_coords)
        assert_allclose(dummy_segment.data['iterh5/trajectory'].time, self.ref_time)

    def test_netcdf_trajectory_loader(self, traj_setup):
        dummy_segment = Segment()
        netcdf_trajectory_loader('dummy', self.current_path, dummy_segment, False)

        assert_allclose(dummy_segment.data['iterh5/trajectory'].xyz, self.ref_coords)
        assert_allclose(dummy_segment.data['iterh5/trajectory'].time, self.ref_time)

    def test_mdtraj_trajectory_loader(self, traj_setup):
        dummy_segment = Segment()
        mdtraj_trajectory_loader('dummy', self.current_path, dummy_segment, False)

        assert_allclose(dummy_segment.data['iterh5/trajectory'].xyz, self.ref_coords)
        assert_allclose(dummy_segment.data['iterh5/trajectory'].time, self.ref_time)

    def test_update_iterh5_file(self, ref_initialized, monkeypatch, tmp_path):
        """
        Test for checking if update_iter_h5file actually finds the correct relative paths
        relative to `west.h5` and not your current working directory.
        The external link would previously be something like `../abc/iter_000000.h5`
        instead of `iter_000000.h5`, even if your `west.h5` file is right next to it.

        Conditions necessary for this to happen:
        a) $WEST_SIM_ROOT is an absolute path.
        b) $WEST_SIM_ROOT must exist.
        c) You are working in a location other than $WEST_SIM_ROOT.
        d) `WESTDataManager.we_h5filename` (usually read from `west.cfg`
           under `west/data/west_data_file`) either doesn't exist (so the `west.h5`
           default is used) or is a relative path.

        EXAMPLE TREE
        |
        |---abc  => $WEST_SIM_ROOT
        |   |--- west.h5
        |   |--- iter_000000.h5
        |
        |---def  => $PWD

        """
        # Setup folder structure and move to abc/
        os.makedirs('abc', exist_ok=True)
        os.makedirs('def', exist_ok=True)
        os.chdir('abc')

        # Setup rc
        rc = WESTRC()
        args = argparse.Namespace(
            verbosity='debug',
            rcfile='../' + self.cfg_filepath,
            we_h5filename=self.h5_filepath,
        )
        rc.process_args(args)

        # Setup data manager and create west.h5 file
        dm = WESTDataManager(rc=rc)
        dm.iter_h5_path_template = '$WEST_SIM_ROOT/iter_{n_iter:06d}.h5'
        dm.we_h5filename = 'west.h5'
        dm.store_h5 = True
        dm.prepare_backing()

        with monkeypatch.context() as m:
            # map $WEST_SIM_ROOT to abc, move to def/
            m.setenv('WEST_SIM_ROOT', f'{tmp_path}/abc')
            os.chdir('../def')

            # Write to west.h5, create iter_000000.h5
            dm.require_iter_group(0)
            dm.update_iter_h5file(0, [])

            # Check files existence (should be in $WEST_SIM_ROOT!) and link location
            assert os.path.exists(f'{tmp_path}/abc/iter_000000.h5')
            assert not os.path.exists(f'{tmp_path}/def/iter_000000.h5')
            with h5py.File('../abc/west.h5', 'r') as h5_iterfile:
                link_path = h5_iterfile.get('iterations/iter_00000000/trajectories', getlink=True).filename
                assert link_path == 'iter_000000.h5'  # Before v2022.16 it would be '../abc/iter_000000.h5'
