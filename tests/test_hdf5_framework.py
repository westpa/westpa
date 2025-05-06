from numpy.testing import assert_allclose
from mdtraj import Trajectory

from westpa.core.segment import Segment
from westpa.core.trajectory import WESTTrajectory, load_mda, load_mdtraj, load_netcdf
from westpa.core.propagators.executable import mdtraj_trajectory_loader, mda_trajectory_loader, netcdf_trajectory_loader


class TestHDF5Framework:
    '''Class to test HDF5 Framework'''

    def test_load_mda(self, traj_setup):
        test_traj = load_mda(self.current_path)

        assert isinstance(test_traj, WESTTrajectory)

        # load_mda automatically converts to nm
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

    def test_mda_trajectory_loader(self, traj_setup):
        dummy_segment = Segment()
        mda_trajectory_loader('dummy', self.current_path, dummy_segment, False)

        assert_allclose(dummy_segment.data['iterh5/trajectory'].xyz / 10, self.ref_coords)
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
