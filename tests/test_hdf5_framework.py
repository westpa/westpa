from numpy.testing import assert_allclose
from mdtraj import Trajectory

from westpa.core.trajectory import WESTTrajectory, load_mda, load_mdtraj, load_netcdf


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
