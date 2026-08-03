import pytest
import os

import numpy as np
import h5py

mda = pytest.importorskip("MDAnalysis")
from westpa.core import mdacrawl  # noqa

here = os.path.dirname(os.path.abspath(__file__))
hdf5_file = os.path.join(here, 'refs', 'west_mdacrawl.h5')


class Test_WESTPAParser:
    """Class to test the WESTPAParser topology reading module"""

    @pytest.fixture(scope="class")
    def mda_universe(self):
        u = mda.Universe(hdf5_file, format='WESTPA')
        yield u

        if hasattr(u, 'trajectory'):
            u.trajectory.close()

    def test_parser_atom_count(self, mda_universe):
        """Parser extracts the correct number of atoms from Tutorial 7.5"""

        expected_atoms = 4010  # Known value

        assert len(mda_universe.atoms) == expected_atoms, f"Expected {expected_atoms} atoms, got {len(mda_universe.atoms)}"
        assert mda_universe.atoms.ids[-1] == expected_atoms - 1, "Atom IDs were not indexed properly"

    def test_parser_residue_names(self, mda_universe):
        """Parser extracts the correct residue names and counts from Tutorial 7.5"""

        # Known values
        expected_residues = 1338
        expected_first_five = ['Na+', 'Cl-', 'HOH', 'HOH', 'HOH']

        assert (
            len(mda_universe.residues) == expected_residues
        ), f"Expected {expected_residues} residues, got {len(mda_universe.residues)}"

        actual_first_five = mda_universe.residues.resnames[:5]
        assert np.array_equal(actual_first_five, expected_first_five), "First five residues did not match"

        # Testing internal search
        water_selection = mda_universe.select_atoms("resname HOH")
        assert len(water_selection) > 0, "Parser failed to map HOH residues correctly."

    def test_parser_masses(self, mda_universe):
        """Parser has correct mass values from mass dictionary"""

        masses = mda_universe.atoms.masses

        assert np.isclose(masses[0], 22.98977), "Sodium mass is mapped incorrectly"
        assert np.isclose(masses[1], 35.45), "Chlorine mass is mapped incorrectly"
        assert np.all(masses > 0.0), "One or more atoms failed mass mapping and defaulted to 0.0"

    def test_parser_bonds(self, mda_universe):
        """Parser converts JSON bond arrays to MDAnalysis bonds correctly"""

        expected_bonds = 2672  # Known value
        assert len(mda_universe.bonds) == expected_bonds, "Number of bonds did not match"

    def test_parser_segments(self, mda_universe):
        """Parser correctly populates segments"""

        assert len(mda_universe.segments) >= 1, "No segments were found in the topology"
        assert len(mda_universe.segments[0].atoms) > 0, "The generated segment is empty"

    def test_parser_elements(self, mda_universe):
        """Parser correctly maps atomic elements"""

        elements = mda_universe.atoms.elements

        assert len(elements) == len(mda_universe.atoms), "Elements array length mismatch"

        water_oxygens = mda_universe.select_atoms("resname HOH and name O")
        if len(water_oxygens) > 0:
            assert np.all(water_oxygens.elements == 'O'), "Water oxygen element is not 'O'"

    def test_parser_bond_connectivity(self, mda_universe):
        """Parser maps bonds to the correct atom indices"""

        water = mda_universe.select_atoms("resname HOH").residues[0]
        assert len(water.atoms.bonds) == 2, "Water molecule does not have exactly 2 bonds"

        # Check that the oxygen is bonded to hydrogens
        oxygen = water.atoms.select_atoms("name O")
        if len(oxygen) == 1:
            bonded_atoms = oxygen[0].bonded_atoms
            assert np.all(np.isin(bonded_atoms.names, ['H1', 'H2', 'HW1', 'HW2'])), "Oxygen bonded to non-hydrogen atom"


class Test_WESTPAReader:
    """Class to test the WESTPAReader dynamic trajectory reading module"""

    @pytest.fixture(scope="class")
    def mda_universe(self):
        u = mda.Universe(hdf5_file, format='WESTPA')
        yield u

        if hasattr(u, 'trajectory'):
            u.trajectory.close()

    def test_reader_frame_count(self, mda_universe):
        """Reader reads correct frame count"""
        frames = mda_universe.trajectory
        assert len(frames) > 0, "Trajectory should not be empty"
        assert len(frames) == 100, "Frame count does not match"

    def test_reader_coordinates(self, mda_universe):
        """Ensure coordinates are correctly extracted and scaled to Ångströms"""
        ts = mda_universe.trajectory[0]
        iter_num, seg_idx, actual_pos, path, local_frame = mda_universe.trajectory.frame_index[0]

        with h5py.File(path, 'r') as f:
            coords = f['coordinates'][actual_pos] * 10

        assert ts.positions is not None
        assert ts.positions.dtype == np.float32
        assert ts.positions.shape == (mda_universe.atoms.n_atoms, 3)
        np.testing.assert_allclose(ts.positions, coords, atol=1e-3, err_msg="Coordinates do not match the expected values")

    def test_reader_metadata(self, mda_universe):
        """Ensure ts.data is correctly populated and no bleeding between frames happen"""
        frame_0 = mda_universe.trajectory[0]

        assert 'pcoord' in frame_0.data
        assert 'coord' in frame_0.data

        assert 'iteration' in frame_0.data
        assert 'weight' in frame_0.data
        cputime_0 = frame_0.data['cputime']

        frame_last = mda_universe.trajectory[-1]
        assert 'weight' in frame_last.data
        assert 'parent_id' in frame_last.data
        cputime_last = frame_last.data['cputime']

        # Check bleeding using cputime since if more than 1 walker exists, the cputime should be different
        assert cputime_0 != cputime_last, "Data bleed occurred, cputime did not update"

    def test_reader_is_picklable(self, mda_universe):
        """Ensure Reader can be pickled and unpickled, which is required for MDAnalysis parallel backends"""
        import pickle

        reader = mda_universe.trajectory
        reader_copy = pickle.loads(pickle.dumps(reader))
        assert reader_copy.n_frames == reader.n_frames, "Unpickled reader has wrong frame count"
        assert reader_copy.n_atoms == reader.n_atoms, "Unpickled reader has wrong atom count"

    def test_parallel_rmsd_matches_serial(self):
        """Check RMSD computed with multiprocessing backend is numerically identical to serial"""
        from MDAnalysis.analysis import rms

        u1 = mda.Universe(hdf5_file, format='WESTPA')
        R_serial = rms.RMSD(u1, u1, select="index 0:50")
        R_serial.run(backend='serial')
        u1.trajectory.close()

        u2 = mda.Universe(hdf5_file, format='WESTPA')
        R_par = rms.RMSD(u2, u2, select="index 0:50")
        R_par.run(backend='multiprocessing', n_workers=2)
        u2.trajectory.close()

        np.testing.assert_allclose(
            R_serial.results.rmsd, R_par.results.rmsd, atol=1e-5, err_msg="Serial and parallel RMSD results differ"
        )

    def test_reader_pbc_dimensions(self, mda_universe):
        """Ensures that periodic boundary conditions (cell dimensions and angles) are correctly extracted, scaled (nm to Å) and usable by MDAnalysis tools"""
        from MDAnalysis.lib.distances import distance_array

        ts = mda_universe.trajectory[0]
        iter_num, seg_idx, actual_pos, path, local_frame = mda_universe.trajectory.frame_index[0]

        # Read raw data
        with h5py.File(path, 'r') as f:
            if 'cell_lengths' in f and 'cell_angles' in f:
                raw_lengths = f['cell_lengths'][actual_pos]
                raw_angles = f['cell_angles'][actual_pos]

                # Expected
                expected_dimensions = np.array(
                    [
                        raw_lengths[0] * 10.0,
                        raw_lengths[1] * 10.0,
                        raw_lengths[2] * 10.0,
                        raw_angles[0],
                        raw_angles[1],
                        raw_angles[2],
                    ],
                    dtype=np.float32,
                )
            else:
                expected_dimensions = np.zeros(6, dtype=np.float32)

        assert ts.dimensions is not None
        assert ts.dimensions.shape == (6,)
        assert ts.dimensions.dtype == np.float32

        np.testing.assert_allclose(
            ts.dimensions,
            expected_dimensions,
            atol=1e-4,
            err_msg="Box dimensions do not match expected",
        )

        if expected_dimensions[0] > 0:

            assert ts.volume > 0

            positions = mda_universe.atoms.positions
            pbc_distances = distance_array(positions, positions, box=ts.dimensions)

            # Verify the distance array computed successfully
            assert pbc_distances.shape == (mda_universe.atoms.n_atoms, mda_universe.atoms.n_atoms)
            assert np.min(pbc_distances) >= 0.0

    def test_save_auxdata_roundtrip(self, mda_universe, tmp_path):
        """Ensures that analysis data can be generated from the universe, saved back to a WESTPA HDF5 file as auxdata"""
        import shutil
        from MDAnalysis.analysis import rms
        from westpa.core.mdacrawl import save_to_west_h5

        # Temp file in order to not change the actual ref file
        test_h5 = tmp_path / "test_west_mdacrawl.h5"
        shutil.copy(hdf5_file, test_h5)

        R = rms.RMSD(mda_universe, mda_universe, select="resname Na+ or resname Cl-", ref_frame=0)
        R.run()

        flat_results = R.results.rmsd[:, 2]

        dataset_name = "test_rmsd"
        save_to_west_h5(mda_universe, flat_results, dataset_name, west_h5_path=str(test_h5), overwrite=True)

        with h5py.File(test_h5, 'r+') as f:
            iter_prec = getattr(mda_universe.trajectory, 'iter_prec', 8)
            first_iter = list(mda_universe.trajectory.frame_index)[0][0]
            iter_name = f'iter_{first_iter:0{iter_prec}d}'
            iter_group = f[f'iterations/{iter_name}']

            # Verify the dataset was created in the correct place
            assert 'auxdata' in iter_group, "auxdata group was not created!"
            assert dataset_name in iter_group['auxdata'], f"{dataset_name} was not saved!"

            saved_data = iter_group['auxdata'][dataset_name][:]
            assert len(saved_data.shape) == 2, "Data was not reshaped to 2D for WESTPA!"

            # Verify the values map correctly for the very first frame
            flat_idx = 0
            test_iter, test_seg, actual_pos, test_path, test_local_frame = mda_universe.trajectory.frame_index[flat_idx]
            expected_val = flat_results[flat_idx]
            actual_val = iter_group['auxdata'][dataset_name][test_seg, test_local_frame]

            assert np.isclose(expected_val, actual_val, equal_nan=True), "Saved HDF5 data does not match the results!"

            del iter_group['auxdata'][dataset_name]
            assert dataset_name not in iter_group['auxdata'], "Dataset was not successfully deleted!"

    def test_saved_auxdata_compatability(self, mda_universe, tmp_path):
        """Ensures that saved auxdata is compatible with other tools"""
        import shutil
        import subprocess
        import h5py
        from MDAnalysis.analysis import rms
        from westpa.core.mdacrawl import save_to_west_h5

        test_h5 = tmp_path / "test_west_mdacrawl.h5"
        shutil.copy(hdf5_file, test_h5)

        R = rms.RMSD(mda_universe, mda_universe, select="resname Na+ or resname Cl-", ref_frame=0)
        R.run()

        # Slice as [:, 2:3] instead of [:, 2] so that save_to_west_h5() is forced to build (Segment, Frames, 1), satisfying w_assign
        flat_results = R.results.rmsd[:, 2:3]

        dataset_name = "test_rmsd"
        save_to_west_h5(mda_universe, flat_results, dataset_name, west_h5_path=str(test_h5), overwrite=True)

        state_yaml = tmp_path / "states.yaml"
        state_yaml.write_text("""
---
states:
  - label: A
    coords:
      - [0.0]
  - label: B
    coords:
      - [10.0]
        """)
        assign_h5 = tmp_path / "assign.h5"

        cmd = [
            "w_assign",
            "-W",
            str(test_h5),
            "--dsspecs",
            "auxdata/test_rmsd",
            "--bins-from-expr",
            "[[0.0, 5.0, 10.0, 100.0, inf]]",
            "--states-from-file",
            str(state_yaml),
            "-o",
            str(assign_h5),
        ]

        try:
            # check=True makes the test to instantly fail if w_assign crashes
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)  # noqa
        except subprocess.CalledProcessError as e:
            import pytest

            pytest.fail(f"w_assign failed to process the auxdata!\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")

        assert assign_h5.exists(), "w_assign did not generate the output file!"

        with h5py.File(assign_h5, 'r') as f:
            assert 'assignments' in f, "The 'assignments' dataset was not found in assign.h5"
            assignments = f['assignments'][:]

            assert assignments.size > 0, "The w_assign assignments array is empty"
            assert len(assignments.shape) == 3, f"Expected a 3D array from w_assign, got {len(assignments.shape)}D"
