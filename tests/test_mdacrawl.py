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
        iter_num, seg_idx, actual_pos, path = mda_universe.trajectory.frame_index[0]

        with h5py.File(path, 'r') as f:
            coords = f['coordinates'][actual_pos] * 10

        assert ts.positions is not None
        assert ts.positions.dtype == np.float32
        assert ts.positions.shape == (mda_universe.atoms.n_atoms, 3)
        np.testing.assert_allclose(ts.positions, coords, atol=1e-3, err_msg="Coordinates do not match the expected values")

    def test_reader_metadata(self, mda_universe):
        """Ensure ts.data is correctly populated and no bleeding between frames happen"""
        frame_0 = mda_universe.trajectory[0]
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
