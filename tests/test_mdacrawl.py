import pytest
import os

import numpy as np

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
