import numpy as np
import os
from functools import cache

from mdtraj import Trajectory


def parseResidueAtoms(residue, map):
    """Parse all atoms from residue. Taken from OpenMM 8.2.0."""

    for atom in residue.findall('Atom'):
        name = atom.attrib['name']
        for id in atom.attrib:
            map[atom.attrib[id]] = name


@cache
def loadNameReplacementTables():
    """Load the list of atom and residue name replacements. Taken from OpenMM 8.2.0."""

    # importing things here because they're only used in this function
    try:
        from importlib.resources import files
    except ImportError:
        from importlib_resources import files

    import xml.etree.ElementTree as etree
    from copy import copy

    residueNameReplacements = {}
    atomNameReplacements = {}

    # This XML file is a to map all sorts of atom names/ residue names to the PDB 3.0 convention.
    tree = etree.parse(files('westpa') / 'data/pdbNames.xml')
    allResidues = {}
    proteinResidues = {}
    nucleicAcidResidues = {}
    for residue in tree.getroot().findall('Residue'):
        name = residue.attrib['name']
        if name == 'All':
            parseResidueAtoms(residue, allResidues)
        elif name == 'Protein':
            parseResidueAtoms(residue, proteinResidues)
        elif name == 'Nucleic':
            parseResidueAtoms(residue, nucleicAcidResidues)
    for atom in allResidues:
        proteinResidues[atom] = allResidues[atom]
        nucleicAcidResidues[atom] = allResidues[atom]
    for residue in tree.getroot().findall('Residue'):
        name = residue.attrib['name']
        for id in residue.attrib:
            if id == 'name' or id.startswith('alt'):
                residueNameReplacements[residue.attrib[id]] = name
        if 'type' not in residue.attrib:
            atoms = copy(allResidues)
        elif residue.attrib['type'] == 'Protein':
            atoms = copy(proteinResidues)
        elif residue.attrib['type'] == 'Nucleic':
            atoms = copy(nucleicAcidResidues)
        else:
            atoms = copy(allResidues)
        parseResidueAtoms(residue, atoms)
        atomNameReplacements[name] = atoms

    return residueNameReplacements, atomNameReplacements


def convert_mdanalysis_top_to_mdtraj(universe):
    """Convert a MDAnalysis Universe object's topology to a ``mdtraj.Topology`` object."""

    from mdtraj import Topology
    from mdtraj.core.element import get_by_symbol
    from MDAnalysis.exceptions import NoDataError

    top = Topology()  # Empty topology object
    residueNameReplacements, atomNameReplacements = loadNameReplacementTables()

    # Add in all the chains (called segments in MDAnalysis)
    for chain_segment in universe.segments:
        top.add_chain()

    all_chains = list(top.chains)

    # Add in all the residues
    for residue in universe.residues:
        try:
            resname = residueNameReplacements[residue.resname]
        except KeyError:
            resname = residue.resname

        top.add_residue(name=resname, chain=all_chains[residue.segindex], resSeq=residue.resid)

    all_residues = list(top.residues)

    # Add in all the atoms
    for atom, resid in zip(universe.atoms, universe.atoms.resindices):
        try:
            atomname = residueNameReplacements[atom.resname][atom.name]
        except (KeyError, TypeError):
            atomname = atom.name

        top.add_atom(name=atomname, element=get_by_symbol(atom.element), residue=all_residues[resid])

    all_atoms = list(top.atoms)

    # Add in all the bonds.  Depending on the topology type (e.g., pdb), there might not be bond information.
    try:
        for b_idx in universe.bonds._bix:
            top.add_bond(all_atoms[b_idx[0]], all_atoms[b_idx[1]])
    except NoDataError:
        top.create_standard_bonds()

    return top


class WESTTrajectory(Trajectory):
    """A subclass of ``mdtraj.Trajectory`` that contains the trajectory of atom coordinates with
    pointers denoting the iteration number and segment index of each frame."""

    def __init__(
        self,
        coordinates,
        topology=None,
        time=None,
        iter_labels=None,
        seg_labels=None,
        pcoords=None,
        parent_ids=None,
        unitcell_lengths=None,
        unitcell_angles=None,
    ):
        if isinstance(coordinates, Trajectory):
            xyz = coordinates.xyz
            topology = coordinates.topology if topology is None else topology
            time = coordinates.time if time is None else time
            unitcell_lengths = coordinates.unitcell_lengths if unitcell_lengths is None else unitcell_lengths
            unitcell_angles = coordinates.unitcell_angles if unitcell_angles is None else unitcell_angles
        else:
            xyz = coordinates

        super(WESTTrajectory, self).__init__(xyz, topology, time, unitcell_lengths, unitcell_angles)
        self._shape = None
        self.iter_labels = iter_labels
        self.seg_labels = seg_labels
        self.pcoords = pcoords
        self.parent_ids = parent_ids

    def _string_summary_basic(self):
        """Basic summary of WESTTrajectory in string form."""

        unitcell_str = 'and unitcells' if self._have_unitcell else 'without unitcells'
        value = "%s with %d frames, %d atoms, %d residues, %s" % (
            self.__class__.__name__,
            self.n_frames,
            self.n_atoms,
            self.n_residues,
            unitcell_str,
        )
        return value

    def _check_labels(self, value):
        if value is None:
            value = 0
        elif isinstance(value, list):
            value = np.array(value)

        if np.isscalar(value):
            value = np.array([value] * self.n_frames, dtype=int)
        elif value.shape != (self.n_frames,):
            raise ValueError('Wrong shape. Got %s, should be %s' % (value.shape, (self.n_frames,)))

        return value

    def _check_pcoords(self, value):
        if value is None:
            value = 0.0
        elif isinstance(value, list):
            value = np.array(value)

        if np.isscalar(value):
            value = np.array([(value,)] * self.n_frames, dtype=float)

        if value.ndim == 1:
            value = np.tile(value, (self.n_frames, 1))
        elif value.ndim != 2:
            raise ValueError('pcoords must be a 2-D array')

        elif value.shape[0] != self.n_frames:
            raise ValueError('Wrong length. Got %s, should be %s' % (value.shape[0], self.n_frames))

        return value

    def iter_label_values(self):
        visited_ids = []

        for i in self.iter_labels:
            if i in visited_ids:
                continue
            visited_ids.append(i)
            yield i

    def seg_label_values(self, iteration=None):
        seg_labels = self.seg_labels[self.iter_labels == iteration]
        visited_ids = []

        for j in seg_labels:
            if j in visited_ids:
                continue
            visited_ids.append(j)
            yield j

    @property
    def label_values(self):
        for i in self.iter_label_values():
            for j in self.seg_label_values(i):
                yield i, j

    def _iter_blocks(self):
        for i, j in self.label_values:
            IandJ = np.logical_and(self.iter_labels == i, self.seg_labels == j)
            yield i, j, IandJ

    @property
    def iter_labels(self):
        """Iteration index corresponding to each frame

        Returns
        -------
        time : np.ndarray, shape=(n_frames,)
            The iteration index corresponding to each frame
        """

        return self._iters

    @iter_labels.setter
    def iter_labels(self, value):
        """Set the iteration index corresponding to each frame"""

        self._iters = self._check_labels(value)
        self._shape = None

    @property
    def seg_labels(self):
        """
        Segment index corresponding to each frame

        Returns
        -------
        time : np.ndarray, shape=(n_frames,)
            The segment index corresponding to each frame
        """

        return self._segs

    @seg_labels.setter
    def seg_labels(self, value):
        """Set the segment index corresponding to each frame"""

        self._segs = self._check_labels(value)
        self._shape = None

    @property
    def pcoords(self):
        return self._pcoords

    @pcoords.setter
    def pcoords(self, value):
        self._pcoords = self._check_pcoords(value)

    @property
    def parent_ids(self):
        return self._parent_ids

    @parent_ids.setter
    def parent_ids(self, value):
        self._parent_ids = self._check_labels(value)

    def join(self, other, check_topology=True, discard_overlapping_frames=False):
        """
        Join two ``Trajectory``s. This overrides ``mdtraj.Trajectory.join``
        so that it also handles WESTPA pointers.
        ``mdtraj.Trajectory.join``'s documentation for more details.
        """

        if isinstance(other, Trajectory):
            other = [other]

        new_traj = super(WESTTrajectory, self).join(
            other, check_topology=check_topology, discard_overlapping_frames=discard_overlapping_frames
        )

        trajectories = [self] + other
        if discard_overlapping_frames:
            for i in range(len(trajectories) - 1):
                x0 = trajectories[i].xyz[-1]
                x1 = trajectories[i + 1].xyz[0]

                if np.all(np.abs(x1 - x0) < 2e-3):
                    trajectories[i] = trajectories[i][:-1]

        iter_labels = []
        seg_labels = []
        parent_ids = []
        pshape = self.pcoords.shape
        pcoords = []

        for t in trajectories:
            if hasattr(t, "iter_labels"):
                iters = t.iter_labels
            else:
                iters = np.zeros(len(t)) - 1  # default iter label: -1

            iter_labels.append(iters)

            if hasattr(t, "seg_labels"):
                segs = t.seg_labels
            else:
                segs = np.zeros(len(t)) - 1  # default seg label: -1

            seg_labels.append(segs)

            if hasattr(t, "parent_ids"):
                pids = t.parent_ids
            else:
                pids = np.zeros(len(t)) - 1  # default parent_id: -1

            parent_ids.append(pids)

            if hasattr(t, "pcoords"):
                p = t.pcoords
            else:
                p = np.zeros((len(t), pshape[-1]), dtype=float)  # default pcoord: 0.0

            pcoords.append(p)

        iter_labels = np.concatenate(iter_labels)
        seg_labels = np.concatenate(seg_labels)
        parent_ids = np.concatenate(parent_ids)
        pcoords = np.concatenate(pcoords)

        new_westpa_traj = WESTTrajectory(
            new_traj, iter_labels=iter_labels, seg_labels=seg_labels, pcoords=pcoords, parent_ids=parent_ids
        )

        return new_westpa_traj

    def slice(self, key, copy=True):
        """
        Slice the ``Trajectory``. This overrides ``mdtraj.Trajectory.slice``
        so that it also handles WESTPA pointers. Please see
        ``mdtraj.Trajectory.slice``'s documentation for more details.
        """

        if isinstance(key, tuple):
            if self._shape is None:
                uniq_iters = np.unique(self.iter_labels)
                max_iter = uniq_iters.max()
                max_seg = self.seg_labels.max()
                max_n_trajs = 0
                for _, _, block in self._iter_blocks():
                    n_trajs = block.sum()
                    if n_trajs > max_n_trajs:
                        max_n_trajs = n_trajs

                self._shape = (max_iter, max_seg, max_n_trajs)
            else:
                max_iter, max_seg, max_n_trajs = self._shape

            M = np.full((max_iter + 1, max_seg + 1, max_n_trajs), -1, dtype=int)
            all_traj_indices = np.arange(self.n_frames, dtype=int)
            for i, j, block in self._iter_blocks():
                traj_indices = all_traj_indices[block]

                for k, traj_idx in enumerate(traj_indices):
                    M[i, j, k] = traj_idx

            selected_indices = M[key].flatten()
            if np.isscalar(selected_indices):
                selected_indices = np.array([selected_indices])
            key = selected_indices[selected_indices != -1]

        iters = self.iter_labels[key]
        segs = self.seg_labels[key]
        pcoords = self.pcoords[key, :]
        parent_ids = self.parent_ids[key]

        traj = super(WESTTrajectory, self).slice(key, copy)
        traj.iter_labels = iters
        traj.seg_labels = segs
        traj.pcoords = pcoords
        traj.parent_ids = parent_ids

        return traj


def get_extension(filename):
    """A function to get the format extension of a file."""

    (base, extension) = os.path.splitext(filename)

    # Return the other part of the extension as well if it's a gzip.
    if extension == '.gz':
        return os.path.splitext(base)[1] + extension

    return extension


def find_top_traj_file(folder, eligible_top, eligible_traj):
    """
    A general (reusable) function for identifying and returning the appropriate
    file names in ``folder`` which are toplogy and trajectory. Useful when writing custom loaders.
    Note that it's possible that the topology_file and trajectory_file are identical.

    Parameters
    ----------
    folder : str or os.Pathlike
        A string or Pathlike to the folder to search.

    eligible_top : list of strings
        A list of accepted topology file extensions.

    eligible_traj : list of strings
        A list of accepted trajectory file extensions.


    Returns
    -------
    top_file : str
        Path to topology file

    traj_file : str
        Path to trajectory file
    """

    # Setting up the return variables
    top_file = traj_file = None

    # Extract a list of all files, ignoring hidden files that start with a '.'
    file_list = [f_name for f_name in os.listdir(folder) if not f_name.startswith('.')]

    for filename in file_list:
        filepath = os.path.join(folder, filename)
        if not os.path.isfile(filepath):
            continue

        ext = get_extension(filename).lower()
        # Catching trajectory formats that can be topology and trajectories at the same time.
        # Only activates when there is a single file in the folder.
        if len(file_list) < 2 and ext in eligible_top and ext in eligible_traj:
            top_file = filename
            traj_file = filename

        # Assuming topology file is copied first.
        if ext in eligible_top and top_file is None:
            top_file = filename
        elif ext in eligible_traj and traj_file is None:
            traj_file = filename

        if top_file is not None and traj_file is not None:
            break

    if traj_file is None:
        raise ValueError('trajectory file not found')

    traj_file = os.path.join(folder, traj_file)

    if top_file is not None:
        top_file = os.path.join(folder, top_file)

    return top_file, traj_file


@cache
def mdtraj_supported_extensions():
    from mdtraj import FormatRegistry, formats as mdformats
    from mdtraj.core.trajectory import _TOPOLOGY_EXTS

    FormatRegistry.loaders['.rst'] = mdformats.amberrst.load_restrt
    FormatRegistry.fileobjects['.rst'] = mdformats.AmberRestartFile
    FormatRegistry.loaders['.ncrst'] = mdformats.amberrst.load_ncrestrt
    FormatRegistry.fileobjects['.ncrst'] = mdformats.AmberRestartFile

    TRAJECTORY_EXTS = list(FormatRegistry.loaders.keys())
    TOPOLOGY_EXTS = list(_TOPOLOGY_EXTS)

    for ext in [".h5", ".hdf5", ".lh5"]:
        TOPOLOGY_EXTS.remove(ext)

    return TOPOLOGY_EXTS, TRAJECTORY_EXTS


@cache
def mdanalysis_supported_extensions():
    import MDAnalysis as mda

    TRAJECTORY_EXTS = [reader.format if isinstance(reader.format, list) else [reader.format] for reader in mda._READERS.values()]
    TRAJECTORY_EXTS = list(set(f'.{ext.lower()}' for ilist in TRAJECTORY_EXTS for ext in ilist))

    TOPOLOGY_EXTS = [parser.format if isinstance(parser.format, list) else [parser.format] for parser in mda._PARSERS.values()]
    TOPOLOGY_EXTS = list(set(f'.{ext.lower()}' for ilist in TOPOLOGY_EXTS for ext in ilist))

    return TOPOLOGY_EXTS, TRAJECTORY_EXTS


def load_mdtraj(folder):
    """
    Load trajectory from ``folder`` using ``mdtraj`` and return a ``mdtraj.Trajectory``
    object. The folder should contain a trajectory and a topology file (with a recognizable
    extension) that is supported by ``mdtraj``. The topology file is optional if the
    trajectory file contains topology data (e.g., HDF5 format).
    """

    from mdtraj import load as load_traj

    TOPOLOGY_EXTS, TRAJECTORY_EXTS = mdtraj_supported_extensions()

    top_file, traj_file = find_top_traj_file(folder, TOPOLOGY_EXTS, TRAJECTORY_EXTS)

    # MDTraj likes the (optional) topology part to be provided within a dictionary
    traj = load_traj(traj_file, **{'top': top_file})

    return traj


def load_netcdf(folder):
    """
    Load netcdf file from ``folder`` using ``scipy.io`` and return a ``mdtraj.Trajectory``
    object. The folder should contain a Amber trajectory file with extensions `.nc` or `.ncdf`.

    Note coordinates and box lengths are all divided by 10 to change from Angstroms to nanometers.
    """

    from scipy.io import netcdf_file

    _, traj_file = find_top_traj_file(folder, [], ['.nc', '.ncdf', '.ncrst'])

    # Extracting these datasets
    datasets = {'coordinates': None, 'cell_lengths': None, 'cell_angles': None, 'time': None}
    convert = ['coordinates', 'cell_lengths']  # Length-based datasets that need to be converted from Å to nm
    optional = ['cell_lengths', 'cell_angles']

    with netcdf_file(traj_file) as rootgrp:
        for key, val in datasets.items():
            if key in optional:
                pass
            elif key in convert and key in rootgrp.variables:
                datasets[key] = rootgrp.variables[key][()].copy() / 10  # From Å to nm
            else:
                datasets[key] = rootgrp.variables[key][()].copy()  # noqa: F841

    map_dataset = {
        'coordinates': datasets['coordinates'],
        'unitcell_lengths': datasets['cell_lengths'],
        'unitcell_angles': datasets['cell_angles'],
        'time': datasets['time'],
    }

    return WESTTrajectory(**map_dataset)


def load_mdanalysis(folder):
    """
    Load a file from ``folder`` using ``MDAnalysis`` and return a ``mdtraj.Trajectory``
    object. The folder should contain a trajectory and a topology file (with a recognizable
    extension) that is supported by ``MDAnalysis``. The topology file is optional if the
    trajectory file contains topology data (e.g., H5MD format).

    Note coordinates and box lengths are all divided by 10 to change from Angstroms to nanometers.
    """

    import MDAnalysis as mda

    TOPOLOGY_EXTS, TRAJECTORY_EXTS = mdanalysis_supported_extensions()

    top_file, traj_file = find_top_traj_file(folder, TOPOLOGY_EXTS, TRAJECTORY_EXTS)

    u = mda.Universe(top_file, traj_file)

    tot_frames = len(u.trajectory)
    coords = np.zeros((tot_frames, len(u.atoms), 3))
    time = np.zeros((tot_frames))

    # Periodic Boundary Conditions
    periodic = u.trajectory.periodic
    cell_lengths = np.zeros((tot_frames, 3)) if periodic else None
    cell_angles = np.zeros((tot_frames, 3)) if periodic else None

    for iframe, frame in enumerate(u.trajectory):
        coords[iframe] = frame._pos
        time[iframe] = frame.time

        if periodic:
            cell_lengths[iframe] = frame.dimensions[:3]
            cell_angles[iframe] = frame.dimensions[3:]

    # Length-based datasets that need to be converted
    convert = [coords, cell_lengths] if periodic else [coords]
    for dset in convert:
        dset = mda.units.convert(dset, 'angstrom', 'nanometer')

    traj = WESTTrajectory(coordinates=coords, unitcell_lengths=cell_lengths, unitcell_angles=cell_angles, time=time)

    return traj
