'''Parser and Reader to expose WESTPA simulation data as MDAnalysis Universe'''

import json
import warnings

import h5py
import numpy as np

try:
    from MDAnalysis.topology.base import TopologyReaderBase

    MDATopologyBase = TopologyReaderBase
except ImportError:
    MDATopologyBase = object


class WESTPAParser(MDATopologyBase):
    format = 'WESTPA'

    def parse(self, **kwargs):
        with h5py.File(self.filename, 'r') as f:
            iter_prec = f.attrs['west_iter_prec']
            first_iter_name = f'iter_{1:0{iter_prec}d}'

            iter_group = f[f'iterations/{first_iter_name}']
            if 'trajectories' in iter_group:
                traj_filename = iter_group['trajectories'].file.filename
            else:
                raise ValueError(f"Could not find 'trajectories' link in {self.filename}")

        with h5py.File(traj_filename, 'r') as f:
            topo_raw = f['topology'][()]

        if isinstance(topo_raw, np.ndarray):
            topo_str = topo_raw.item().decode('utf-8')
        elif isinstance(topo_raw, bytes):
            topo_str = topo_raw.decode('utf-8')
        else:
            topo_str = str(topo_raw)

        topo_str = topo_str.strip()
        if topo_str[0] == '{':
            return self._parse_json(topo_str)
        else:
            raise ValueError("Unknown topology format inside HDF5.")

    def _parse_json(self, topo_str):
        from MDAnalysis.core.topology import Topology
        from MDAnalysis.core.topologyattrs import Atomnames, Atomids, Resids, Resnames, Elements, Segids, Masses, Bonds
        from MDAnalysis.guesser.tables import masses as mass_table

        data = json.loads(topo_str)
        atom_names = []
        elements = []
        resnames = []
        resids = []
        atom_resindex = []
        residue_segindex = []
        bonds = []

        res_idx = 0
        seg_idx = 0
        for chain in data['chains']:
            for residue in chain['residues']:
                resnames.append(residue['name'])
                resids.append(residue['resSeq'])
                residue_segindex.append(seg_idx)
                for atom in residue['atoms']:
                    atom_names.append(atom['name'])
                    elements.append(atom['element'])
                    atom_resindex.append(res_idx)
                res_idx += 1
            seg_idx += 1

        for bond in data['bonds']:
            bonds.append(tuple(bond))

        n_atoms = len(atom_names)
        n_res = len(resnames)
        n_seg = seg_idx
        # Try both uppercase and capitalized since MDAnalysis's internal dictionary could be inconsistent (checked)
        mass_values = np.array([mass_table.get(e.upper(), mass_table.get(e.capitalize(), 0.0)) for e in elements], dtype=np.float64)

        return Topology(
            n_atoms=n_atoms,
            n_res=n_res,
            n_seg=n_seg,
            attrs=[
                Atomnames(np.array(atom_names, dtype=object)),
                Atomids(np.arange(n_atoms, dtype=np.int32)),
                Resids(np.array(resids, dtype=np.int32)),
                Resnames(np.array(resnames, dtype=object)),
                Elements(np.array(elements, dtype=object)),
                Segids(np.array([str(i) for i in range(n_seg)], dtype=object)),
                Masses(mass_values),
                Bonds(bonds),
            ],
            atom_resindex=np.array(atom_resindex, dtype=np.int32),
            residue_segindex=np.array(residue_segindex, dtype=np.int32),
        )


try:
    from MDAnalysis.coordinates.base import ReaderBase

    MDAReaderBase = ReaderBase
except ImportError:
    MDAReaderBase = object


class WESTPAReader(MDAReaderBase):
    format = 'WESTPA'

    def __init__(self, filename, n_atoms=None, **kwargs):
        self._n_atoms = n_atoms
        super().__init__(filename, **kwargs)
        self.filename = filename

        self._h5 = h5py.File(filename, 'r')
        # Will be used for caching later
        self._current_h5 = None
        self._current_path = None

        try:
            self._build_frame_index()
        except Exception:
            self._h5.close()
            raise

    def _build_frame_index(self):
        self.iter_prec = self._h5.attrs['west_iter_prec']

        self.frame_index = []
        n_iters = np.count_nonzero([self._h5['summary']['walltime'][:] > 0])

        for i in range(1, n_iters + 1):
            iter_name = f'iter_{i:0{self.iter_prec}d}'
            iter_group = self._h5[f'iterations/{iter_name}']

            if 'trajectories' not in iter_group:
                warnings.warn(
                    f"No trajectory group found in {iter_name}"
                )  # Sometimes the simulation ends before the trajectory file is written
                break

            seg_idx = iter_group['seg_index'][:]
            n_segs = len(seg_idx)

            traj_file = iter_group['trajectories']

            ptr = traj_file['pointer'][:]

            for seg_id in range(n_segs):
                # Filter using pointer dataset
                valid = ptr[:, 0] > 0
                seg_mask = ptr[:, 1] == seg_id
                actual_positions = np.where(valid & seg_mask)[0]
                for actual_pos in actual_positions:
                    self.frame_index.append((i, seg_id, actual_pos, traj_file.file.filename))

            self.ts = self._Timestep(self.n_atoms)

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def n_frames(self):
        return len(self.frame_index)

    def _read_frame(self, i):
        iter_num, seg_idx, actual_pos, path = self.frame_index[i]

        # Cache file handles to prevent massive file open/close
        if self._current_path != path:
            if self._current_h5:
                self._current_h5.close()
            self._current_h5 = h5py.File(path, 'r')
            self._current_path = path

        coords = self._current_h5['coordinates'][actual_pos] * 10.0  # MDTraj normalizes units to nm but MDAnalysis uses Ångströms
        self.ts.positions = coords.astype(np.float32)
        self.ts.frame = i
        self.ts.data.clear()  # Prevents bleeding of ts.data to each frame

        # Periodic boundary box info
        if 'cell_lengths' in self._current_h5 and 'cell_angles' in self._current_h5:
            lengths = self._current_h5['cell_lengths'][actual_pos] * 10.0
            angles = self._current_h5['cell_angles'][actual_pos]

            self.ts.dimensions = np.array([lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]], dtype=np.float32)
        else:
            # Fallback if no PBC data exists
            self.ts.dimensions = np.zeros(6, dtype=np.float32)

        # Energies
        if 'kineticEnergy' in self._current_h5:
            self.ts.data['kinetic_energy'] = self._current_h5['kineticEnergy'][actual_pos]
        if 'potentialEnergy' in self._current_h5:
            self.ts.data['potential_energy'] = self._current_h5['potentialEnergy'][actual_pos]

        # WESTPA specific metadata
        iter_name = f'iter_{iter_num:0{self.iter_prec}d}'
        si = self._h5[f'iterations/{iter_name}/seg_index'][seg_idx]
        self.ts.data['iteration'] = iter_num
        for name in si.dtype.names:
            if name not in self.ts.data:
                self.ts.data[name] = si[name]

        return self.ts

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        return self._read_frame(self.ts.frame + 1)

    def _reopen(self):
        self.ts.frame = -1

    def close(self):
        if self._current_h5 is not None:
            self._current_h5.close()
            self._current_h5 = None
        if hasattr(self, '_h5') and self._h5 is not None:
            self._h5.close()
            self._h5 = None

    # Parallelization support
    def __getstate__(self):
        if self._current_h5 is not None:
            self._current_h5.close()
        if self._h5 is not None:
            self._h5.close()

        state = self.__dict__.copy()
        state['_h5'] = None
        state['_current_h5'] = None
        state['_current_path'] = None
        return state


def __setstate__(self, state):
    self.__dict__.update(state)
    self._h5 = h5py.File(self.filename, 'r')
    self._current_h5 = None
    self._current_path = None
