'''Parser and Reader to expose WESTPA simulation data as MDAnalysis Universe'''

import json

import h5py
import numpy as np

from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import Atomnames, Atomids, Resids, Resnames, Elements, Segids, Masses, Bonds
from MDAnalysis.guesser.tables import masses as mass_table


class WESTPAParser(TopologyReaderBase):
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
