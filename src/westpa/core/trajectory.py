import numpy as np
import os

from mdtraj import Trajectory, load as load_traj, FormatRegistry, formats as mdformats
from mdtraj.core.trajectory import _TOPOLOGY_EXTS, _get_extension as get_extension

FormatRegistry.loaders['.rst'] = mdformats.amberrst.load_restrt
FormatRegistry.fileobjects['.rst'] = mdformats.AmberRestartFile

TRAJECTORY_EXTS = list(FormatRegistry.loaders.keys())
TOPOLOGY_EXTS = list(_TOPOLOGY_EXTS)
for ext in [".h5", ".hdf5", ".lh5"]:
    TOPOLOGY_EXTS.remove(ext)


class WESTTrajectory(Trajectory):
    '''A subclass of ``mdtraj.Trajectory`` that contains the trajectory of atom coordinates with
    pointers denoting the iteration number and segment index of each frame.'''

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
        "Set the iteration index corresponding to each frame"

        self._iters = self._check_labels(value)
        self._shape = None

    @property
    def seg_labels(self):
        """Segment index corresponding to each frame

        Returns
        -------
        time : np.ndarray, shape=(n_frames,)
            The segment index corresponding to each frame
        """
        return self._segs

    @seg_labels.setter
    def seg_labels(self, value):
        "Set the segment index corresponding to each frame"

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
        """Join two ``Trajectory``s. This overrides ``mdtraj.Trajectory.join``
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
        """Slice the ``Trajectory``. This overrides ``mdtraj.Trajectory.slice``
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


def load_trajectory(folder):
    '''Load trajectory from ``folder`` using ``mdtraj`` and return a ``mdtraj.Trajectory``
    object. The folder should contain a trajectory and a topology file (with a recognizable
    extension) that is supported by ``mdtraj``. The topology file is optional if the
    trajectory file contains topology data (e.g., HDF5 format).
    '''
    traj_file = top_file = None
    for filename in os.listdir(folder):
        filepath = os.path.join(folder, filename)
        if not os.path.isfile(filepath):
            continue

        ext = get_extension(filename).lower()
        if ext in TOPOLOGY_EXTS and top_file is None:
            top_file = filename
        elif ext in TRAJECTORY_EXTS and traj_file is None:
            traj_file = filename

        if top_file is not None and traj_file is not None:
            break

    if traj_file is None:
        raise ValueError('trajectory file not found')

    traj_file = os.path.join(folder, traj_file)

    kwargs = {}
    if top_file is not None:
        top_file = os.path.join(folder, top_file)
        kwargs['top'] = top_file

    traj = load_traj(traj_file, **kwargs)
    return traj
