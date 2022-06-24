import concurrent.futures as cf
import functools
import inspect
import operator
import os

import mdtraj
import numpy as np

from tqdm import tqdm
from typing import Callable
from westpa.analysis.core import Walker, Trace
from westpa.core.states import InitialState
from westpa.core.h5io import WESTIterationFile


class Trajectory:
    """A callable that returns the trajectory of a walker or trace.

    Parameters
    ----------
    fget : callable
        Function for retrieving a single trajectory segment. Must take a
        :class:`Walker` instance as its first argument and accept a boolean
        keyword argument `include_initpoint`. The function should return a
        sequence (e.g., a list or ndarray) representing the trajectory of
        the walker. If `include_initpoint` is True, the trajectory segment
        should include its initial point. Otherwise, the trajectory segment
        should exclude its initial point.
    fconcat : callable, optional
        Function for concatenating trajectory segments. Must take a sequence
        of trajectory segments as input and return their concatenation. The
        default concatenation function is :func:`concatenate`.

    """

    def __init__(self, fget=None, *, fconcat=None):
        if fget is None:
            return functools.partial(self.__init__, fconcat=fconcat)

        if 'include_initpoint' not in inspect.signature(fget).parameters:
            raise ValueError("'fget' must accept a parameter 'include_initpoint'")

        self._fget = fget
        self.fconcat = fconcat

        self._segment_collector = SegmentCollector(self)

    @property
    def segment_collector(self):
        """SegmentCollector: Segment retrieval manager."""
        return self._segment_collector

    @property
    def fget(self):
        """callable: Function for getting trajectory segments."""
        return self._fget

    @property
    def fconcat(self):
        """callable: Function for concatenating trajectory segments."""
        return self._fconcat

    @fconcat.setter
    def fconcat(self, value):
        if value is None:
            value = concatenate
        elif not isinstance(value, Callable):
            raise TypeError("'fconcat' must be a callable object")
        self._fconcat = value

    def __call__(self, obj, include_initpoint=True, **kwargs):
        if isinstance(obj, Walker):
            value = self.fget(obj, include_initpoint=include_initpoint, **kwargs)
            self._validate_segment(value)
            return value
        if isinstance(obj, Trace):
            initpoint_mask = np.full(len(obj), False)
            initpoint_mask[0] = include_initpoint
            segments = self.segment_collector.get_segments(obj, initpoint_mask, **kwargs)
            return self.fconcat(segments)
        raise TypeError('argument must be a Walker or Trace instance')

    def _validate_segment(self, value):
        if not hasattr(value, '__getitem__'):
            msg = f"{type(value).__name__!r} object can't be concatenated"
            raise TypeError(msg)


class SegmentCollector:
    """An object that manages the retrieval of trajectory segments.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory to which the segment collector is attached.
    use_threads : bool, default False
        Whether to use a pool of threads to retrieve trajectory segments
        asynchronously. Setting this parameter to True may be may be
        useful when segment retrieval is an I/O bound task.
    max_workers : int, optional
        Maximum number of threads to use. The default value is specified in the
        `ThreadPoolExecutor <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor>`_
        documentation.
    show_progress : bool, default False
        Whether to show a progress bar when retrieving multiple segments.

    """

    def __init__(self, trajectory, use_threads=False, max_workers=None, show_progress=False):
        self.trajectory = trajectory
        self.use_threads = use_threads
        self.max_workers = max_workers
        self.show_progress = show_progress

    @property
    def trajectory(self):
        return self._trajectory

    @trajectory.setter
    def trajectory(self, value):
        if not isinstance(value, Trajectory):
            msg = f'trajectory must be an instance of {Trajectory}'
            raise TypeError(msg)
        self._trajectory = value

    @property
    def use_threads(self):
        return self._use_threads

    @use_threads.setter
    def use_threads(self, value):
        if not isinstance(value, bool):
            raise TypeError('use_threads must be True or False')
        self._use_threads = value

    @property
    def max_workers(self):
        return self._max_workers

    @max_workers.setter
    def max_workers(self, value):
        if value is None:
            self._max_workers = None
            return
        if value <= 0:
            raise ValueError('max_workers must be greater than 0')
        self._max_workers = value

    @property
    def show_progress(self):
        return self._show_progress

    @show_progress.setter
    def show_progress(self, value):
        if not isinstance(value, bool):
            raise ValueError('show_progress must be True or False')
        self._show_progress = value

    def get_segments(self, walkers, initpoint_mask=None, **kwargs):
        """Retrieve the trajectories of multiple walkers.

        Parameters
        ----------
        walkers : sequence of Walker
            The walkers for which to retrieve trajectories.
        initpoint_mask : sequence of bool, optional
            A Boolean mask indicating whether each trajectory segment should
            include (True) or exclude (False) its initial point. Default is
            all True.

        Returns
        -------
        list of sequences
            The trajectory of each walker.

        """
        if initpoint_mask is None:
            initpoint_mask = np.full(len(walkers), True)
        else:
            initpoint_mask = np.asarray(initpoint_mask, dtype=bool)

        get_segment = functools.partial(self.trajectory, **kwargs)

        tqdm_kwargs = dict(
            desc='Retrieving segments',
            disable=(not self.show_progress),
            position=0,
            total=len(walkers),
        )

        if self.use_threads:
            with cf.ThreadPoolExecutor(self.max_workers) as executor:
                future_to_key = {
                    executor.submit(get_segment, walker, include_initpoint=i): key
                    for key, (walker, i) in enumerate(zip(walkers, initpoint_mask))
                }
                futures = list(tqdm(cf.as_completed(future_to_key), **tqdm_kwargs))
                futures.sort(key=future_to_key.get)
                segments = (future.result() for future in futures)
        else:
            it = (get_segment(walker, include_initpoint=i) for walker, i in zip(walkers, initpoint_mask))
            segments = tqdm(it, **tqdm_kwargs)

        return list(segments)


class BasicMDTrajectory(Trajectory):
    """Trajectory reader for MD trajectories stored as in the
    `Basic Tutorial <https://github.com/westpa/westpa_tutorials/tree/main/basic_nacl>`_.

    Parameters
    ----------
    top : str or mdtraj.Topology, default 'bstate.pdb'
    traj_ext : str, default '.dcd'
    state_ext : str, default '.xml'
    sim_root : str, default '.'

    """

    def __init__(self, top='bstate.pdb', traj_ext='.dcd', state_ext='.xml', sim_root='.'):
        self.top = top
        self.traj_ext = traj_ext
        self.state_ext = state_ext
        self.sim_root = sim_root

        def fget(walker, include_initpoint=True, atom_indices=None, sim_root=None):
            sim_root = sim_root or self.sim_root

            if isinstance(self.top, str):
                top = os.path.join(sim_root, 'common_files', self.top)
            else:
                top = self.top

            path = os.path.join(
                sim_root,
                'traj_segs',
                format(walker.iteration.number, '06d'),
                format(walker.index, '06d'),
                'seg' + self.traj_ext,
            )
            if top is not None:
                traj = mdtraj.load(path, top=top)
            else:
                traj = mdtraj.load(path)

            if include_initpoint:
                parent = walker.parent

                if isinstance(parent, InitialState):
                    if parent.istate_type == InitialState.ISTATE_TYPE_BASIS:
                        path = os.path.join(
                            sim_root,
                            'bstates',
                            parent.basis_state.auxref,
                        )
                    else:
                        path = os.path.join(
                            sim_root,
                            'istates',
                            str(parent.iter_created),
                            str(parent.state_id) + self.state_ext,
                        )
                else:
                    path = os.path.join(
                        sim_root,
                        'traj_segs',
                        format(walker.iteration.number - 1, '06d'),
                        format(parent.index, '06d'),
                        'seg' + self.state_ext,
                    )

                frame = mdtraj.load(path, top=traj.top)
                traj = frame.join(traj, check_topology=False)

            if atom_indices is not None:
                traj.atom_slice(atom_indices, inplace=True)

            return traj

        super().__init__(fget)

        self.segment_collector.use_threads = True
        self.segment_collector.show_progress = True


class HDF5MDTrajectory(Trajectory):
    """Trajectory reader for MD trajectories stored by the HDF5 framework."""

    def __init__(self):
        def fget(walker, include_initpoint=True, atom_indices=None):
            iteration = walker.iteration

            try:
                link = iteration.h5group['trajectories']
            except KeyError:
                msg = 'the HDF5 framework does not appear to have been used to store trajectories for this run'
                raise ValueError(msg)

            with WESTIterationFile(link.file.filename) as traj_file:
                traj = traj_file.read_as_traj(
                    iteration=iteration.number,
                    segment=walker.index,
                    atom_indices=atom_indices,
                )

            if include_initpoint:
                parent = walker.parent

                if isinstance(parent, InitialState):
                    if parent.istate_type == InitialState.ISTATE_TYPE_BASIS:
                        link = walker.run.h5file.get_iter_group(0)['trajectories']
                        with WESTIterationFile(link.file.filename) as traj_file:
                            frame = traj_file.read_as_traj(
                                iteration=0,
                                segment=parent.basis_state_id,
                                atom_indices=atom_indices,
                            )
                    elif parent.istate_type == InitialState.ISTATE_TYPE_GENERATED:
                        link = walker.run.h5file.get_iter_group(0)['trajectories']
                        istate_iter = -int(parent.iter_created)  # the conversion to int is because iter_created is uint
                        with WESTIterationFile(link.file.filename) as traj_file:
                            frame = traj_file.read_as_traj(
                                iteration=istate_iter,
                                segment=parent.state_id,
                                atom_indices=atom_indices,
                            )
                    else:
                        raise ValueError('unsupported initial state type: %d' % parent.istate_type)
                else:
                    frame = fget(parent, include_initpoint=False, atom_indices=atom_indices)[-1]
                traj = frame.join(traj, check_topology=False)

            return traj

        super().__init__(fget)

        self.segment_collector.use_threads = False
        self.segment_collector.show_progress = True


def concatenate(segments):
    """Return the concatenation of a sequence of trajectory segments.

    Parameters
    ----------
    segments : sequence of sequences
        A sequence of trajectory segments.

    Returns
    -------
    sequence
        The concatenation of `segments`.

    """
    if isinstance(segments[0], np.ndarray):
        return np.concatenate(segments)
    if isinstance(segments[0], mdtraj.Trajectory):
        return segments[0].join(segments[1:], check_topology=False)
    return functools.reduce(operator.concat, segments)
