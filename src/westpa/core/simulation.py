import io
import itertools
import logging
import math
import operator
import os
import time
from datetime import timedelta

import numpy as np

from .state import State
from .segment import Segment
from .binning import NopMapper
from .resamplers import HuberKimResampler
from .source_sink import Source, Sink
from ._data_manager import DataManager
from ..work_managers import SerialWorkManager
from ..work_managers.core import WorkManager

logger = logging.getLogger(__name__)


# copied from https://docs.python.org/3/library/itertools.html#itertools.batched
# TODO: Replace with itertools.batched() when we require Python >=3.12.
def batched(iterable, n, *, strict=False):
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


def default_pcoord(segment):
    if (state := segment.initial_state.coord) is None:
        raise ValueError(f"can't use default progress coordinate: {state} doesn't have a 'coord' value")
    return np.stack((segment.initial_state.coord, segment.final_state.coord))


def _report_bin_statistics(bins):
    bin_counts = np.fromiter(map(len, bins), dtype=int, count=len(bins))
    bin_probs = np.fromiter(map(operator.attrgetter('weight'), bins), dtype=float, count=len(bins))

    min_bin_prob = bin_probs[bin_probs != 0].min()
    max_bin_prob = bin_probs.max()
    bin_drange = math.log(max_bin_prob / min_bin_prob)
    n_pop = np.count_nonzero(bin_counts)

    logger.info('{:d} of {:d} ({:%}) bins are populated'.format(n_pop, len(bins), n_pop / len(bins)))
    logger.info('minimum non-zero bin probability:       {:g}'.format(min_bin_prob))
    logger.info('maximum bin probability:                {:g}'.format(max_bin_prob))
    logger.info('bin probability dynamic range (kT):     {:g}'.format(bin_drange))


class Simulation:
    """Interface for initializing and running a weighted ensemble (WE) simulation.

    Parameters
    ----------
    datafile : path-like or io.BytesIO
        HDF5 file used to store simulation data.
    propagator : Propagator
        Routine for propagating trajectories forward in time.
        See the :class:`~westpa.Propagator` protocol for details.
    pcoord_calculator : callable, optional
        Routine that computes the progress coordinate time series for
        a given trajectory segment. It must accept arguments
        ``(segment, parent)`` and return either a 2-D array (``pcoord``) or a
        tuple of the form ``(pcoord, auxdata)``, where ``auxdata`` is a
        dictionary of named arrays to store as auxiliary data.
        If ``segment`` continues a trajectory, ``parent`` is the
        preceding segment; otherwise it is None.
        If `pcoord_calculator` is not specified, either `propagator` must set
        the ``pcoord`` segment attribute, or else the raw
        coordinates are used as progress coordinates, equivalent to passing::

            def default_pcoord(segment, parent):
                return numpy.stack(
                    (segment.initial_state.coord, segment.final_state.coord)
                )

    bin_mapper : BinMapper, optional
        Routine for assigning trajectories to bins. See the
        :class:`~westpa.BinMapper` protocol for details. By default, all the
        trajectories are grouped into a single bin.
    bin_target_counts : int or sequence of int, default 1
        Target number of trajectories (allocation) for each bin. If an integer
        is provided, the value is applied to all the bins. If a sequence
        is provided, its length must match the output length of `bin_mapper`.
    resampler : Resampler, optional
        Routine for resampling the trajectories in each bin. Defaults to
        ``HuberKimResampler()``.
    source : Source, optional
        Set of states from which to reinitiate trajectories that reach a sink.
        Must be provided together with `sink`. See the
        :meth:`configure_recycling` method for more details.
    sink : Sink or iterable of Sink, optional
        One or more sink (target) regions. Must be provided together with `source`.
    istate_generator : callable, optional
        Routine for modifying the source distribution on the fly. It must
        accept a state from `source` as input and return a new state.
    work_manager : WorkManager, optional
        Work manager for executing calls to `propagator`, `pcoord_calculator`, and
        `istate_generator`. By default, calls are executed serially.

    Attributes
    ----------
    datafile : path-like or io.BytesIO
        HDF5 data file.
    propagator : Propagator
        Propagator.
    pcoord_calculator : callable or None
        Progress coordinate calculator.
    bin_mapper : BinMapper
        Bin mapper.
    bin_target_counts : numpy.ndarray
        Bin target counts.
    resampler : Resampler
        Resampler.
    source : Source or None
        Source distribution.
    sinks : tuple of Sink
        Sink region(s).
    istate_generator : callable or None
        Initial state generator.
    work_manager : WorkManager
        Work manager.
    n_iter : int or None
        Current iteration number.
    initialized : bool
        Whether the simulation has been initialized.
    segments : list of Segment
        Current segment information.

    Methods
    -------
    initialize
    run
    configure_recycling
    disable_recycling
    register_callback

    """

    def __init__(
        self,
        datafile,
        propagator,
        pcoord_calculator=None,
        bin_mapper=None,
        bin_target_counts=1,
        resampler=None,
        source=None,
        sink=None,
        istate_generator=None,
        work_manager=None,
    ):
        self._data_manager = DataManager(datafile)

        self._propagator = None
        self._pcoord_calculator = None
        self._bin_mapper = None
        self._bin_target_counts = None
        self._resampler = None

        self._source = None
        self._sinks = ()
        self._istate_generator = None

        self._work_manager = None

        self._n_iter = None
        self._prev_iter_segments = []
        self._segments = []
        self._resampled_segments = []  # populated by _run_we()
        self._next_iter_segments = []  # populated by _prepare_new_iteration()

        self.propagator = propagator
        self.pcoord_calculator = pcoord_calculator
        self.bin_mapper = bin_mapper
        self.bin_target_counts = bin_target_counts
        self.resampler = resampler or HuberKimResampler()

        if source or sink:
            if not (source and sink):
                raise ValueError("'source' and 'sink' must be provided together")
            self.configure_recycling(source, sink, istate_generator)

        self.work_manager = work_manager or SerialWorkManager()

        if isinstance(datafile, io.BytesIO):
            initialized = bool(datafile.getbuffer().nbytes)
        else:
            initialized = os.path.exists(datafile)

        if initialized:  # open existing simulation
            self._data_manager.open_backing()
            self._n_iter = self._data_manager.current_iteration
            self._segments = self._data_manager.get_segments()
            if self._n_iter > 1:
                self._prev_iter_segments = self._data_manager.get_segments(self._n_iter - 1)
            self._data_manager.close_backing()

    @property
    def datafile(self):
        return self._data_manager.we_h5filename

    @property
    def propagator(self):
        return self._propagator

    @propagator.setter
    def propagator(self, value):
        if not callable(value):
            raise TypeError("'propagator' must be callable")
        self._propagator = value

    @property
    def pcoord_calculator(self):
        return self._pcoord_calculator

    @pcoord_calculator.setter
    def pcoord_calculator(self, value):
        if value is not None and not callable(value):
            raise TypeError("'pcoord_calculator' must be callable or None")
        self._pcoord_calculator = value

    @property
    def bin_mapper(self):
        return self._bin_mapper

    @bin_mapper.setter
    def bin_mapper(self, value):
        if value is None:
            value = NopMapper()
        elif not callable(value):
            raise TypeError("'bin_mapper' must be callable or None")
        self._bin_mapper = value

    @property
    def bin_target_counts(self):
        return self._bin_target_counts

    @bin_target_counts.setter
    def bin_target_counts(self, value):
        value = np.asarray(value, dtype=int)
        if value.ndim > 1:
            raise TypeError("'bin_target_counts' must be an integer or a sequence of integers")
        if (value < 1).any():
            raise ValueError("'bin_target_counts' must be positive")
        self._bin_target_counts = value.astype(np.min_scalar_type(value.max()))

    @property
    def resampler(self):
        return self._resampler

    @resampler.setter
    def resampler(self, value):
        if not callable(value):
            raise TypeError("'resampler' must be callable")
        self._resampler = value

    @property
    def source(self):
        return self._source

    @property
    def sinks(self):
        return self._sinks

    @property
    def istate_generator(self):
        return self._istate_generator

    @property
    def work_manager(self):
        return self._work_manager

    @work_manager.setter
    def work_manager(self, value):
        if not isinstance(value, WorkManager):
            raise TypeError("'work_manager' must be a WorkManager object")
        self._work_manager = value

    @property
    def n_iter(self):
        return self._n_iter

    @property
    def initialized(self):
        return self._n_iter is not None

    @property
    def segments(self):
        return self._segments

    def initialize(
        self,
        states,
        weights=None,
    ):
        """Initialize the simulation.

        Parameters
        ----------
        states : State or iterable of State
            States from which to initiate trajectories (one per trajectory).
        weights : sequence of float, optional
            Weight to assign each trajectory. Defaults to a uniform distribution.

        """
        if self.initialized:
            raise RuntimeError("can't initialize the simulation: already initialized")

        self._data_manager.prepare_backing()
        logger.info(f'Created HDF5 file {self.datafile!r}')

        if isinstance(states, State):
            states = [states]
        else:
            states = list(states)

        if weights is None:
            weights = np.ones(len(states))
        else:
            weights = np.array(weights, dtype=float)
            if (weights <= 0).any():
                raise ValueError("'weights' must be positive")
            if len(weights) != len(states):
                raise ValueError("length of 'weights' must match number of initial states")
        weights /= weights.sum()

        self._segments = [
            Segment(
                n_iter=1,
                weight=weight,
                parent_id=-1,
                wtg_parent_ids={-1},
                initial_state=state,
                status=Segment.Status.PREPARED,
            )
            for state, weight in zip(states, weights)
        ]

        self._data_manager.prepare_iteration(1, self._segments)
        self._data_manager.current_iteration = 1
        self._n_iter = 1

        logger.info('Simulation prepared.')
        self._report_segment_statistics()
        self._data_manager.flush_backing()
        self._data_manager.close_backing()

    def run(self, n_iters=1, max_walltime=None):
        """Run the simulation.

        Parameters
        ----------
        n_iters : int, default 1
            Number of iterations to run.
        max_walltime : float, optional
            Maximum wall-clock time in seconds. If provided, the simulation
            will be stopped if it is estimated that the next iteration would
            cause the total run time to exceed `max_walltime`.

        """
        if not self.initialized:
            raise RuntimeError('simulation must be initialized before calling run()')

        with self.work_manager as work_manager:
            if work_manager.is_master:
                work_manager.install_sigint_handler()
                self._run(n_iters, max_walltime)
            else:
                work_manager.run()

    def configure_recycling(self, source, sink, istate_generator=None):
        """Configure source-sink boundary conditions (recycling).

        Parameters
        ----------
        source : Source, optional
            Set of states from which to reinitiate trajectories that reach a
            sink.
        sink : Sink or iterable of Sink, optional
            One or more sink (target) regions.
        istate_generator : callable, optional
            Routine for modifying the source distribution on the fly. It must
            accept a state from `source` as input and return a new state.

        """
        if not isinstance(source, Source):
            raise TypeError("'source' must be a Source object")

        if isinstance(sink, Sink):
            sinks = (sink,)
        else:
            sinks = sink
            if not all(isinstance(item, Sink) for item in sinks):
                raise TypeError("'sink' must be a Sink object or an iterable of Sink objects")

        if istate_generator is not None and not callable(istate_generator):
            raise TypeError("'istate_generator' must be callable")

        self._source = source
        self._sinks = sinks
        self._istate_generator = istate_generator

    def disable_recycling(self):
        """Disable recycling."""
        self._source = None
        self._sinks = ()
        self._istate_generator = None

    def _run(self, n_iters, max_walltime):
        self._prepare_run()

        start_time = time.time()
        stop_time = None
        if max_walltime:
            stop_time = start_time + max_walltime
            logger.info(f'Maximum wallclock time: {timedelta(seconds=max_walltime or 0)}')

        max_iter = self._n_iter + n_iters - 1

        iter_elapsed = 0
        while self._n_iter <= max_iter:
            if max_walltime and time.time() + 1.1 * iter_elapsed >= stop_time:
                logger.info(f'Iteration {self._n_iter} would require more than the alloted time. Ending run.')
                return
            try:
                iter_start_time = time.time()

                logger.info('\n' + time.asctime())
                logger.info(f'Iteration {self._n_iter} (of {max_iter})')

                self._prepare_iteration()
                self._propagate()

                cputime = sum(segment.cputime for segment in self._segments)

                self._next_iteration()

                iter_elapsed = time.time() - iter_start_time
                iter_summary = self._data_manager.get_iter_summary(self._n_iter - 1)
                iter_summary['walltime'] += iter_elapsed
                iter_summary['cputime'] = cputime
                self._data_manager.update_iter_summary(iter_summary, self._n_iter - 1)

                walltime = timedelta(seconds=iter_summary['walltime'].item())
                cputime = timedelta(seconds=cputime)

                logger.info('Iteration completed successfully')
                logger.info(f'Iteration walltime: {walltime}' + (f', cputime: {cputime}' if cputime else ''))
            finally:
                self._data_manager.flush_backing()

        self._finalize_run()

        logger.info(time.asctime())
        logger.info('WESTPA run complete.')

    def _report_segment_statistics(self, save_summary=True):
        seg_probs = np.fromiter(
            map(operator.attrgetter('weight'), self._segments),
            dtype=float,
            count=len(self._segments),
        )
        norm = seg_probs.sum()

        min_seg_prob = seg_probs.min()
        max_seg_prob = seg_probs.max()
        seg_drange = np.log(max_seg_prob / min_seg_prob)

        eps = np.finfo(float).eps

        logger.info('number of segments:                         {:d}'.format(len(self._segments)))
        logger.info('per-segment minimum non-zero probability:   {:g}'.format(min_seg_prob))
        logger.info('per-segment maximum non-zero probability:   {:g}'.format(max_seg_prob))
        logger.info('per-segment probability dynamic range (kT): {:g}'.format(seg_drange))
        logger.info('norm = {:g}, error in norm = {:g} ({:.2g}*eps)'.format(norm, (norm - 1), (norm - 1) / eps))

        if min_seg_prob < 1e-100:
            logger.warning(
                'Minimum segment weight is < 1e-100 and might not be physically relevant. '
                'Please reconsider your progress coordinate or binning scheme.'
            )

        if save_summary:
            iter_summary = self._data_manager.get_iter_summary()
            iter_summary['n_particles'] = len(self._segments)
            iter_summary['norm'] = norm
            iter_summary['min_seg_prob'] = min_seg_prob
            iter_summary['max_seg_prob'] = max_seg_prob
            if np.isnan(iter_summary['cputime']):
                iter_summary['cputime'] = 0.0
            if np.isnan(iter_summary['walltime']):
                iter_summary['walltime'] = 0.0
            self._data_manager.update_iter_summary(iter_summary)

    def _prepare_run(self):
        self._data_manager.prepare_run()
        self._invoke_callbacks('prepare_run')

    def _finalize_run(self):
        self._invoke_callbacks('finalize_run')
        self._data_manager.finalize_run()

    def _prepare_iteration(self):
        self._invoke_callbacks('prepare_iteration')
        self._report_segment_statistics()

    def _finalize_iteration(self):
        self._data_manager.update_seg_index(self._n_iter, self._segments)
        self._data_manager.write_auxdata(self._n_iter, self._segments)

    def _propagate(self):
        if value := getattr(self.propagator, 'block_size', None):
            propagator_block_size = value
        else:
            propagator_block_size = len(self._segments)

        istate_futures = set()
        propagator_futures = set()
        pcoord_future_map = {}

        # partition segments by status
        unprepared_segments = []
        prepared_segments = []
        complete_segments = []
        for segment in self._segments:
            match segment.status:
                case Segment.Status.UNSET:
                    unprepared_segments.append(segment)
                case Segment.Status.PREPARED:
                    prepared_segments.append(segment)
                case Segment.Status.COMPLETE:
                    complete_segments.append(segment)

        # dispatch pending istate generation tasks
        if unprepared_segments:
            states = self.source.random_choice(len(unprepared_segments))
            if self.istate_generator is not None:
                for state in states:
                    future = self.work_manager.submit(self.istate_generator, args=(state,))
                    istate_futures.add(future)
            else:
                segments = []
                for state in states:
                    segment = unprepared_segments.pop()
                    segment.initial_state = state
                    segments.append(segment)
                self._data_manager.write_initial_states(self._n_iter, segments)
                prepared_segments += segments

        # dispatch pending propagation tasks
        for segments in batched(prepared_segments, propagator_block_size):
            future = self.work_manager.submit(self.propagator, args=(segments,))
            propagator_futures.add(future)
        prepared_segments.clear()

        # dispatch pending pcoord calculation tasks
        if segments := [s for s in complete_segments if s.pcoord is None]:
            pcoord_future_map |= self._calculate_pcoords(segments)

        logger.info('Waiting for segments to complete...')

        while futures := istate_futures | propagator_futures | pcoord_future_map.keys():
            future = self.work_manager.wait_any(futures)

            if future in istate_futures:
                istate_futures.remove(future)
                try:
                    state = future.get_result()
                except Exception as e:
                    raise RuntimeError("call to 'istate_generator' failed") from e

                segment = unprepared_segments.pop()
                segment.initial_state = state
                segment.status = Segment.Status.PREPARED
                prepared_segments.append(segment)

                self._segments[segment.seg_id] = segment
                self._data_manager.write_initial_states(self._n_iter, segments=[segment])

                if len(prepared_segments) == propagator_block_size or not istate_futures:
                    future = self.work_manager.submit(self.propagator, args=(prepared_segments,))
                    propagator_futures.add(future)
                    prepared_segments.clear()

            elif future in propagator_futures:
                propagator_futures.remove(future)
                try:
                    segments = future.get_result()
                except Exception as e:
                    raise RuntimeError("call to 'propagator' failed") from e

                for segment in segments:
                    segment.status = Segment.Status.COMPLETE
                    self._segments[segment.seg_id] = segment
                self._data_manager.write_final_states(self._n_iter, segments)

                pcoord_future_map |= self._calculate_pcoords(segments)

            elif future in pcoord_future_map:
                segment = pcoord_future_map.pop(future)
                try:
                    result = future.get_result()
                except Exception as e:
                    raise RuntimeError("call to 'pcoord_calculator' failed") from e

                if isinstance(result, tuple):
                    pcoord, auxdata = result
                    segment.pcoord = pcoord
                    segment.data.update(auxdata)
                else:
                    segment.pcoord = result

                self._segments[segment.seg_id] = segment
                self._data_manager.write_pcoords(self._n_iter, segments=[segment])

        self._data_manager.flush_backing()

    def _calculate_pcoords(self, segments):
        future_map = {}

        if self.pcoord_calculator is not None:
            for segment in segments:
                if segment.initpoint_type == segment.InitPoint.CONTINUES:
                    parent = self._prev_iter_segments[segment.parent_id]
                else:
                    parent = None

                future = self.work_manager.submit(self.pcoord_calculator, args=(segment, parent))
                future_map[future] = segment
        else:
            for segment in segments:
                segment.pcoord = default_pcoord(segment)

            self._segments[segment.seg_id] = segment
            self._data_manager.write_pcoords(self._n_iter, segments)

        return future_map

    def _run_we(self):
        self._invoke_callbacks('pre_we')

        # Initialize the weight transfer graph with self-loops.
        segments = [s.copy(wtg_parent_ids=[s.seg_id]) for s in self._segments]

        # Assign walkers to bins.
        bins = self.bin_mapper(segments)
        _report_bin_statistics(bins)

        if self.bin_target_counts.ndim == 0:
            target_counts = np.repeat(self.bin_target_counts, len(bins))
        else:
            if len(self.bin_target_counts) != len(bins):
                raise ValueError("length of 'bin_target_counts' must match the number of bins")
            target_counts = self.bin_target_counts.copy()

        # Recycle walkers.
        # When a walker is recycled, it is removed from its bin (withheld
        # from resampling), and the target count of the bin is decreased by one,
        # down to a floor of one (the count must remain positive in case there
        # are other walkers in the bin).
        for n, sink in enumerate(self.sinks):
            p_recycled = 0
            n_recycled = 0

            for idx, bin in enumerate(bins):
                if segments := {segment for segment in bin if segment in sink}:
                    bin -= segments
                    target_counts[idx] = max(1, target_counts[idx].item() - len(segments))

                    for segment in segments:
                        self._segments[segment.seg_id].endpoint_type = segment.EndPoint.RECYCLED

                    p_recycled += sum(map(operator.attrgetter('weight'), segments))
                    n_recycled += len(segments)

            if n_recycled > 0:
                label = sink.label if sink.label else n
                logger.info(f'Recycled {p_recycled} probability ({n_recycled} walkers) from sink {label!r}')

        # Resample the remaining walkers.
        resampled_bins = list(map(self.resampler, bins, target_counts))

        self._resampled_segments = list(itertools.chain(*resampled_bins))

    def _prepare_new_iteration(self):
        # recycled walkers
        for segment in self._segments:
            if segment.endpoint_type != Segment.EndPoint.RECYCLED:
                continue

            new_segment = Segment(
                n_iter=self._n_iter + 1,
                weight=segment.weight,
                wtg_parent_ids=[segment.seg_id],
                parent_id=-1,
                status=Segment.Status.UNSET,
            )
            self._next_iter_segments.append(new_segment)

        # continuing walkers
        for segment in self._resampled_segments:
            self._segments[segment.seg_id].endpoint_type = Segment.EndPoint.CONTINUES

            new_segment = Segment(
                n_iter=self._n_iter + 1,
                weight=segment.weight,
                wtg_parent_ids=segment.wtg_parent_ids,
                parent_id=segment.seg_id,
                initial_state=segment.final_state,
                status=Segment.Status.PREPARED,
            )
            self._next_iter_segments.append(new_segment)

        # pruned walkers
        for segment in self._segments:
            if segment.endpoint_type == Segment.EndPoint.UNSET:
                segment.endpoint_type = Segment.EndPoint.MERGED

        self._data_manager.prepare_iteration(self._n_iter + 1, self._next_iter_segments)

    def _next_iteration(self):
        self._prev_iter_segments.clear()

        self._run_we()
        self._prepare_new_iteration()
        self._finalize_iteration()

        self._data_manager.current_iteration += 1
        self._n_iter += 1

        self._prev_iter_segments = self._segments
        self._segments = self._next_iter_segments.copy()
        self._resampled_segments.clear()
        self._next_iter_segments.clear()

    def register_callback(self, hook, function, priority=0):
        pass

    def _invoke_callbacks(self, hook):
        pass
