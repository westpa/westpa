import logging

from westpa.core.binning.mab import MABBinMapper
from westpa.core.sim_manager import WESimManager, grouper
from westpa.core.states import InitialState, pare_basis_initial_states
from westpa.core import wm_ops
from westpa.core.segment import Segment
import numpy as np

log = logging.getLogger(__name__)


class MABSimManager(WESimManager):
    def initialize_simulation(self, basis_states, target_states, start_states, segs_per_state=1, suppress_we=False):
        if len(target_states) > 0:
            if isinstance(self.system.bin_mapper, MABBinMapper):
                log.error("MABBinMapper cannot be an outer binning scheme with a target state\n")

        super().initialize_simulation(
            basis_states, target_states, start_states, segs_per_state=segs_per_state, suppress_we=suppress_we
        )

    def propagate(self):
        log.debug("MABSimManager in use")
        segments = list(self.incomplete_segments.values())
        log.debug('iteration {:d}: propagating {:d} segments'.format(self.n_iter, len(segments)))

        # all futures dispatched for this iteration
        futures = set()
        segment_futures = set()

        # Immediately dispatch any necessary initial state generation
        istate_gen_futures = self.get_istate_futures()
        futures.update(istate_gen_futures)

        # Dispatch propagation tasks using work manager
        for segment_block in grouper(self.propagator_block_size, segments):
            segment_block = [_f for _f in segment_block if _f]
            pbstates, pistates = pare_basis_initial_states(
                self.current_iter_bstates, list(self.current_iter_istates.values()), segment_block
            )
            future = self.work_manager.submit(wm_ops.propagate, args=(pbstates, pistates, segment_block))
            futures.add(future)
            segment_futures.add(future)

        while futures:
            # TODO: add capacity for timeout or SIGINT here
            future = self.work_manager.wait_any(futures)
            futures.remove(future)

            if future in segment_futures:
                segment_futures.remove(future)
                incoming = future.get_result()
                self.n_propagated += 1

                self.segments.update({segment.seg_id: segment for segment in incoming})
                self.completed_segments.update({segment.seg_id: segment for segment in incoming})

                new_istate_futures = self.get_istate_futures()
                istate_gen_futures.update(new_istate_futures)
                futures.update(new_istate_futures)

                with self.data_manager.expiring_flushing_lock():
                    self.data_manager.update_segments(self.n_iter, incoming)

            elif future in istate_gen_futures:
                istate_gen_futures.remove(future)
                _basis_state, initial_state = future.get_result()
                log.debug('received newly-prepared initial state {!r}'.format(initial_state))
                initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                with self.data_manager.expiring_flushing_lock():
                    self.data_manager.update_initial_states([initial_state], n_iter=self.n_iter + 1)
                self.we_driver.avail_initial_states[initial_state.state_id] = initial_state
            else:
                log.error('unknown future {!r} received from work manager'.format(future))
                raise AssertionError('untracked future {!r}'.format(future))

        self.we_driver.assign(self.segments.values())
        self.get_istate_futures()
        log.debug('done with propagation')
        self.save_bin_data()
        self.data_manager.flush_backing()

    def prepare_iteration(self):
        log.debug('beginning iteration {:d}'.format(self.n_iter))

        # the WE driver needs a list of all target states for this iteration
        # along with information about any new weights introduced (e.g. by recycling)
        target_states = self.data_manager.get_target_states(self.n_iter)
        new_weights = self.data_manager.get_new_weight_data(self.n_iter)

        self.we_driver.new_iteration(target_states=target_states, new_weights=new_weights)

        # Get basis states used in this iteration
        self.current_iter_bstates = self.data_manager.get_basis_states(self.n_iter)

        # Get the segments for this iteration and separate into complete and incomplete
        if self.segments is None:
            segments = self.segments = {segment.seg_id: segment for segment in self.data_manager.get_segments()}
            log.debug('loaded {:d} segments'.format(len(segments)))
        else:
            segments = self.segments
            log.debug('using {:d} pre-existing segments'.format(len(segments)))

        completed_segments = self.completed_segments = {}
        incomplete_segments = self.incomplete_segments = {}
        for segment in segments.values():
            if segment.status == Segment.SEG_STATUS_COMPLETE:
                completed_segments[segment.seg_id] = segment
            else:
                incomplete_segments[segment.seg_id] = segment
        log.debug('{:d} segments are complete; {:d} are incomplete'.format(len(completed_segments), len(incomplete_segments)))

        if len(incomplete_segments) == len(segments):
            # Starting a new iteration
            self.rc.pstatus('Beginning iteration {:d}'.format(self.n_iter))
        elif incomplete_segments:
            self.rc.pstatus('Continuing iteration {:d}'.format(self.n_iter))
        self.rc.pstatus(
            '{:d} segments remain in iteration {:d} ({:d} total)'.format(len(incomplete_segments), self.n_iter, len(segments))
        )

        # Get the initial states active for this iteration (so that the propagator has them if necessary)
        self.current_iter_istates = {
            state.state_id: state for state in self.data_manager.get_segment_initial_states(list(segments.values()))
        }
        log.debug('This iteration uses {:d} initial states'.format(len(self.current_iter_istates)))

        n_segments = len(segments)
        pcoords_with_weights = np.empty((n_segments, self.system.pcoord_ndim + 2), dtype=self.system.pcoord_dtype)

        for iseg, segment in enumerate(segments.values()):
            pcoords_with_weights[iseg] = np.append(segment.pcoord[0, :], [segment.weight, 1.0])

        # Assign this iteration's segments' initial points to bins and report on bin population
        initial_binning = self.system.bin_mapper.construct_bins()
        initial_assignments = self.system.bin_mapper.assign(pcoords_with_weights)
        for (segment, assignment) in zip(iter(segments.values()), initial_assignments):
            initial_binning[assignment].add(segment)
        self.report_bin_statistics(initial_binning, [], save_summary=True)
        del pcoords_with_weights, initial_binning

        self.rc.pstatus("MAB binning in use")

        self.rc.pstatus("Bottleneck bin occupancy may not be accurately reported")

        self.rc.pstatus("Waiting for segments to complete...")

        # Let the WE driver assign completed segments
        if completed_segments and len(incomplete_segments) == 0:
            self.we_driver.assign(list(completed_segments.values()))

        # load restart data
        self.data_manager.prepare_segment_restarts(
            incomplete_segments.values(), self.current_iter_bstates, self.current_iter_istates
        )

        # Get the basis states and initial states for the next iteration, necessary for doing on-the-fly recycling
        self.next_iter_bstates = self.data_manager.get_basis_states(self.n_iter + 1)
        self.next_iter_bstate_cprobs = np.add.accumulate([bstate.probability for bstate in self.next_iter_bstates])

        self.we_driver.avail_initial_states = {
            istate.state_id: istate for istate in self.data_manager.get_unused_initial_states(n_iter=self.n_iter + 1)
        }
        log.debug('{:d} unused initial states found'.format(len(self.we_driver.avail_initial_states)))

        # Invoke callbacks
        self.invoke_callbacks(self.prepare_iteration)

        # dispatch and immediately wait on result for prep_iter
        log.debug('dispatching propagator prep_iter to work manager')
        self.work_manager.submit(wm_ops.prep_iter, args=(self.n_iter, segments)).get_result()
