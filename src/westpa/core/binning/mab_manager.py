import logging

from westpa.core.sim_manager import WESimManager, grouper
from westpa.core.states import InitialState, pare_basis_initial_states
from westpa.core import wm_ops

log = logging.getLogger(__name__)


class MABSimManager(WESimManager):
    def report_bin_statistics(self, bins, save_summary=False):
        self.rc.pstatus("MAB binning in use.")
        super().report_bin_statistics(bins, save_summary)

    def propagate(self):
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
