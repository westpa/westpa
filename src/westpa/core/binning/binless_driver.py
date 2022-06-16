import logging

import numpy as np
from westpa.core.we_driver import WEDriver

log = logging.getLogger(__name__)


class BinlessDriver(WEDriver):
    def assign(self, segments, initializing=False):
        '''Assign segments to initial and final bins, and update the (internal) lists of used and available
        initial states. This function is adapted to the MAB scheme, so that the inital and final segments are
        sent to the bin mapper at the same time, otherwise the inital and final bin boundaries can be inconsistent.'''

        log.debug("BinlessDriver in use.")
        # collect initial and final coordinates into one place
        n_segments = len(segments)
        all_pcoords = np.empty((n_segments * 2, self.system.pcoord_ndim + 2), dtype=self.system.pcoord_dtype)

        for iseg, segment in enumerate(segments):
            all_pcoords[iseg] = np.append(segment.pcoord[0, :], [segment.weight, 0.0])
            all_pcoords[n_segments + iseg] = np.append(segment.pcoord[-1, :], [segment.weight, 1.0])

        # assign based on initial and final progress coordinates
        assignments = self.bin_mapper.assign(all_pcoords)
        initial_assignments = assignments[:n_segments]
        if initializing:
            final_assignments = initial_assignments
        else:
            final_assignments = assignments[n_segments:]

        initial_binning = self.initial_binning
        final_binning = self.final_binning
        flux_matrix = self.flux_matrix
        transition_matrix = self.transition_matrix
        for (segment, iidx, fidx) in zip(segments, initial_assignments, final_assignments):
            initial_binning[iidx].add(segment)
            final_binning[fidx].add(segment)
            flux_matrix[iidx, fidx] += segment.weight
            transition_matrix[iidx, fidx] += 1

        n_recycled_total = self.n_recycled_segs
        n_new_states = n_recycled_total - len(self.avail_initial_states)

        log.debug(
            '{} walkers scheduled for recycling, {} initial states available'.format(
                n_recycled_total, len(self.avail_initial_states)
            )
        )

        if n_new_states > 0:
            return n_new_states
        else:
            return 0
