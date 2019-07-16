import numpy as np
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

pcoord_len = 2
pcoord_dtype = np.float32


class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_len = pcoord_len
        self.pcoord_dtype = pcoord_dtype

        binbounds = [0.0] + np.arange(3.0, 8.1, 0.1).tolist() + [float('inf')]

        self.bin_mapper = RectilinearBinMapper([binbounds])
        self.bin_target_counts = np.empty((self.bin_mapper.nbins,), np.int)
        self.bin_target_counts[...] = 12


def gen_state_labels(mapper):
    dtest = np.linspace(2.5, 8.7, 1000)[:,None]

    state_list = []
    assignments = mapper.assign(dtest)
    uassign, indx = np.unique(assignments, return_index=True)
    pcoords = dtest[indx]

    state_a = np.where(pcoords < 4.2)[0]
    state_b = np.where(pcoords > 7.0)[0]

    state_list.append({'label': 'state_a', 'coords': pcoords[state_a]})
    state_list.append({'label': 'state_b', 'coords': pcoords[state_b]})

    return state_list
