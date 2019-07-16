import numpy as np
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

pcoord_len = 11
pcoord_dtype = np.float32


class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 1
        self.pcoord_len = pcoord_len
        self.pcoord_dtype = pcoord_dtype
        binbounds = ([0.0] + [2.8, 2.88, 3.0, 3.10, 3.29, 3.79, 3.94, 4.12, 4.39, 5.43] + [5.90+1.0*i for i in range(0,11)] + [30,float('inf')]) 
        self.bin_mapper = RectilinearBinMapper([binbounds])
        self.bin_target_counts = np.empty((self.bin_mapper.nbins,), np.int)
        self.bin_target_counts[...] = 24
