from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy as np
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper, RecursiveBinMapper, FuncBinMapper
import westpa

from westpa_wexplore import wexplore, wex_utils
from scipy.spatial.distance import cdist

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

def eucl_dist(p, centers):
    d = np.array([np.linalg.norm(np.array(p - c,dtype=float)) for c in centers], dtype=np.float32)

    # the above is the Frobenius norm, for RMSD divide by N^1/2, where N is number of atoms
    # (Note that p has length 3*N)
    
    d /= math.sqrt(len(p)/3)
    
    return d

class System(WESTSystem):
    def initialize(self):
        # The number of dimensions should be the number of atoms that we have multipled by 3.
        self.pcoord_ndim = 3
        self.pcoord_len = 11
        self.pcoord_dtype = np.float32

        self.bin_mapper = wexplore.WExploreBinMapper(n_regions=[10,10,10], d_cut=[5, 2.0, 0.8], dfunc=eucl_dist)
        # The initial center is on the coordinates of one of the basis states.
        init_struct = np.loadtxt('18-crown-6-K+.pdb', dtype=str)
        atom_coords = init_struct[5:8]
        atom_coords = atom_coords.astype(float).flatten()
        self.bin_mapper.centers = [atom_coords]
        self.bin_mapper.add_bin(None, 0)
        self.max_replicas = 48
        self.bin_target_counts = self.bin_mapper.balance_replicas(self.max_replicas,
                                np.array([0,], np.int_))


def pcoord_loader(fieldname, pcoord_return_filename, destobj, single_point):
    """Read progress coordinate data into the ``pcoord`` field on ``destobj``. 
    An exception will be raised if the data is malformed.  If ``single_point`` is true,
    then only one (N-dimensional) point will be read, otherwise system.pcoord_len points
    will be read.
    """
    
    system = westpa.rc.get_system_driver()
    natoms = 1
    
    assert fieldname == 'pcoord'
    
    init_struct = np.loadtxt(pcoord_return_filename, dtype=str)
    # We're pulling in columns 5, 6, and 7 because this is where the X,Y,Z coords are in the pdb.
    try:
        atom_coords = init_struct[:,5:8]
    except:
        atom_coords = init_struct[5:8]
    pcoord = atom_coords.astype(float).flatten()
    
    if single_point:
        expected_shape = (system.pcoord_ndim,)
        if pcoord.ndim == 0:
            pcoord.shape = (1,)
    else:
        # We want to reshape the progress coordinate so that each row is a frame,
        # and each dimension is the number of atoms * 3.
        pcoord.shape = (11, natoms*3)
        expected_shape = (system.pcoord_len, system.pcoord_ndim)
        if pcoord.ndim == 1:
            pcoord.shape = (len(pcoord),1)
    if pcoord.shape != expected_shape:
        raise ValueError('progress coordinate data has incorrect shape {!r} [expected {!r}]'.format(pcoord.shape,
                                                                                                    expected_shape))
    destobj.pcoord = pcoord
    

