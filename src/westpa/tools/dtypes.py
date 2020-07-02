'''Numpy/HDF5 data types shared among several WESTPA tools'''

import numpy as np

# Pick up a few data types from the WEST core if possible
try:
    from westpa.core.data_manager import n_iter_dtype, seg_id_dtype, weight_dtype
except ImportError:
    n_iter_dtype = np.uint32
    seg_id_dtype = np.int64
    weight_dtype = np.float64

# A quantity averaged over iterations
iter_block_ci_dtype = np.dtype(
    [
        ('iter_start', n_iter_dtype),
        ('iter_stop', n_iter_dtype),
        ('expected', np.float64),
        ('ci_lbound', np.float64),
        ('ci_ubound', np.float64),
        ('sterr', np.float64),
        ('corr_len', n_iter_dtype),
    ]
)

# A quantity to store event duration distribution stuff.
# Comes from the old w_kinetics.

ed_list_dtype = np.dtype(
    [('istate', np.uint16), ('fstate', np.uint16), ('duration', np.float64), ('weight', np.float64), ('seg_id', seg_id_dtype)]
)
