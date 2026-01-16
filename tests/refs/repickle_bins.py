#!/usr/bin/env python
#
# repickle_bins.py
#
#
import h5py

from westpa.tools.binning import mapper_from_hdf5


with h5py.File('west_ref_new.h5', 'a') as h5file:
    i = 0
    for iiter in range(1, len(h5file['summary'])):
        try:
            binhash = h5file[f'iterations/iter_{iiter:>08}'].attrs['binhash']
            (a, b, c) = mapper_from_hdf5(h5file['bin_topologies'], binhash)
            new_pickle, new_hash = a.pickle_and_hash()

            index_row = h5file['bin_topologies/index'][i]
            index_row['hash'] = new_hash
            index_row['pickle_len'] = len(new_pickle)
            h5file['bin_topologies/index'][i] = index_row
            h5file['bin_topologies/pickles'][i, : len(new_pickle)] = memoryview(bytes(new_pickle))
            h5file[f'iterations/iter_{iiter:>08}'].attrs['binhash'] = new_hash
            i += 1
        except KeyError:
            pass
