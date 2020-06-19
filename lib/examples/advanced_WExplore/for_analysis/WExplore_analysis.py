import numpy as np
import h5py

import os
import sys

# adding some westpa things to the path

WEST_ROOT = os.environ['WEST_ROOT']
for lib in ['lib/wwmgr', 'src', 'lib/west_tools']:
    path = os.path.join(WEST_ROOT, lib)
    if path not in sys.path:
        sys.path.append(path)

import west
        
# the data manager object helps access simulation data written to the west.h5 file
data_manager = west.rc.get_data_manager()
data_manager.open_backing(mode='r')

string_hashes = data_manager.we_h5file['bin_topologies']['index'][:]['hash']
all_mappers = []
for si, shash in enumerate(string_hashes):
    bin_mapper = data_manager.get_bin_mapper(shash)
    all_mappers.append(bin_mapper)

# each mapper in all_mappers holds the coordinates of all the images

# the code below fetches the coordinates of all the regions from the 
# latest mapper

latest_mapper = all_mappers[-1]
level_inds = latest_mapper.level_indices

# get region indices
biggest_regions_inds = level_inds[0]
medium_regions_inds = level_inds[1]
small_regions_inds = level_inds[2]

# find region centers
big_regions = latest_mapper.fetch_centers(biggest_regions_inds)
med_regions = latest_mapper.fetch_centers(medium_regions_inds)
sml_regions = latest_mapper.fetch_centers(small_regions_inds)

print("There are {0} big regions, {1} medium regions, and {2} small regions"
      .format(len(big_regions),len(med_regions),len(sml_regions)))

# Since these coordinates are simply the pcoord values used in system.py,
# we can import our distance function and apply it to analyze the images

from system import eucl_dist

n_med = len(med_regions)
dists = np.zeros((n_med,n_med))
for i, x in enumerate(med_regions):
    for j, y in enumerate(med_regions):
        dists[i,j] = eucl_dist(x,[y])

print("The minimum non-zero distance between images is: {0}".format(dists[np.nonzero(dists)].min()))

# We can also write these to trajectories using mdtraj

import mdtraj as mdj

ion_pdb = mdj.load_pdb('18-crown-6-K+.pdb')

big_reg_reshape = big_regions.reshape((len(big_regions),ion_pdb.n_atoms,3))/10
med_reg_reshape = med_regions.reshape((len(med_regions),ion_pdb.n_atoms,3))/10
sml_reg_reshape = sml_regions.reshape((len(sml_regions),ion_pdb.n_atoms,3))/10

mdj.Trajectory(big_reg_reshape,ion_pdb.top).save_dcd('big_regions.dcd')
mdj.Trajectory(med_reg_reshape,ion_pdb.top).save_dcd('med_regions.dcd')
mdj.Trajectory(sml_reg_reshape,ion_pdb.top).save_dcd('sml_regions.dcd')

# Great!  We hope you had fun analyzing your WExplore results!
