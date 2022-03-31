import numpy as np
from westpa.core.binning import FuncBinMapper
from westpa.core.extloader import get_object


def map_binless(coords, mask, output, *args, **kwargs):
    '''Binning which adaptively places bins based on the positions of extrema segments and
    bottleneck segments, which are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.'''

    groups_per_dim = kwargs.get("groups_per_dim")
    group_function = get_object(kwargs.get("group_function"))
    ndim = len(groups_per_dim)

    if not np.any(mask):
        return output

    allcoords = np.copy(coords)
    allmask = np.copy(mask)

    isfinal = None
    splitting = False

    # the segments should be sent in by the driver as half initial segments and half final segments
    # allcoords contains all segments
    # coords should contain ONLY final segments
    if coords.shape[1] > ndim:
        if coords.shape[1] > ndim + 1:
            isfinal = allcoords[:, ndim + 1].astype(np.bool_)
        else:
            isfinal = np.ones(coords.shape[0], dtype=np.bool_)
        coords = coords[isfinal, :ndim]
        mask = mask[isfinal]
        splitting = True

    # in case where there is no final segments but initial ones in range
    if not np.any(mask):
        coords = allcoords[:, :ndim]
        mask = allmask
        splitting = False

    # filter the list of coordinates (which contains coordinates outside of the binless region)
    # to obtain only the ones we want to cluster
    # we do this per pcoord dimension
    for n in range(ndim):
        binless_coords = coords[mask, n]
        nsegs_binless = len(binless_coords)

        # we need to make sure that the number of segments in the binless region is greater than
        # the number of clusters we request
        # if only one segment in the binless region, assign it to a single cluster
        if nsegs_binless == 1:
            clusters = [0]
        # if there are more than one segment in the binless region but still the total is less than
        # our target number, adjust our target number to be the number of segments in the binless
        # region minus one
        elif nsegs_binless <= groups_per_dim[n]:
            clusters = group_function(binless_coords, nsegs_binless - 1, splitting)
        # if there are enough segments in the binless region, proceed as planned
        elif nsegs_binless > groups_per_dim[n]:
            clusters = group_function(binless_coords, groups_per_dim[0], splitting)

        # this is a good place to say this... output is a list which matches the length of allcoords
        # allcoords is a collection of all initial and final segment coords for that iteration
        # we first filtered those to only contain the final data points, since those are the ones we care
        # about clustering
        # we then filtered to only have the coords in the binless region, since, again, those are what we care about
        # we then assigned each to a cluster which is essentially a slitting index
        # all that's left is to find where each binless segment is in the output and insert the cluster index there
        for idx, val in enumerate(binless_coords):
            output[allcoords[:, n] == val] = clusters[idx]

        return output


class BinlessMapper(FuncBinMapper):
    '''Adaptively place bins in between minimum and maximum segments along
    the progress coordinte. Extrema and bottleneck segments are assigned
    to their own bins.'''

    def __init__(self, ngroups, group_function):
        kwargs = dict(groups_per_dim=ngroups, group_function=group_function)
        n_total_groups = np.prod(ngroups)
        super().__init__(map_binless, n_total_groups, kwargs=kwargs)
