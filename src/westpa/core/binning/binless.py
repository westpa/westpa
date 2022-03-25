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

    #    weights = None
    isfinal = None
    #    splitting = False

    # the segments should be sent in by the driver as half initial segments and half final segments
    # allcoords contains all segments
    # coords should contain ONLY final segments
    if coords.shape[1] > ndim:
        if coords.shape[1] > ndim + 1:
            isfinal = allcoords[:, ndim + 1].astype(np.bool_)
        else:
            isfinal = np.ones(coords.shape[0], dtype=np.bool_)
        coords = coords[isfinal, :ndim]
        #        weights = allcoords[isfinal, ndim + 0]
        mask = mask[isfinal]
    #        splitting = True

    # in case where there is no final segments but initial ones in range
    if not np.any(mask):
        coords = allcoords[:, :ndim]
        mask = allmask
    #        weights = None
    #        splitting = False

    output = group_function(coords, groups_per_dim)

    return output


class BinlessMapper(FuncBinMapper):
    '''Adaptively place bins in between minimum and maximum segments along
    the progress coordinte. Extrema and bottleneck segments are assigned
    to their own bins.'''

    def __init__(self, ngroups, group_function):
        kwargs = dict(groups_per_dim=ngroups, group_function=group_function)
        n_total_groups = np.prod(ngroups)
        super().__init__(map_binless, n_total_groups, kwargs=kwargs)
