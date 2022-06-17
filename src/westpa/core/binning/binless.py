import numpy as np
from westpa.core.binning import FuncBinMapper
from westpa.core.extloader import get_object
import logging

log = logging.getLogger(__name__)


def map_binless(coords, mask, output, *args, **kwargs):
    '''Adaptively groups walkers according to a user-defined grouping function
    that is defined externally. Very general implementation but limited to
    only a two dimensional progress coordinate (for now).'''

    n_groups = kwargs.get("n_groups")
    n_dims = kwargs.get("n_dims")
    group_function = get_object(kwargs.get("group_function"))
    log.debug(f'binless arguments: {kwargs}')
    try:
        group_function_kwargs = kwargs.get('group_function_kwargs')['group_arguments']
    except KeyError:
        group_function_kwargs = {}
    ndim = n_dims

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
    # this is done with all dimensions at once
    binless_coords = coords[mask]
    nsegs_binless = len(binless_coords)

    # we need to make sure that the number of segments in the binless region is greater than
    # the number of clusters we request
    # if only one segment in the binless region, assign it to a single cluster
    if nsegs_binless == 1:
        clusters = [0]

    # if there are more than one segment in the binless region but still the total is less than
    # our target number, adjust our target number to be the number of segments in the binless
    # region minus one
    elif nsegs_binless < n_groups:
        clusters = group_function(binless_coords, nsegs_binless, splitting, **group_function_kwargs)
    # if there are enough segments in the binless region, proceed as planned
    elif nsegs_binless >= n_groups:
        clusters = group_function(binless_coords, n_groups, splitting, **group_function_kwargs)

    # this is a good place to say this... output is a list which matches the length of allcoords
    # allcoords is a collection of all initial and final segment coords for that iteration
    # we first filtered those to only contain the final data points, since those are the ones we care
    # about clustering
    # we then filtered to only have the coords in the binless region, since, again, those are what we care about
    # we then assigned each to a cluster which is essentially a slitting index
    # all that's left is to find where each binless segment is in the output and insert the cluster index there
    for idx, val in enumerate(binless_coords):
        if ndim > 1:
            mask2 = np.logical_and(allcoords[:, 0] == val[0], allcoords[:, 1] == val[1])
        else:
            mask2 = allcoords[:, 0] == val[0]
        output[mask2] = clusters[idx]

    return output


class BinlessMapper(FuncBinMapper):
    '''Adaptively group walkers according to a user-defined grouping
    function that is defined externally.'''

    def __init__(self, ngroups, ndims, group_function, **group_function_kwargs):
        kwargs = dict(n_groups=ngroups, n_dims=ndims, group_function=group_function, group_function_kwargs=group_function_kwargs)
        n_total_groups = ngroups
        super().__init__(map_binless, n_total_groups, kwargs=kwargs)
