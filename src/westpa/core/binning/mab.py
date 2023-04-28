import numpy as np
from westpa.core.binning import FuncBinMapper
import logging
import westpa

log = logging.getLogger(__name__)


def map_mab(coords, mask, output, *args, **kwargs):
    '''Binning which adaptively places bins based on the positions of extrema segments and
    bottleneck segments, which are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.'''

    pca = kwargs.get("pca")
    bottleneck = kwargs.get("bottleneck")
    direction = kwargs.get("direction")
    skip = kwargs.get("skip")
    nbins_per_dim = kwargs.get("nbins_per_dim")
    mab_log = kwargs.get("mab_log")
    ndim = len(nbins_per_dim)

    if not np.any(mask):
        return output

    if skip is None:
        skip = [0] * ndim

    allcoords = np.copy(coords)
    allmask = np.copy(mask)

    weights = None
    isfinal = None
    splitting = False
    report = False

    # the segments should be sent in by the driver as half initial segments and half final segments
    # allcoords contains all segments
    # coords should contain ONLY final segments
    if coords.shape[1] > ndim:
        if coords[0, -1] == 0:
            report = True
        if coords.shape[1] > ndim + 1:
            isfinal = allcoords[:, ndim + 1].astype(np.bool_)
        else:
            isfinal = np.ones(coords.shape[0], dtype=np.bool_)
        coords = coords[isfinal, :ndim]
        weights = allcoords[isfinal, ndim + 0]
        mask = mask[isfinal]
        splitting = True

    if not np.any(mask):
        coords = allcoords[:, :ndim]
        mask = allmask
        weights = None
        splitting = False

    varcoords = np.copy(coords)
    originalcoords = np.copy(coords)
    if pca and len(output) > 1:
        colavg = np.mean(coords, axis=0)
        for i in range(len(coords)):
            for j in range(len(coords[i])):
                varcoords[i][j] = coords[i][j] - colavg[j]
        covcoords = np.cov(np.transpose(varcoords), aweights=weights)
        eigval, eigvec = np.linalg.eigh(covcoords)
        eigvec = eigvec[:, np.argmax(np.absolute(eigvec), axis=1)]
        for i in range(len(eigvec)):
            if eigvec[i, i] < 0:
                eigvec[:, i] = -1 * eigvec[:, i]
        for i in range(ndim):
            for j in range(len(output)):
                coords[j][i] = np.dot(varcoords[j], eigvec[:, i])

    maxlist = []
    minlist = []
    difflist = []
    flipdifflist = []
    for n in range(ndim):
        # identify the boundary segments
        maxcoord = np.max(coords[mask, n])
        mincoord = np.min(coords[mask, n])
        maxlist.append(maxcoord)
        minlist.append(mincoord)

        # detect the bottleneck segments, this uses the weights
        if splitting:
            temp = np.column_stack((originalcoords[mask, n], weights[mask]))
            sorted_indices = temp[:, 0].argsort()
            temp = temp[sorted_indices]
            for p in range(len(temp)):
                if temp[p][1] == 0:
                    temp[p][1] = 10**-323
            fliptemp = np.flipud(temp)

            difflist.append(None)
            flipdifflist.append(None)
            maxdiff = 0
            flipmaxdiff = 0
            for i in range(1, len(temp) - 1):
                comprob = 0
                flipcomprob = 0
                j = i + 1
                while j < len(temp):
                    comprob = comprob + temp[j][1]
                    flipcomprob = flipcomprob + fliptemp[j][1]
                    j = j + 1
                diff = -np.log(comprob) + np.log(temp[i][1])
                if diff > maxdiff:
                    difflist[n] = temp[i][0]
                    maxdiff = diff
                flipdiff = -np.log(flipcomprob) + np.log(fliptemp[i][1])
                if flipdiff > flipmaxdiff:
                    flipdifflist[n] = fliptemp[i][0]
                    flipmaxdiff = flipdiff

    if mab_log and report:
        westpa.rc.pstatus("################ MAB stats ################")
        westpa.rc.pstatus("minima in each dimension:      {}".format(minlist))
        westpa.rc.pstatus("maxima in each dimension:      {}".format(maxlist))
        westpa.rc.pstatus("direction in each dimension:   {}".format(direction))
        westpa.rc.pstatus("skip in each dimension:        {}".format(skip))
        westpa.rc.pstatus("###########################################")
        westpa.rc.pflush()

    # assign segments to bins
    # the total number of linear bins is the boundary base
    boundary_base = np.prod(nbins_per_dim)

    # the bottleneck base is offset by the number of boundary walkers,
    # which is two per dimension unless there is a direction specified
    # in a particluar dimension, then it's just one
    bottleneck_base = boundary_base

    for i in range(0, ndim):
        if direction[i] != 0:
            bottleneck_base += 1
        else:
            bottleneck_base += 2

    # if a dimension is being "skipped", leave only one bin total as
    # the offset
    for i in range(0, ndim):
        if skip[i] != 0:
            boundary_base -= nbins_per_dim[i] - 1

    for i in range(len(output)):
        if not allmask[i]:
            continue

        # special means either a boundary or bottleneck walker (not a walker in the linear space)
        special = False
        # this holder is the bin number, which only needs to be unique for different walker groups
        holder = 0
        if splitting:
            for n in range(ndim):
                coord = allcoords[i][n]

                # if skipped, just assign the walkers to the same bin (offset of boundary base)
                if skip[n] != 0:
                    holder = boundary_base + n
                    break

                # assign bottlenecks, taking directionality into account
                if bottleneck:
                    if direction[n] < 0:
                        if coord == flipdifflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            break

                    if direction[n] > 0:
                        if coord == difflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            break

                    if direction[n] == 0:
                        if coord == difflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            break
                        elif coord == flipdifflist[n]:
                            holder = bottleneck_base + n + 1
                            special = True
                            break

                # assign boundary walkers, taking directionality into account
                if direction[n] < 0:
                    if coord == minlist[n]:
                        holder = boundary_base + n
                        special = True
                        break

                elif direction[n] > 0:
                    if coord == maxlist[n]:
                        holder = boundary_base + n
                        special = True
                        break

                elif direction[n] == 0:
                    if coord == minlist[n]:
                        holder = boundary_base + n
                        special = True
                        break
                    elif coord == maxlist[n]:
                        holder = boundary_base + n + 1
                        special = True
                        break

        # the following are for the "linear" portion
        if not special:
            for n in range(ndim):

                # if skipped, it's added to the same bin as the special walkers above
                if skip[n] != 0:
                    holder = boundary_base + n
                    break

                coord = allcoords[i][n]
                nbins = nbins_per_dim[n]
                minp = minlist[n]
                maxp = maxlist[n]

                bins = np.linspace(minp, maxp, nbins + 1)
                bin_number = np.digitize(coord, bins) - 1

                if isfinal is None or not isfinal[i]:
                    if bin_number >= nbins:
                        bin_number = nbins - 1
                    elif bin_number < 0:
                        bin_number = 0
                elif bin_number >= nbins or bin_number < 0:
                    if np.isclose(bins[-1], coord):
                        bin_number = nbins - 1
                    elif np.isclose(bins[0], coord):
                        bin_number = 0
                    else:
                        raise ValueError("Walker out of boundary")

                holder += bin_number * np.prod(nbins_per_dim[:n])

        # output is the main list that, for each segment, holds the bin assignment
        output[i] = holder

    return output


class MABBinMapper(FuncBinMapper):
    '''Adaptively place bins in between minimum and maximum segments along
    the progress coordinte. Extrema and bottleneck segments are assigned
    to their own bins.'''

    def __init__(self, nbins, direction=None, skip=None, bottleneck=True, pca=False, mab_log=False):
        # Verifying parameters
        if nbins is None:
            raise ValueError("nbins_per_dim is missing")
        ndim = len(nbins)

        if direction is None:
            direction = [0] * ndim
        elif len(direction) != ndim:
            direction = [0] * ndim
            log.warn("Direction list is not the correct dimensions, setting to defaults.")

        if skip is None:
            skip = [0] * ndim
        elif len(skip) != ndim:
            skip = [0] * ndim
            log.warn("Skip list is not the correct dimensions, setting to defaults.")

        kwargs = dict(nbins_per_dim=nbins, direction=direction, skip=skip, bottleneck=bottleneck, pca=pca, mab_log=mab_log)
        # the following is neccessary because functional bin mappers need to "reserve"
        # bins and tell the sim manager how many bins they will need to use, this is
        # determined by taking all direction/skipping info into account
        n_total_bins = np.prod(nbins)
        for i in range(0, ndim):
            if skip[i] == 0:
                if direction[i] != 0:
                    n_total_bins += 1 + 1 * bottleneck
                else:
                    n_total_bins += 2 + 2 * bottleneck
            else:
                n_total_bins -= nbins[i] - 1
                n_total_bins += 1 * ndim  # or else it will be one bin short
        super().__init__(map_mab, n_total_bins, kwargs=kwargs)
