import numpy as np
from westpa.core.binning import FuncBinMapper


def map_mab(coords, mask, output, *args, **kwargs):
    '''Binning which adaptively places bins based on the positions of extrema segments and
    bottleneck segments, which are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.'''

    pca = kwargs.pop("pca", False)
    bottleneck = kwargs.pop("bottleneck", True)
    nbins_per_dim = kwargs.get("nbins_per_dim")
    ndim = len(nbins_per_dim)

    if not np.any(mask):
        return output

    allcoords = np.copy(coords)
    allmask = np.copy(mask)

    weights = None
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
        weights = allcoords[isfinal, ndim + 0]
        mask = mask[isfinal]
        splitting = True

    # in case where there is no final segments but initial ones in range
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
        covcoords = np.cov(np.transpose(varcoords))
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
                    temp[p][1] = 10 ** -39
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

    # assign segments to bins
    # the total number of linear bins + 2 boundary bins each dim
    boundary_base = np.prod(nbins_per_dim)
    bottleneck_base = boundary_base + 2 * ndim
    for i in range(len(output)):
        if not allmask[i]:
            continue

        special = False
        holder = 0
        if splitting:
            for n in range(ndim):
                coord = allcoords[i][n]

                if bottleneck:
                    if coord == difflist[n]:
                        holder = bottleneck_base + 2 * n
                        special = True
                        break
                    elif coord == flipdifflist[n]:
                        holder = bottleneck_base + 2 * n + 1
                        special = True
                        break
                if coord == minlist[n]:
                    holder = boundary_base + 2 * n
                    special = True
                    break
                elif coord == maxlist[n]:
                    holder = boundary_base + 2 * n + 1
                    special = True
                    break

        if not special:
            for n in range(ndim):
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
                    raise ValueError("Walker out of boundary")

                holder += bin_number * np.prod(nbins_per_dim[:n])

        output[i] = holder

    return output


class MABBinMapper(FuncBinMapper):
    '''Adaptively place bins in between minimum and maximum segments along
    the progress coordinte. Extrema and bottleneck segments are assigned
    to their own bins.'''

    def __init__(self, nbins, bottleneck=True, pca=False):
        kwargs = dict(nbins_per_dim=nbins, bottleneck=bottleneck, pca=pca)
        ndim = len(nbins)
        n_total_bins = np.prod(nbins) + ndim * (2 + 2 * bottleneck)
        super().__init__(map_mab, n_total_bins, kwargs=kwargs)
