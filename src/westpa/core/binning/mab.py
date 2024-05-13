import logging
import numpy as np
import westpa
from westpa.core.binning import FuncBinMapper
from os.path import expandvars


log = logging.getLogger(__name__)


class MABBinMapper(FuncBinMapper):
    """
    Adaptively place bins in between minimum and maximum segments along
    the progress coordinte. Extrema and bottleneck segments are assigned
    to their own bins.

    """

    def __init__(
        self,
        nbins,
        direction=None,
        skip=None,
        bottleneck=True,
        pca=False,
        mab_log=False,
        bin_log=False,
        bin_log_path="$WEST_SIM_ROOT/binbounds.log",
    ):
        """
        Parameters
        ----------
        nbins : list of int
            List of int for nbins in each dimension.
        direction : Union(list of int, None), default: None
            List of int for 'direction' in each dimension.
            Direction options are as follows:
                0   : default split at leading and lagging boundaries
                1   : split at leading boundary only
                -1  : split at lagging boundary only
                86  : no splitting at either leading or lagging boundary
        skip : Union(list of int, None), default: None
            List of int for each dimension. Default None for skip=0.
            Set to 1 to 'skip' running mab in a dimension.
        bottleneck : bool, default: True
            Whether to turn on or off bottleneck walker splitting.
        pca : bool, default: False
            Can be True or False (default) to run PCA on pcoords before bin assignment.
        mab_log : bool, default: False
            Whether to output mab info to west.log.
        bin_log : bool, default: False
            Whether to output mab bin boundaries to bin_log_path file.
        bin_log_path : str, default: "$WEST_SIM_ROOT/binbounds.log"
            Path to output bin boundaries.

        """
        # Verifying parameters
        if nbins is None:
            raise ValueError("nbins_per_dim is missing")
        ndim = len(nbins)

        if direction is None:
            direction = [0] * ndim
        elif len(direction) != ndim:
            direction = [0] * ndim
            log.warning("Direction list is not the correct dimensions, setting to defaults.")

        if skip is None:
            skip = [0] * ndim
        elif len(skip) != ndim:
            skip = [0] * ndim
            log.warning("Skip list is not the correct dimensions, setting to defaults.")

        kwargs = dict(
            nbins_per_dim=nbins,
            direction=direction,
            skip=skip,
            bottleneck=bottleneck,
            pca=pca,
            mab_log=mab_log,
            bin_log=bin_log,
            bin_log_path=bin_log_path,
        )

        n_total_bins = self.determine_total_bins(**kwargs)

        super().__init__(map_mab, n_total_bins, kwargs=kwargs)

    def determine_total_bins(self, nbins_per_dim, direction, skip, bottleneck, **kwargs):
        """
        The following is neccessary because functional bin mappers need to "reserve"
        bins and tell the sim manager how many bins they will need to use, this is
        determined by taking all direction/skipping info into account.

        Parameters
        ----------
        nbins_per_dim : int
            Number of total bins in each direction.
        direction : list of int
            Direction in each dimension. See __init__ for more information.
        skip : list of int
            List of 0s and 1s indicating whether to skip each dimension.
        bottleneck : bool
            Whether to include separate bin for bottleneck walker(s).
        **kwargs : dict
            Arbitary keyword arguments. Contains unneeded MAB parameters.

        Returns
        -------
        n_total_bins : int
            Number of total bins.

        """
        n_total_bins = np.prod(nbins_per_dim)
        ndim = len(nbins_per_dim)
        for i in range(ndim):
            if skip[i] == 0:
                if direction[i] != 0:
                    n_total_bins += 1 + 1 * bottleneck
                else:
                    n_total_bins += 2 + 2 * bottleneck
            else:
                n_total_bins -= nbins_per_dim[i] - 1
                n_total_bins += 1 * ndim  # or else it will be one bin short
        return n_total_bins


def map_mab(coords, mask, output, *args, **kwargs):
    """
    Binning which adaptively places bins based on the positions of extrema segments and
    bottleneck segments, which are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.

    Parameters
    ----------
    coords : ndarray
        An array with pcoord and weight info.
    mask : ndarray
        Array of 1 (True) and 0 (False), to filter out unwanted segment info.
    output : list
        The main list that, for each segment, holds the bin assignment.
    *args : list
        Variable length arguments.
    **kwargs : dict
        Arbitary keyword arguments. Contains most of the MAB-needed parameters.

    Returns
    ------
    output : list
        The main list that, for each segment, holds the bin assignment.

    """

    # Argument Processing
    nbins_per_dim = kwargs.get("nbins_per_dim")
    ndim = len(nbins_per_dim)
    pca = kwargs.get("pca", False)
    bottleneck = kwargs.get("bottleneck", True)
    direction = kwargs.get("direction", ([0] * ndim))
    skip = kwargs.get("skip", ([0] * ndim))
    mab_log = kwargs.get("mab_log", False)
    bin_log = kwargs.get("bin_log", False)
    bin_log_path = kwargs.get("bin_log_path", "$WEST_SIM_ROOT/binbounds.log")

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
    n_bottleneck_filled = 0

    for i in range(0, ndim):
        # for single direction, 1 boundary walker
        if direction[i] == 1 or direction[i] == -1:
            bottleneck_base += 1
        # 2 boundary walkers with 0 direction
        elif direction[i] == 0:
            bottleneck_base += 2
        # for 86 direction, no boundary walkers so offset of 0
        elif direction[i] == 86:
            bottleneck_base += 0

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
                    if direction[n] == -1:
                        if coord == flipdifflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            n_bottleneck_filled += 1
                            break

                    if direction[n] == 1:
                        if coord == difflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            n_bottleneck_filled += 1
                            break

                    # both directions when using 0 or with
                    # special value of 86 for no lead/lag split
                    if direction[n] == 0 or direction[n] == 86:
                        if coord == difflist[n]:
                            holder = bottleneck_base + n
                            special = True
                            n_bottleneck_filled += 1
                            break
                        elif coord == flipdifflist[n]:
                            holder = bottleneck_base + n + 1
                            special = True
                            n_bottleneck_filled += 1
                            break

                # assign boundary walkers, taking directionality into account
                if direction[n] == -1:
                    if coord == minlist[n]:
                        holder = boundary_base + n
                        special = True
                        break

                elif direction[n] == 1:
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

                # special value for direction with no lead/lag split
                elif direction[n] == 86:
                    # westpa.rc.pstatus(f"No lead/lag split for dim {n}")
                    # westpa.rc.pflush()
                    # nornmally adds to special bin but here just leaving it forever empty
                    # holder = boundary_base + n
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

    if bin_log and report:
        if westpa.rc.sim_manager.n_iter:
            with open(expandvars(bin_log_path), 'a') as bb_file:
                # Iteration Number
                bb_file.write(f'iteration: {westpa.rc.sim_manager.n_iter}\n')
                bb_file.write('bin boundaries: ')
                for n in range(ndim):
                    # Write binbounds per dim
                    bb_file.write(f'{np.linspace(minlist[n], maxlist[n], nbins_per_dim[n] + 1)}\t')
                    # Min/Max pcoord
                bb_file.write(f'\nmin/max pcoord: {minlist} {maxlist}\n')
                bb_file.write(f'bottleneck bins: {n_bottleneck_filled}\n')
                if n_bottleneck_filled > 0:
                    # Bottlenecks bins exist (passes any of the if bottleneck: checks)
                    bb_file.write(f'bottleneck pcoord: {flipdifflist} {difflist}\n\n')
                else:
                    bb_file.write('\n')

    return output
