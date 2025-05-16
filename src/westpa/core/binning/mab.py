import logging
from typing import List, Optional
import numpy as np
import westpa
from westpa.core.binning import FuncBinMapper
from os.path import expandvars

log = logging.getLogger(__name__)


class MABBinMapper(FuncBinMapper):
    """
    Adaptively place bins between minimum and maximum segments along
    the progress coordinate. Extrema and bottleneck segments are assigned
    to their own bins.
    """

    def __init__(
        self,
        nbins: List[int],
        direction: Optional[List[int]] = None,
        skip: Optional[List[int]] = None,
        bottleneck: bool = True,
        pca: bool = False,
        mab_log: bool = False,
        bin_log: bool = False,
        bin_log_path: str = "$WEST_SIM_ROOT/binbounds.log",
    ):
        """
        Parameters
        ----------
        nbins : list of int
            List of number of bins in each dimension.
        direction : Optional[list of int], default: None
            List of directions in each dimension. Direction options:
                0   : default split at leading and lagging boundaries
                1   : split at leading boundary only
                -1  : split at lagging boundary only
                86  : no splitting at either leading or lagging boundary (both bottlenecks included)
        skip : Optional[list of int], default: None
            List of skip flags for each dimension. Default None (no skipping).
        bottleneck : bool, default: True
            Whether to enable bottleneck walker splitting.
        pca : bool, default: False
            Whether to perform PCA on progress coordinates before bin assignment.
        mab_log : bool, default: False
            Whether to output MAB info to west.log.
        bin_log : bool, default: False
            Whether to output MAB bin boundaries to a log file.
        bin_log_path : str, default: "$WEST_SIM_ROOT/binbounds.log"
            Path to output bin boundaries.
        """
        # Verifying parameters
        if nbins is None:
            raise ValueError("nbins is missing")
        ndim = len(nbins)

        direction = direction or [0] * ndim
        if len(direction) != ndim:
            direction = [0] * ndim
            log.warning("Direction list is not the correct dimensions, setting to defaults.")

        skip = skip or [0] * ndim
        if len(skip) != ndim:
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

    def determine_total_bins(
        self, nbins_per_dim: List[int], direction: List[int], skip: List[int], bottleneck: bool, **kwargs
    ) -> int:
        """
        Calculate the total number of bins needed, taking direction and skipping into account.
        This function is necessary because functional bin mappers need to "reserve"
        bins and tell the sim manager how many bins they will need to use, this is
        determined by taking all direction/skipping info into account.

        Parameters
        ----------
        nbins_per_dim : list of int
            Number of total bins in each dimension within the linear portion.
        direction : list of int
            Direction in each dimension.
        skip : list of int
            List indicating whether to skip each dimension.
        bottleneck : bool
            Whether to include a separate bin for bottleneck walker(s).
        **kwargs : dict
            Additional MAB parameters (unused).

        Returns
        -------
        n_total_bins : int
            Number of total bins.
        """
        # Update nbins_per_dim with any skipped dimensions, setting number of bins along skipped dimensions to 1
        skip = np.array([bool(s) for s in skip])
        nbins_per_dim = np.array(nbins_per_dim)
        nbins_per_dim[skip] = 1

        # Total bins is product of all linear bins plus and special bins
        n_total_bins = nbins_per_dim.prod()
        for direct, skip_dim in zip(direction, skip):
            if not skip_dim:
                if direct in [-1, 1]:
                    # 1 lead or lag bin + 1 bottleneck bin
                    n_total_bins += 1 + 1 * bottleneck
                elif direct == 0:
                    # 2 lead/lag bins + 2 bottleneck bins
                    n_total_bins += 2 + 2 * bottleneck
                elif direct == 86:
                    # 0 lead/lag + 2 bottleneck bins
                    n_total_bins += 2 * bottleneck
        return n_total_bins


def map_mab(coords: np.ndarray, mask: np.ndarray, output: List[int], *args, **kwargs) -> List[int]:
    """
    Adaptively place bins based on extrema and bottleneck segments along the progress coordinate.

    Bottleneck segments are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension (unless skipped) and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.

    Parameters
    ----------
    coords : np.ndarray
        An array with pcoord and weight info.
    mask : np.ndarray
        Boolean array to filter out unwanted segments.
    output : list
        The main list that, for each segment, holds the bin assignment.
    *args : list
        Additional arguments.
    **kwargs : dict
        Additional keyword arguments. Contains most of the MAB-needed parameters.

    Returns
    ------
    output : list
        List with bin assignments for each segment.
    """

    # Argument Processing
    nbins_per_dim = kwargs.get("nbins_per_dim")
    ndim = len(nbins_per_dim)
    pca = kwargs.get("pca", False)
    bottleneck = kwargs.get("bottleneck", True)
    direction = kwargs.get("direction", [0] * ndim)
    skip = kwargs.get("skip", [0] * ndim)
    mab_log = kwargs.get("mab_log", False)
    bin_log = kwargs.get("bin_log", False)
    bin_log_path = kwargs.get("bin_log_path", "$WEST_SIM_ROOT/binbounds.log")

    if not np.any(mask):
        return output

    if skip is None:
        skip = [0] * ndim

    allcoords = coords.copy()
    allmask = mask.copy()

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
            isfinal = allcoords[:, ndim + 1].astype(bool)
        else:
            isfinal = np.ones(coords.shape[0], dtype=bool)
        coords = coords[isfinal, :ndim]
        weights = allcoords[isfinal, ndim]
        mask = mask[isfinal]
        splitting = True

    if not np.any(mask):
        coords = allcoords[:, :ndim]
        mask = allmask
        weights = None
        splitting = False

    originalcoords = np.copy(coords)
    if pca and len(output) > 1:
        coords = apply_pca(coords, weights)

    # Computing special bins (bottleneck and boundary bins)
    minlist, maxlist, bottlenecks_forward, bottlenecks_reverse = calculate_bin_boundaries(
        originalcoords, weights, mask, skip, splitting, bottleneck
    )

    if mab_log and report:
        log_mab_stats(minlist, maxlist, direction, skip)

    # Assign segments to bins
    n_bottleneck_filled = bin_assignment(
        allcoords,
        allmask,
        minlist,
        maxlist,
        bottlenecks_forward,
        bottlenecks_reverse,
        nbins_per_dim,
        direction,
        skip,
        splitting,
        bottleneck,
        output,
    )

    # Report MAB bin statistics
    if bin_log and report and westpa.rc.sim_manager.n_iter:
        log_bin_boundaries(
            skip,
            bottleneck,
            direction,
            bin_log_path,
            minlist,
            maxlist,
            nbins_per_dim,
            n_bottleneck_filled,
            bottlenecks_forward,
            bottlenecks_reverse,
        )

    return output


def apply_pca(coords, weights):
    colavg = np.mean(coords, axis=0)
    varcoords = coords - colavg
    covcoords = np.cov(varcoords.T, aweights=weights)
    eigval, eigvec = np.linalg.eigh(covcoords)
    eigvec = eigvec[:, np.argmax(np.abs(eigvec), axis=1)]
    eigvec[:, np.diag(eigvec) < 0] *= -1
    return np.dot(varcoords, eigvec)


def calculate_bin_boundaries(coords, weights, mask, skip, splitting, bottleneck):
    """
    This function calculates minima, maxima, and bottleneck segments.
    """
    skip = np.array([bool(s) for s in skip])

    # Initialize lists to hold minima and maxima along each dimension
    minlist, maxlist = [], []
    # Initialize lists to hold bottleneck segments along each dimension
    bottlenecks_forward, bottlenecks_reverse = [None] * len(coords[0]), [None] * len(coords[0])
    # number of unmasked coords
    n_coords = mask.sum()
    # Grabbing all unmasked coords and weights
    unmasked_coords = coords[mask, :]
    unmasked_weights = weights[mask] if weights is not None else None
    # Replace any zero weights with non-zero values so that log(weight) is well-defined
    if unmasked_weights is not None:
        unmasked_weights[unmasked_weights == 0] = 10**-323
    # Looping over each dimension of progress coordinate, even those being skipped
    for n in range(len(coords[0])):
        # We calculate the min and max pcoord along each dimension (boundary segments) even if skipping
        maxlist.append(np.max(coords[mask, n]))
        minlist.append(np.min(coords[mask, n]))
        # Now we calculate the bottleneck segments
        if splitting and bottleneck and not skip[n]:
            bottlenecks_forward[n], bottlenecks_reverse[n] = detect_bottlenecks(unmasked_coords, unmasked_weights, n_coords, n)

    return minlist, maxlist, bottlenecks_forward, bottlenecks_reverse


def detect_bottlenecks(unmasked_coords, unmasked_weights, n_coords, n):
    """
    Detect the bottleneck segments along the given coordinate n, this uses the weights
    """
    # Grabbing all unmasked coords in current dimension, plus corresponding weights
    # Sort by current dimension in coord, smallest to largest
    sorted_indices = unmasked_coords[:, n].argsort(kind='stable')

    # Grab sorted coords and weights
    coords_srt = unmasked_coords[sorted_indices, :]
    weights_srt = unmasked_weights[sorted_indices]

    # Also sort in reverse order for opposite direction
    coords_srt_flip = np.flipud(coords_srt)
    weights_srt_flip = np.flipud(weights_srt)

    # Initialize the max directional differences along current dimension as None (these may not be updated)
    bottleneck_coords, bottleneck_coords_flip = None, None
    maxdiff, maxdiff_flip = -np.inf, -np.inf

    # Looping through all non-boundary coords
    # Compute the cumulative weight on either side of each non-boundary walker
    for i in range(1, n_coords - 1):
        # Summing up weights of all walkers ahead of current walker along current dim in both directions
        cumulative_prob = np.sum(weights_srt[i + 1 :])
        cumulative_prob_flip = np.sum(weights_srt_flip[i + 1 :])
        # Compute the difference of log cumulative weight of current walker and all walkers ahead of it (Z im the MAB paper)
        # We use the log as weights vary over many orders of magnitude
        # Note a negative Z indicates the cumulative weight ahead of the current walker is larger than the weight of the current walker,
        # while a positive Z indicates the cumulative weight ahead of the current walker is smaller, indicating a barrier
        Z = np.log(weights_srt[i]) - np.log(cumulative_prob)
        Z_flip = np.log(weights_srt_flip[i]) - np.log(cumulative_prob_flip)
        # Update ALL coords of the current walker into bottlenecks_forward if it is largest
        # This way we uniquely identify a walker by its full set of coordinates
        if Z > maxdiff:
            bottleneck_coords = coords_srt[i, :]
            maxdiff = Z
        if Z_flip > maxdiff_flip:
            bottleneck_coords_flip = coords_srt_flip[i, :]
            maxdiff_flip = Z_flip
    return bottleneck_coords, bottleneck_coords_flip


def log_mab_stats(minlist, maxlist, direction, skip):
    westpa.rc.pstatus("################ MAB stats ################")
    westpa.rc.pstatus(f"minima in each dimension:      {minlist}")
    westpa.rc.pstatus(f"maxima in each dimension:      {maxlist}")
    westpa.rc.pstatus(f"direction in each dimension:   {direction}")
    westpa.rc.pstatus(f"skip in each dimension:        {skip}")
    westpa.rc.pstatus("###########################################")
    westpa.rc.pflush()


def bin_assignment(
    coords,
    mask,
    minlist,
    maxlist,
    bottlenecks_forward,
    bottlenecks_reverse,
    nbins_per_dim,
    direction,
    skip,
    splitting,
    bottleneck,
    output,
):
    """
    Assign segments to bins based on the minima, maxima, and
    bottleneck segments along the progress coordinate.
    """
    # Update nbins_per_dim with any skipped dimensions, setting number of bins along skipped dimensions to 1
    skip = np.array([bool(s) for s in skip])
    nbins_per_dim = np.array(nbins_per_dim)
    nbins_per_dim[skip] = 1
    direction = np.array(direction)

    ndim = len(nbins_per_dim)
    n_bottleneck_filled = 0

    # Boolean arrays that track use of special bins along each dimension
    skip_bneck_fwd = np.array([d == -1 if bottleneck else True for d in direction]) + skip
    skip_bneck_rev = np.array([d == 1 if bottleneck else True for d in direction]) + skip
    skip_lead = np.array([d in [86, -1] for d in direction]) + skip
    skip_lag = np.array([d in [86, 1] for d in direction]) + skip

    # List of dimensions that are not skipped
    active_dims = np.array([n for n in range(ndim) if not skip[n]])

    # Compute the boundary bin ID offsets
    # In forward direction, this is all the linear bins
    boundary_bin_id_offset_fwd = nbins_per_dim.prod()
    # In reverse, we add the number of forward boundary bins to the offset
    boundary_bin_id_offset_rev = boundary_bin_id_offset_fwd + (~skip_lead).sum()

    # Compute the bottleneck bin ID offsets
    # In forward direction, bin IDs are offset by all linear and boundary bins
    bneck_bin_id_offset_fwd = boundary_bin_id_offset_rev + (~skip_lag).sum()
    # In reverse, we add the number of forward bottleneck bins to the offset
    bneck_bin_id_offset_rev = bneck_bin_id_offset_fwd + (~skip_bneck_fwd).sum()

    # Bin assignment loop over all walkers
    for i in range(len(output)):
        # Skip masked walkers, these walkers bin IDs are unchanged
        if not mask[i]:
            continue
        # Initialize bin ID and special tracker for current coord
        # The special variable indicates a boundary or bottleneck walker (not assigned to the linear space)
        bin_id, special = 0, False

        # Searching for bottleneck bins first
        if splitting and bottleneck:
            for n in active_dims:
                # Grab coord(s) of current walker
                coord = coords[i][:ndim]
                # Assign bottlenecks, taking directionality into account
                # Check both directions when using 0 or 86
                # Note: 86 implies no leading or lagging bins, but does add bottlenecks for *both* directions when bottleneck is enabled
                # Note: All bottleneck bins will typically be filled unless a walker is simultaneously in bottleneck bins along multiple dimensions
                # or there are too few walkers to compute free energy barriers
                if (coord == bottlenecks_forward[n]).all() and not skip_bneck_fwd[n]:
                    bin_id = bneck_bin_id_offset_fwd + n - skip_bneck_fwd[:n].sum()
                    special = True
                    n_bottleneck_filled += 1
                    break
                elif (coord == bottlenecks_reverse[n]).all() and not skip_bneck_rev[n]:
                    bin_id = bneck_bin_id_offset_rev + n - skip_bneck_rev[:n].sum()
                    special = True
                    n_bottleneck_filled += 1
                    break

        # Now check for boundary walkers, taking directionality into account
        # This should only be done after fully checking for bottleneck walkers
        if splitting and not special:
            for n in active_dims:
                # Grab coord of current walker along current dimension
                coord = coords[i, n]
                if (coord == maxlist[n]) and not skip_lead[n]:
                    bin_id = boundary_bin_id_offset_fwd + n - skip_lead[:n].sum()
                    special = True
                    break
                elif (coord == minlist[n]) and not skip_lag[n]:
                    bin_id = boundary_bin_id_offset_rev + n - skip_lag[:n].sum()
                    special = True
                    break

        # Now check for linear bin walkers
        if not special:
            # Again we loop over the dimensions
            # Note: no need to worry about skipping as we've already set all skipped dimensions to 1 bin
            for n in range(ndim):
                coord = coords[i][n]
                nbins = nbins_per_dim[n]
                minp = minlist[n]
                maxp = maxlist[n]

                # Generate the bins along this dimension
                bins = np.linspace(minp, maxp, nbins + 1)

                # Assign walker to a bin along this dimension
                bin_number = np.digitize(coord, bins) - 1  # note np.digitize is 1-indexed

                # Sometimes the walker is exactly at the max/min value,
                # which would put it in the next bin
                if bin_number == nbins:
                    bin_number -= 1
                elif bin_number == -1:
                    bin_number = 0
                elif bin_number > nbins or bin_number < -1:
                    raise ValueError("Walker out of boundary.")

                # Assign to bin within the full dimensional space
                bin_id += bin_number * np.prod(nbins_per_dim[:n])

        # Output is the main list that, for each segment, holds the bin assignment
        output[i] = bin_id
    return n_bottleneck_filled


def log_bin_boundaries(
    skip,
    bottleneck,
    direction,
    bin_log_path,
    minlist,
    maxlist,
    nbins_per_dim,
    n_bottleneck_filled,
    bottlenecks_forward,
    bottlenecks_reverse,
):
    ndim = len(nbins_per_dim)
    skip = np.array([bool(s) for s in skip])
    active_dims = np.array([n for n in range(ndim) if not skip[n]])
    max_bottleneck = np.sum([1 if direction[n] in [-1, 1] else 2 for n in active_dims]) if bottleneck else 0
    with open(expandvars(bin_log_path), 'a') as bb_file:
        # Iteration Number
        bb_file.write(f'Iteration: {westpa.rc.sim_manager.n_iter}\n')
        bb_file.write('MAB linear bin boundaries: ')
        for n in range(ndim):
            # Write binbounds per dim
            bb_file.write(f'{np.linspace(minlist[n], maxlist[n], nbins_per_dim[n] + 1)}\t')
        # Min/Max pcoord
        bb_file.write(f'\nLagging pcoord in each dimension: {minlist}\n')
        bb_file.write(f'Leading pcoord in each dimension: {maxlist}\n')
        # Bottlenecks bins exist
        if bottleneck:
            bb_file.write(f'Number of bottleneck bins filled: {n_bottleneck_filled} / {max_bottleneck}\n')
            for n in active_dims:
                if direction[n] in [0, 1, 86]:
                    bb_file.write(f'Dimension {n} forward bottleneck walker at: {[bottlenecks_forward[n]]}\n')
                if direction[n] in [0, -1, 86]:
                    bb_file.write(f'Dimension {n} backward bottleneck walker at: {[bottlenecks_reverse[n]]}\n')
            bb_file.write('\n')
        else:
            bb_file.write('\n')
