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
        ndim = len(nbins_per_dim)
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
    ncoords = len(output)

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

    varcoords = np.copy(coords)
    originalcoords = np.copy(coords)
    if pca and ncoords > 1:
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
            for j in range(ncoords):
                coords[j][i] = np.dot(varcoords[j], eigvec[:, i])

    maxlist = []
    minlist = []
    difflist = []
    difflist_flip = []
    # Looping over each dimension of progress coordinate, even those being skipped
    n_coords = mask.sum() # number of unmasked coords
    for n in range(ndim):
        # We calculate the min and max pcoord along each dimension (boundary segments) even if skipping
        maxcoord = np.max(coords[mask, n])
        mincoord = np.min(coords[mask, n])
        maxlist.append(maxcoord)
        minlist.append(mincoord)

        # Detect the bottleneck segments, this uses the weights
        if splitting:
            # Grabbing all unmasked coords in current dimension, plus corresponding weights
            coords_weights = np.column_stack((originalcoords[mask, n], weights[mask]))
            # Sort by coord smallest to largest
            sorted_indices = coords_weights[:, 0].argsort()
            coords_weights = coords_weights[sorted_indices]
            # Replace any zero weights with non-zero values so that log(weight) is well-defined
            coords_weights[:, 1][coords_weights[:, 1] == 0] = 10**-323
            coords_weights_flip = np.flipud(coords_weights)
            # Initialize the max directional differences along current dimension as None (these may not be updated)
            difflist.append(None)
            difflist_flip.append(None)
            maxdiff = -np.inf
            maxdiff_flip = -np.inf

            # Looping through all non-boundary coords
            # Compute the cumulative weight on either side of each non-boundary walker
            for i in range(1, n_coords - 1):
                cumulative_prob = 0
                cumulative_prob_flip = 0
                # Simultaneously looping through all walkers above/below walkers at position j from the top/bottom of extrema
                # Summing up weights of all walkers ahead of current walker in both directions
                j = i + 1
                while j < n_coords:
                    cumulative_prob = cumulative_prob + coords_weights[j][1]
                    cumulative_prob_flip = cumulative_prob_flip + coords_weights_flip[j][1]
                    j += 1
                # Compute the difference of log cumulative weight of current walker and all walkers ahead of it (Z im the MAB paper)
                # We use the log as weights vary over many orders of magnitude
                # Note a negative Z indicates the cumulative weight ahead of the current walker is larger than the weight of the current walker,
                # while a positive Z indicates the cumulative weight ahead of the current walker is smaller, indicating a barrier
                Z = np.log(coords_weights[i][1]) - np.log(cumulative_prob)
                if Z > maxdiff:
                    # Update the current coord into difflist if it is largest
                    difflist[n] = coords_weights[i][0]
                    maxdiff = Z
                # Compute the same in the reverse direction
                Z_flip = np.log(coords_weights_flip[i][1]) - np.log(cumulative_prob_flip)
                if Z_flip > maxdiff_flip:
                    difflist_flip[n] = coords_weights_flip[i][0]
                    maxdiff_flip = Z_flip

    if mab_log and report:
        westpa.rc.pstatus("################ MAB stats ################")
        westpa.rc.pstatus("minima in each dimension:      {}".format(minlist))
        westpa.rc.pstatus("maxima in each dimension:      {}".format(maxlist))
        westpa.rc.pstatus("direction in each dimension:   {}".format(direction))
        westpa.rc.pstatus("skip in each dimension:        {}".format(skip))
        westpa.rc.pstatus("###########################################")
        westpa.rc.pflush()

    # Update nbins_per_dim with any skipped dimensions, setting number of bins along skipped dimensions to 1
    skip = np.array([bool(s) for s in skip])
    nbins_per_dim = np.array(nbins_per_dim)
    nbins_per_dim[skip] = 1

    # Assigning segments to bins. 
    # First we compute the offset for boundary bin IDs
    # This will just be the product of all non-skipped static bins
    boundary_bin_id_offset = nbins_per_dim.prod()

    # Now compute the bottleneck bin ID offset
    # The bottleneck IDs are offset by the static bins plus the boundary bins, which exist along all non-skipped dimensions.
    bottleneck_bin_id_offset = boundary_bin_id_offset
    for n in range(ndim):
        if not skip[n]:
            # for single direction, 1 boundary walker
            if direction[n] in [1,-1]:
                bottleneck_bin_id_offset += 1
            # 2 boundary walkers with 0 direction
            elif direction[n] == 0:
                bottleneck_bin_id_offset += 2
            # for 86 direction, no boundary walkers so no change in offset
            elif direction[n] == 86:
                continue
    
    # Bin assignment loop over all walkers
    n_bottleneck_filled = 0 # Tracks number of bottleneck bins filled
    for i in range(ncoords):
        if not allmask[i]:
            # Skip masked walkers
            continue
        # The special variable indicates a boundary or bottleneck walker (not assigned to the linear space)
        # Initialize special as false for current walker
        special = False
        # Initialize bin ID
        bin_id = 0
        # Searching for special split bins from pre-computed lists
        if splitting:
            # Loop over dimensions
            for n in range(ndim):
                if skip[n]:
                    # Skip this dimension
                    continue
                # Grab coord of current walker along current dimension
                coord = allcoords[i][n]
                # Assign bottlenecks, taking directionality into account
                # Check both directions when using 0 or 86
                # Note: 86 implies no leading or lagging bins, but does add bottlenecks for *both* directions when bottleneck is enabled
                # Note: All bottleneck bins will typically be filled unless a walker is simultaneously in bottleneck bins along multiple dimensions
                # or there are too few walkers to compute free energy barriers
                if bottleneck:
                    if ( (direction[n] == -1) and (coord == difflist_flip[n]) ) or ( (direction[n] in [1,0,86]) and (coord == difflist[n]) ):
                        bin_id = bottleneck_bin_id_offset + n - skip[:n].sum()
                        special = True
                        n_bottleneck_filled += 1
                        break
                    elif ( (direction[n] in [0,86]) and (coord == difflist_flip[n]) ):
                        bin_id = bottleneck_bin_id_offset + n + 1 - skip[:n].sum()
                        special = True
                        n_bottleneck_filled += 1
                        break
                # Now check for boundary walkers, taking directionality into account
                if direction[n] == 86:
                    # 86 uses no lead/lag splitting
                    # But do not *break* the loop as we still need to check along other dimensions for bottleneck walkers
                    continue
                elif ( (direction[n] in [0,-1]) and (coord == minlist[n]) ) or ( (direction[n] == 1) and (coord == maxlist[n]) ):
                    bin_id = boundary_bin_id_offset + n - skip[:n].sum()
                    special = True
                    break
                elif ( (direction[n] == 0) and (coord == maxlist[n]) ):
                    bin_id = boundary_bin_id_offset + n + 1 - skip[:n].sum()
                    special = True
                    break

        # The following loop only applies if the current walker is not assigned to a special bin
        # i.e. this is for the "linear" portion
        if not special:
            # Again we loop over the dimensions
            # Note we do not need to worry about skipping as we've already set all skipped dimensions to have only 1 bin along that dimension
            for n in range(ndim):
                coord = allcoords[i][n]
                nbins = nbins_per_dim[n]
                minp = minlist[n]
                maxp = maxlist[n]

                # Generate the bins along this dimension
                bins = np.linspace(minp, maxp, nbins + 1)

                # Assign walker to a bin along this dimension
                bin_number = np.digitize(coord, bins) - 1 # note np.digitize is 1-indexed

                # Check for problem values
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
                
                # Assign to bin within the full dimensional space
                bin_id += bin_number * np.prod(nbins_per_dim[:n])

        # Output is the main list that, for each segment, holds the bin assignment
        output[i] = bin_id

    # Report MAB bin statistics
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
                    bb_file.write(f'bottleneck pcoord: {difflist_flip} {difflist}\n\n')
                else:
                    bb_file.write('\n')

    return output
