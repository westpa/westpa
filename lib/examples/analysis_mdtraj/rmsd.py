#!/usr/bin/env python
# rmsd.py
import mdtraj as md
import numpy
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', dest='refpath', required=True)
    parser.add_argument('--top', dest='toppath', required=True)
    parser.add_argument('--mob', dest='mobpath', required=True)
    parser.add_argument('--for', dest='FORM', required=True)
    args = parser.parse_args()
    return args.refpath, args.toppath, args.mobpath, args.FORM

def calc_pcoord(refpath, toppath, mobpath, FORM):
    """ Calculate pcoord (RMSD) using MDTraj and save results to file specified
    in get_pcoord.sh/runseg.sh. Here the filename is rmsd.dat, but if you were
    calculating somebody else like a simple distance you could change the filename
    to distance.dat instead. Just make sure to change the filename both in this
    script and in get_pcoord.sh/runseg.sh.

    Parameters:
        refpath (str): path to initial state coordinate file.
        toppath (str): path to topology file.
        mobpath (str): path to trajectory file.
        FORM (str): indicates whether we're evaluating a basis/initial state or not.
            If we are evaluating an initial/basis state (ie. if the script is
            called from get_pcoord.sh) then FORM = 'RESTRT', and we check to
            make sure our pcoord is a numpy array with shape (1,). Otherwise,
            the pcoord is a numpy array with shape = (pcoord_len, pcoord_ndim)
            as specified in west.cfg.
    """

    # Load the reference crystal and the trajectory
    # Use the load_netcdf() function so MDtraj knows it is a netcdf file.
    crystal = md.load_netcdf(refpath, top=toppath)
    traj = md.load_netcdf(mobpath, top=toppath)

    # Get a list of CA indices from the topology file.
    CA_indices = crystal.topology.select("name == CA")

    # Calculate the rmsd of the trajectory relative to the crystal, using only
    # the C-Alpha atoms for the calculation (we must specify this as there is
    # explicit solvent present in the simulation.)
    # The rmsd() function takes an optional third int argument which refers to
    # the frame in the reference to measure distances to. By default, the frame
    # is set to 0. A general form of the function is:
    # MDTraj.rmsd(target, reference, frame=0) which returns a numpy array
    rmsd = md.rmsd(traj, crystal, atom_indices=CA_indices)

    # Write RMSD to output file.
    if FORM == "RESTRT":
	    # We only need the last value in the array.
	    rmsd = numpy.array(rmsd[-1])
	    # WESTPA expects a 1x1 array, so we must correct the shape if needed.
        if rmsd.ndim == 0:
	    rmsd.shape = (1,)
        numpy.savetxt("rmsd.dat", rmsd)
    else:
        numpy.savetxt("rmsd.dat", rmsd)

def main():
    # Get arguments from the caller and pass to calc_pcoord().
    refpath, toppath, mobpath, FORM = parse_arguments()
    calc_pcoord(refpath, toppath, mobpath, FORM)

if __name__ == "__main__":
    main()
