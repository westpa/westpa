#!/usr/bin/env python
# rmsd.py
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
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
    """ Calculate pcoord (RMSD) using MDAnalysis and save results to file specified
    in get_pcoord.sh/runseg.sh. Here the filename is rmsd.dat, but if you were
    calculating something else such as distance you could change the filename to
    distance.dat instead. Just make sure to change the filename both in this
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

    # Create Universe objects for initial structure and segment
    # structure. (args: topology file, trajectory file)
    # If segment file is Amber netCDF trajectory, it must have extension
    # ".ncdf" to be recognized automatically by MDAnalysis. The filetype can
    # also be specified using the optional "format" argument.
    init_u = mda.Universe(toppath, refpath, format="RESTRT")
    seg_u = mda.Universe(toppath, mobpath, format=str(FORM))

    # Create c-alpha AtomGroups.
    init_cAlpha = init_u.select_atoms("name CA")
    seg_cAlpha = seg_u.select_atoms("name CA")

    # Calculate RMSD (relative to initial structure) at each time step.
    R = RMSD(seg_cAlpha, init_cAlpha, select = 'name CA', center=True, superposition=True)
    R.run()

    # Write RMSD to output file.
    if FORM == "RESTRT":
        numpy.savetxt("rmsd.dat", R.rmsd[:,2])
    else:
        numpy.savetxt("rmsd.dat", R.rmsd[:,2])

def main():
    # Get arguments from the caller and pass to calc_pcoord().
    refpath, toppath, mobpath, FORM = parse_arguments()
    calc_pcoord(refpath, toppath, mobpath, FORM)

if __name__ == "__main__":
    main()
