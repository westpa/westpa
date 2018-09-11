#!/bin/bash


# --------------------------------
# GROMACS Trajectory Tool for WESTPA
# --------------------------------
# 
# Written by Anthony Bogetti on 28.08.18
# 
# This script will stitch together a trajectory file from your GROMACS-WESTPA
# simulation that can be viewed in VMD or another molecular dynmaics 
# visualization software.  Run this script with the command ./gmxTraj.sh
# from the same directory where the west.h5 file from your WESTPA simulation
# is located.  The results of this analysis will be stored in a new folder
# called trajAnalysis as the file trace.trr.  Load trace.trr into VMD to 
# visualize the trajectory.  As a note, you will need to have your computer
# configured to run w_succ from the WESTPA software package and cpptraj from 
# the Amber software package.  Though, if the simulation has completed successfully,
# these commands will most likely be ready to run.


# The variables defined below are the name of the new analysis directory that
# will be created and the name of an intermediate file in the process of 
# stitching together the trajectory file.
dir=trajAnalysis
file=path.txt

# w_succ, one of the WESTPA analysis tools, is called here to look through the h5 file
# of the relevant simulation and print the iteration and segment IDs of all
# successful trajectories.  This information is stored and will be used later on.
w_succ > succ.txt

# The analysis directory is then made and the parameter file for the NaCl system is
# copied into it.  All analysis will take place within this directory.
if [ -d "$dir" ]; then
  rm -r $dir
fi
mkdir $dir

# The topology file of the NaCl solvated system is copied into the trajectory analysis
# directory.
cp ./prep/2_solvate/nacl_solvated.gro $dir

# The output from w_succ above is moved into the analysis directory.
mv succ.txt $dir 
cd $dir

# The first six lines of succ.txt are removed, as they are just some notes
# written by the w_succ program and are not relevant for the analysis.
cat succ.txt | tail -n +7 > iters.txt

# The input file for cpptraj is prepared.
if [ -f "cpptraj.in" ]; then
  rm "cpptraj.in"
fi

# The first successful trajectory is chosen and the iteration and segment
# of that trajectory are assigned to variables.
siter=$(head iters.txt -n 1 | awk '{ print $1 }')
sseg=$(head iters.txt -n 1 | awk '{ print $2 }')

# We'll move back to the main simulation directory to prepare to run w_trace.
cd ../

# w_trace is run, generating a history of the successful trajectory specified
# above.  This will create a list of all of the iterations and segments that
# need to be stitched together to create a smooth, viewable, successful trajectory.
w_trace $siter:$sseg

# Output files from w_trace are moved into the trajAnalysis directory.
mv $(echo 'traj_'$siter'_'$sseg'_trace.txt') $dir
mv trajs.h5 $dir
cd $dir

# The first few lines of the output of w_trace are removed (including the 
# initial state of the system, which doesn't have an iter:seg ID)
cat $(echo 'traj_'$siter'_'$sseg'_trace.txt') | tail +9 > path.txt


# Now, the file listing all of the successful trajectory's historical iterations
# and segments is read line by line and the iteration IDs and segment IDs
# are added to a variable that specifies the path to the coordinate file of
# that successful trajectory.  This variable is then called as input to the
# trajcat command.


# Please note that while the iteration and segment IDs here are padded to six
# digits with zeroes, the length of this number is specified in the west.cfg file
# in the main WESTPA simuation directory and can be changed by the user.  If you
# ran the simulation with more than 100000 iterations or segments and adjusted this
# parameter in the west.cfg file you will need to adjust it here too.  For 99% of
# users, however, the following should work just fine.

while read file; do
	iter=$(echo $file | awk '{print $1}')
	seg=$(echo $file | awk '{print $2}')
	filestring='../traj_segs/'$(printf "%06d" $iter)'/'$(printf "%06d" $seg)'/''seg.trr' 
        echo -n "$filestring " >> filestrings.txt	
done < "path.txt"

inputs=$(cat 'filestrings.txt')

# The text normally displayed to the terminal is written
# to the file traj.log.

gmx trjcat -cat -f $inputs -o trace.trr > traj.log

echo Trajectory file creation is complete.
echo To view your trajectory, load the NaCl parameter file into VMD followed by the trace.trr file, both located in the trajAnalysis directory.

# The intermediary files are removed to clean up the analysis directory.
rm succ.txt iters.txt path.txt $(echo 'traj_'$siter'_'$sseg'_trace.txt') trajs.h5 filestrings.txt 
cd ..

