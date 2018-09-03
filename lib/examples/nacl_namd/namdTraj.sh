#!/bin/bash


# --------------------------------
# NAMD Trajectory Tool for WESTPA
# --------------------------------
# 
# Written by Anthony Bogetti on 30.08.18
# 
# This script will stitch together a trajectory file from your NAMD-WESTPA
# simulation that can be viewed in VMD or another molecular dynmaics 
# visualization software.  Run this script with the command ./namdTraj.sh
# from the same directory where the west.h5 file from your WESTPA simulation
# is located.  The results of this analysis will be stored in a new folder
# called trajAnalysis as the file trace.dcd.  Load trace.dcd into VMD to 
# visualize the trajectory.  As a note, you will need to have your computer
# configured to run w_succ from the WESTPA software package and catdcd from 
# the VMD distribution site.  Note that catdcd is NOT included in the VMD
# software package.  You can download the catdcd4.0 files from the VMD site.
# Extract the tarball and you will see a multiple directories specifying
# different operating systems (LINUX, LINUXAMD64, MACOSX etc.).  Find the 
# catdcd executable file for your system and add the containing folder (it
# will be named catdcd4.0) to your path.  Then you will be able to run the 
# command catdcd from the command line.  But you won't really need to do that,
# just make sure that this script can call the catdcd command when run.

# The variables defined below are the name of the new analysis directory that
# will be created and the name of an intermediate file in the process of 
# stitching together the trajectory file.
dir=trajAnalysis
file=iters.txt

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
#cp prep/2_solvate/nacl_solvated.gro $dir

# The output from w_succ above is moved into the analysis directory.
mv succ.txt $dir 
cd $dir

# The first six lines of succ.txt are removed, as they are just some notes
# written by the w_succ program and are not relevant for the analysis.
cat succ.txt | tail -n +7 > iters.txt

# Since catdcd takes as its input a filepath to the trajectory file and we want
# to stitch together multiple trajectories, we want to create a giant variable consisting
# of each filepath separated by a single space.  To do this, the filepaths are appended
# into a text file (filestrings.txt) and then catenated into the variable name later on.
# The variable is then called in the catdcd command as input.

# Please note that while the iteration and segment IDs here are padded to six
# digits with zeroes, the length of this number is specified in the west.cfg file
# in the main WESTPA simuation directory and can be changed by the user.  If you
# ran the simulation with more than 100000 iterations or segments and adjusted this
# parameter in the west.cfg file you will need to adjust it here too.  For 99% of
# users, however, the following should work just fine.

while read file; do
	iter=$(echo $file | awk '{print $1}')
	seg=$(echo $file | awk '{print $2}')
	filestring='../traj_segs/'$(printf "%06d" $iter)'/'$(printf "%06d" $seg)'/''seg.dcd' 
        echo -n "$filestring " >> filestrings.txt	
done < "$file"

inputs=$(cat 'filestrings.txt')

# The text normally displayed to the terminal is written
# to the file traj.log.

catdcd -o trace.dcd $inputs > traj.log

# The intermediary files are removed to clean up the analysis directory.
rm succ.txt iters.txt filestrings.txt
cd ..
