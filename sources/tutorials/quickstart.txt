.. _quickstart_guide:

Quickstart Guide
================

Overview
--------

This tutorial is intended to get you started quickly in running a WESTPA
simulation. For more advanced options, please check out the tutorials
that are available on our Sphinx website.

Installing WESTPA
------------------

To install WESTPA, just follow the three steps below.

1) Clone the WESTPA repo.

::

    # git clone https://github.com/westpa/westpa.git
    Cloning into 'westpa'...
    remote: Counting objects: 6377, done.
    remote: Total 6377 (delta 0), reused 0 (delta 0), pack-reused 6377
    Receiving objects: 100% (6377/6377), 8.85 MiB | 0 bytes/s, done.
    Resolving deltas: 100% (4015/4015), done.
    Checking connectivity... done.

2) Make sure that you have the latest Anaconda python binary (`Anaconda
distribution <https://www.continuum.io/downloads>`__)

::

    #which python
    $HOME/apps/anaconda/bin/python

3) Go to the westpa folder and run setup.sh.

::

    #cd westpa
    #./setup.sh
    $HOME/westpa/lib $HOME/test_qsg/westpa
    $HOME/westpa/lib/blessings $HOME/westpa/lib $HOME/westpa
    From git://github.com/erikrose/blessings
     * branch            master     -> FETCH_HEAD
    Updating d3ba51c..979988c
    ...

Running a WESTPA simulation
----------------------------

Now let's run a quick WESTPA simulation of Na+/Cl- association using the
NAMD dynamics engine. To run this simulation, follow the steps below.

1) Check that you have the binary for the NAMD dynamics engine in your
PATH variable

::

    #which namd2

2) Export the WEST\_ROOT environment variable

::

    #export WEST_ROOT=$(pwd)

3) Head over to the example simulation directory and initialize the
simulation

::

    #cd lib/examples/nacl_namd
    #./init.sh
    simulation nacl_namd root is $HOME/westpa/lib/examples/nacl_namd
    Creating HDF5 file '$HOME/westpa/lib/examples/nacl_namd/west.h5'
    2 target state(s) present
    Calculating progress coordinate values for basis states.
    1 basis state(s) present
    Preparing initial states
            
            Total bins:            22
            Initial replicas:      24 in 1 bins, total weight = 1
            Total target replicas: 528
           
    Simulation prepared.
    1 of 22 (4.545455%) active bins are populated
    per-bin minimum non-zero probability:       1
    per-bin maximum probability:                1
    per-bin probability dynamic range (kT):     0
    per-segment minimum non-zero probability:   0.0416667
    per-segment maximum non-zero probability:   0.0416667
    per-segment probability dynamic range (kT): 0
    norm = 1, error in norm = 0 (0*epsilon)

4) Run the WESTPA simulation

::

    #./run.sh
    simulation nacl_namd root is $HOME/westpa/lib/examples/nacl_namd

Voila! Your simulation results can be found in the west.h5 file.
