HDF5 File Schema
================

WESTPA stores all of its simulation data in the cross-platform, self-describing
`HDF5 <http://www.hdfgroup.org/HDF5>`_ file format. This file format can be
read and written by a variety of languages and toolkits, including C/C++,
Fortran, Python, Java, and `Matlab
<http://www.mathworks.com/help/matlab/ref/hdf5read.html>`_ so that analysis of
weighted ensemble simulations is not tied to using the WESTPA framework. HDF5
files are organized like a filesystem, where arbitrarily-nested groups (i.e.
directories) are used to organize datasets (i.e. files). The excellent `HDFView
<http://www.hdfgroup.org/hdf-java-html/hdfview/>`_ program may be used to
explore WEST data files.

The canonical file format reference for a given version of the WEST code is
described in `src/west/data_manager.py
<https://github.com/westpa/westpa/blob/master/src/west/data_manager.py>`_.

Overall structure
-----------------

::

    /
        #ibstates/
            index
            naming
                bstate_index
                bstate_pcoord
                istate_index
                istate_pcoord
        #tstates/
            index
        bin_topologies/
            index
            pickles
        iterations/
            iter_XXXXXXXX/\|iter_XXXXXXXX/
                auxdata/
                bin_target_counts
                ibstates/
                    bstate_index
                    bstate_pcoord
                    istate_index
                    istate_pcoord
                pcoord
                seg_index
                wtgraph
            ...
        summary

The root group (/)
------------------

The root of the WEST HDF5 file contains the following entries (where a
trailing "/" denotes a group):

=============== ======================= =======================================
Name            Type                    Description
=============== ======================= =======================================
ibstates/       Group                   Initial and basis states for this
                                        simulation
tstates/        Group                   Target (recycling) states for this
                                        simulation; may be empty
bin_topologies/ Group                   Data pertaining to the binning scheme
                                        used in each iteration
iterations/     Group                   Iteration data
summary         Dataset (1-dimensional, Summary data by iteration
                compound)
=============== ======================= =======================================

The iteration summary table (/summary)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=============== ===============================================================
Field           Description
=============== ===============================================================
n_particles     the total number of walkers in this iteration
norm            total probability, for stability monitoring
min_bin_prob    smallest probability contained in a bin
max_bin_prob    largest probability contained in a bin
min_seg_prob    smallest probability carried by a walker
max_seg_prob    largest probability carried by a walker
cputime         total CPU time (in seconds) spent on propagation for this
                iteration
walltime        total wallclock time (in seconds) spent on this iteration
binhash         a hex string identifying the binning used in this iteration
=============== ===============================================================

Per iteration data (/iterations/iter_XXXXXXXX)
-----------------------------------------------

Data for each iteration is stored in its own group, named according to the
iteration number and zero-padded out to 8 digits, as in
``/iterations/iter_00000001`` for iteration 1. This is done solely for
convenience in dealing with the data in external utilities that sort output by
group name lexicographically. The field width is in fact configurable via the
``iter_prec`` configuration entry under ``data`` section of the WESTPA
configuration file.

The HDF5 group for each iteration contains the following elements:

=================== ======================= ===================================
Name                Type                    Description
=================== ======================= ===================================
auxdata/            Group                   All user-defined auxiliary data0
                                            sets
bin_target_counts   Dataset (1-dimensional) The per-bin target count for the
                                            iteration
ibstates/           Group                   Initial and basis state data for
                                            the iteration
pcoord              Dataset (3-dimensional) Progress coordinate data for the
                                            iteration stored as a (num of
                                            segments, pcoord_len, pcoord_ndim)
                                            array
seg_index           Dataset (1-dimensional, Summary data for each segment
                    compound)               
wtgraph             Dataset (1-dimensional)
=================== ======================= ===================================

The segment summary table (/iterations/iter_XXXXXXXX/seg_index)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=============== ===============================================================
Field           Description
=============== ===============================================================
weight          Segment weight
parent_id       Index of parent
wtg_n_parents
wtg_offset
cputime         Total cpu time required to run the segment
walltime        Total walltime required to run the segment
endpoint_type
status
=============== ===============================================================

Bin Topologies group (/bin_topologies)
---------------------------------------

Bin topologies used during a WE simulation are stored as a unique hash
identifier and a serialized ``BinMapper`` object in `python pickle
<http://docs.python.org/2/library/pickle.html>`_ format. This group contains
two datasets:

- ``index``: Compound array containing the bin hash and pickle length
- ``pickle``: The pickled ``BinMapper`` objects for each unique mapper stored
  in a (num unique mappers, max pickled size) array
