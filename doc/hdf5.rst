HDF5 file structure
===================

WESTPA stores simulation data in the cross-platform, self-describing
`HDF5 <http://www.hdfgroup.org/HDF5>`_ file format. This file format can be
read and written by a variety of languages and toolkits, including C/C++,
Fortran, Python, Java, and `Matlab
<http://www.mathworks.com/help/matlab/ref/hdf5read.html>`_ so that analysis of
weighted ensemble simulations is not tied to using the WESTPA framework. HDF5
files are organized like a filesystem, where arbitrarily-nested groups (i.e.
directories) are used to organize datasets (i.e. files). The excellent `HDFView
<http://www.hdfgroup.org/hdf-java-html/hdfview/>`_ program may be used to
explore WEST data files.

The canonical file format reference for a given version of WESTPA is
described in `src/westpa/core/data_manager.py
<https://github.com/westpa/westpa/blob/develop/src/westpa/core/data_manager.py>`_.

Overall structure
-----------------

::

    / (root)
    ├── @westpa_fileformat_version
    ├── @westpa_iter_prec
    │
    ├── iterations/
    │   ├── iter_00000001/
    │   │   ├── auxdata/ (optional)
    │   │   ├── final_states
    │   │   ├── initial_states
    │   │   ├── pcoord
    │   │   ├── seg_index
    │   │   └── wtgraph
    │   │
    │   └── iter_00000002/
    │       ...
    │
    └── summary


Root group (/)
--------------

The root of the WESTPA HDF5 file contains the following members:

=============== ======================= =======================================
Members         Type                    Description
=============== ======================= =======================================
``iterations/`` Group                   Iteration data
``summary``     Dataset (1-D, compound) Summary data by iteration
=============== ======================= =======================================

Iteration summary table (/summary)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================ ===============================================================
Field            Description
================ ===============================================================
``n_particles``  Total number of walkers in this iteration
``norm``         Total probability, for stability monitoring
``min_bin_prob`` Smallest probability contained in a bin
``max_bin_prob`` Largest probability contained in a bin
``min_seg_prob`` Smallest probability carried by a walker
``max_seg_prob`` Largest probability carried by a walker
``cputime``      Total CPU time (in seconds) spent on propagation for this
                 iteration
``walltime``     Total wallclock time (in seconds) spent on this iteration
================ ===============================================================

Per-iteration data (/iterations/iter_N)
---------------------------------------

Data for each iteration is stored in its own group, named according to the
iteration number and zero-padded out to 8 digits, as in
``/iterations/iter_00000001`` for iteration 1. This is done solely for
convenience in dealing with the data in external utilities that sort output by
group name lexicographically. The field width is in fact configurable via the
``iter_prec`` configuration entry under ``data`` section of the WESTPA
configuration file.

The HDF5 group for each iteration contains the following members:

=================== ======================= ===================================
Member              Type                    Description
=================== ======================= ===================================
``auxdata/``        Group                   User-defined auxiliary data sets
``initial_states``  Dataset (1-D, compound) Initial state of each segment
``final_states``    Dataset (1-D, compound) Final state of each segment
``pcoord``          Dataset (3-D)           Array of shape
                                            ``(n_particles, pcoord_len, pcoord_ndim)``
                                            containing progress coordinate data
``seg_index``       Dataset (1-D, compound) Summary data for each segment
``wtgraph``         Dataset (1-D)           Weight transfer graph data
=================== ======================= ===================================

Segment summary table (/iterations/iter_N/seg_index)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================== ===============================================================
Field              Description
================== ===============================================================
``weight``         Segment weight
``parent_id``      Index of parent
``wtg_n_parents``  Number of entries in ``wtgraph``
``wtg_offset``     Offset into ``wtgraph``
``cputime``        Total CPU time required to propagate the segment
``walltime``       Total walltime required to propagate the segment
``initpoint_type`` Constant indicating the segment's origin
``endpoint_type``  Constant indicating the segment's fate
``status``         Propagation status
================== ===============================================================
