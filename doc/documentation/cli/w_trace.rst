w_trace
=======

usage::

    w_trace [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [-W WEST_H5FILE]
               [-d DSNAME] [--output-pattern OUTPUT_PATTERN] [-o OUTPUT]
               N_ITER:SEG_ID [N_ITER:SEG_ID ...]

Trace individual WEST trajectories and emit (or calculate) quantities along the
trajectory.

Trajectories are specified as N_ITER:SEG_ID pairs. Each segment is traced back
to its initial point, and then various quantities (notably n_iter and seg_id)
are printed in order from initial point up until the given segment in the given
iteration.

Output is stored in several files, all named according to the pattern given by
the -o/--output-pattern parameter. The default output pattern is "traj_%d_%d",
where the printf-style format codes are replaced by the iteration number and
segment ID of the terminal segment of the trajectory being traced.

Individual datasets can be selected for writing using the ``-d/--dataset`` option
(which may be specified more than once). The simplest form is ``-d dsname``,
which causes data from dataset ``dsname`` along the trace to be stored to
HDF5.  The dataset is assumed to be stored on a per-iteration basis, with
the first dimension corresponding to seg_id and the second dimension
corresponding to time within the segment.  Further options are specified
as comma-separated key=value pairs after the data set name, as in::

    -d dsname,alias=newname,index=idsname,file=otherfile.h5,slice=[100,...]

The following options for datasets are supported::

    alias=newname
        When writing this data to HDF5 or text files, use ``newname``
        instead of ``dsname`` to identify the dataset. This is mostly of
        use in conjunction with the ``slice`` option in order, e.g., to
        retrieve two different slices of a dataset and store then with
        different names for future use.

    index=idsname
        The dataset is not stored on a per-iteration basis for all
        segments, but instead is stored as a single dataset whose
        first dimension indexes n_iter/seg_id pairs. The index to
        these n_iter/seg_id pairs is ``idsname``.

    file=otherfile.h5
        Instead of reading data from the main WEST HDF5 file (usually
        ``west.h5``), read data from ``otherfile.h5``.

    slice=[100,...]
        Retrieve only the given slice from the dataset. This can be
        used to pick a subset of interest to minimize I/O.


positional arguments
--------------------

::

  N_ITER:SEG_ID         Trace trajectory ending (or at least alive at) N_ITER:SEG_ID.

optional arguments
------------------

::

  -h, --help            show this help message and exit
  -d DSNAME, --dataset DSNAME
                        Include the dataset named DSNAME in trace output. An extended form like
                        DSNAME[,alias=ALIAS][,index=INDEX][,file=FILE][,slice=SLICE] will obtain the
                        dataset from the given FILE instead of the main WEST HDF5 file, slice it by
                        SLICE, call it ALIAS in output, and/or access per-segment data by a
                        n_iter,seg_id INDEX instead of a seg_id indexed dataset in the group for
                        n_iter.

general options
---------------

::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

WEST input data options
-----------------------

::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

output options
--------------

::

  --output-pattern OUTPUT_PATTERN
                        Write per-trajectory data to output files/HDF5 groups whose names begin with
                        OUTPUT_PATTERN, which must contain two printf-style format flags which will be
                        replaced with the iteration number and segment ID of the terminal segment of
                        the trajectory being traced. (Default: traj_%d_%d.)
  -o OUTPUT, --output OUTPUT
                        Store intermediate data and analysis results to OUTPUT (default: trajs.h5).


westpa.cli.tools.w\_trace module
--------------------------------

.. automodule:: westpa.cli.tools.w_trace
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
