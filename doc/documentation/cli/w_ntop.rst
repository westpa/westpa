w_ntop
======

usage::

 w_ntop [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [-W WEST_H5FILE]
              [--first-iter N_ITER] [--last-iter N_ITER] [-a ASSIGNMENTS] [-n COUNT] [-t TIMEPOINT]
              [--highweight | --lowweight | --random] [-o OUTPUT]

Select walkers from bins . An assignment file mapping walkers to
bins at each timepoint is required (see``w_assign --help`` for further
information on generating this file). By default, high-weight walkers are
selected (hence the name ``w_ntop``: select the N top-weighted walkers from
each bin); however, minimum weight walkers and randomly-selected walkers
may be selected instead.

Output format
-------------

The output file (-o/--output, by default "ntop.h5") contains the following
datasets::

  ``/n_iter`` [iteration]
    *(Integer)* Iteration numbers for each entry in other datasets.

  ``/n_segs`` [iteration][bin]
    *(Integer)* Number of segments in each bin/state in the given iteration.
    This will generally be the same as the number requested with
    ``--n/--count`` but may be smaller if the requested number of walkers
    does not exist.

  ``/seg_ids`` [iteration][bin][segment]
    *(Integer)* Matching segments in each iteration for each bin.
    For an iteration ``n_iter``, only the first ``n_iter`` entries are
    valid. For example, the full list of matching seg_ids in bin 0 in the
    first stored iteration is ``seg_ids[0][0][:n_segs[0]]``.

  ``/weights`` [iteration][bin][segment]
    *(Floating-point)* Weights for each matching segment in ``/seg_ids``.

Command-line arguments
----------------------

optional arguments::

  -h, --help            show this help message and exit
  --highweight          Select COUNT highest-weight walkers from each bin.
  --lowweight           Select COUNT lowest-weight walkers from each bin.
  --random              Select COUNT walkers randomly from each bin.

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

WEST input data options::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

iteration range::

  --first-iter N_ITER   Begin analysis at iteration N_ITER (default: 1).
  --last-iter N_ITER    Conclude analysis with N_ITER, inclusive (default: last completed iteration).

input options::

  -a ASSIGNMENTS, --assignments ASSIGNMENTS
                        Use assignments from the given ASSIGNMENTS file (default: assign.h5).

selection options::

  -n COUNT, --count COUNT
                        Select COUNT walkers from each iteration for each bin (default: 1).
  -t TIMEPOINT, --timepoint TIMEPOINT
                        Base selection on the given TIMEPOINT within each iteration. Default (-1)
                        corresponds to the last timepoint.

output options::

  -o OUTPUT, --output OUTPUT
                        Write output to OUTPUT (default: ntop.h5).


westpa.cli.tools.w\_ntop module
-------------------------------

.. automodule:: westpa.cli.tools.w_ntop
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
