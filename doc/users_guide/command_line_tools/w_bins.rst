.. _w_bins:

w_bins
======

``w_bins`` deals with binning modification and statistics

Overview
--------

Usage::

  $WEST_ROOT/bin/w_bins [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
               [-W WEST_H5FILE]
               {info,rebin} ...

Display information and statistics about binning in a WEST simulation, or
modify the binning for the current iteration of a WEST simulation.

Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for
more information on the general options.

Options Under 'info'
~~~~~~~~~~~~~~~~~~~~

Usage::

  $WEST_ROOT/bin/w_bins info [-h] [-n N_ITER] [--detail]
                    [--bins-from-system | --bins-from-expr BINS_FROM_EXPR | --bins-from-function BINS_FROM_FUNCTION | --bins-from-file]

Positional options::

  info
    Display information about binning.

Options for 'info'::

  -n N_ITER, --n-iter N_ITER
    Consider initial points of segment N_ITER (default: current
    iteration).

  --detail
    Display detailed per-bin information in addition to summary
    information.

Binning options for 'info'::

  --bins-from-system
    Bins are constructed by the system driver specified in the WEST
    configuration file (default where stored bin definitions not
    available).

  --bins-from-expr BINS_FROM_EXPR, --binbounds BINS_FROM_EXPR
    Construct bins on a rectilinear grid according to the given BINEXPR.
    This must be a list of lists of bin boundaries (one list of bin
    boundaries for each dimension of the progress coordinate), formatted
    as a Python expression. E.g. "[[0,1,2,4,inf],[-inf,0,inf]]". The
    numpy module and the special symbol "inf" (for floating-point
    infinity) are available for use within BINEXPR.

  --bins-from-function BINS_FROM_FUNCTION, --binfunc BINS_FROM_FUNCTION
    Supply an external function which, when called, returns a properly
    constructed bin mapper which will then be used for bin assignments.
    This should be formatted as "[PATH:]MODULE.FUNC", where the function
    FUNC in module MODULE will be used; the optional PATH will be
    prepended to the module search path when loading MODULE.

  --bins-from-file
    Load bin specification from the data file being examined (default
    where stored bin definitions available).

Options Under 'rebin'
~~~~~~~~~~~~~~~~~~~~~

Usage::

  $WEST_ROOT/bin/w_bins rebin [-h] [--confirm] [--detail]
                     [--bins-from-system | --bins-from-expr BINS_FROM_EXPR | --bins-from-function BINS_FROM_FUNCTION]
                     [--target-counts TARGET_COUNTS | --target-counts-from FILENAME]

Positional option::

  rebin
    Rebuild current iteration with new binning.

Options for 'rebin'::

  --confirm
    Commit the revised iteration to HDF5; without this option, the
    effects of the new binning are only calculated and printed.

  --detail
    Display detailed per-bin information in addition to summary
    information.

Binning options for 'rebin';

Same as the binning options for 'info'.

Bin target count options for 'rebin';::

  --target-counts TARGET_COUNTS
    Use TARGET_COUNTS instead of stored or system driver target counts.
    TARGET_COUNTS is a comma-separated list of integers. As a special
    case, a single integer is acceptable, in which case the same target
    count is used for all bins.

  --target-counts-from FILENAME
    Read target counts from the text file FILENAME instead of using
    stored or system driver target counts. FILENAME must contain a list
    of integers, separated by arbitrary whitespace (including newlines).

Input Options
-------------

::

  -W WEST_H5FILE, --west_data WEST_H5FILE
    Take WEST data from WEST_H5FILE (default: read from the HDF5 file
    specified in west.cfg).

Examples
--------

(TODO: Write up an example)
