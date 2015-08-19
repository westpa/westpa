.. _w_pdist:

w_pdist
=======

``w_pcpdist`` constructs and calculates the progress coordinate probability
distribution's evolution over a user-specified number of simulation iterations.
``w_pcpdist`` supports progress coordinates with dimensionality ≥ 1.

The resulting distribution can be viewed with the :ref:`plothist` tool.

Overview
--------

Usage::

  $WEST_ROOT/bin/w_pdist [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                         [-W WEST_H5FILE] [--first-iter N_ITER] [--last-iter N_ITER]
                         [-b BINEXPR] [-o OUTPUT]
                         [--construct-dataset CONSTRUCT_DATASET | --dsspecs DSSPEC [DSSPEC ...]]
                         [--serial | --parallel | --work-manager WORK_MANAGER]
                         [--n-workers N_WORKERS] [--zmq-mode MODE]
                         [--zmq-info INFO_FILE] [--zmq-task-endpoint TASK_ENDPOINT]
                         [--zmq-result-endpoint RESULT_ENDPOINT]
                         [--zmq-announce-endpoint ANNOUNCE_ENDPOINT]
                         [--zmq-listen-endpoint ANNOUNCE_ENDPOINT]
                         [--zmq-heartbeat-interval INTERVAL]
                         [--zmq-task-timeout TIMEOUT] [--zmq-client-comm-mode MODE]

Note: This tool supports parallelization, which may be more efficient for
especially large datasets.

Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for more
information on the general options.

Input/output options
~~~~~~~~~~~~~~~~~~~~

These arguments allow the user to specify where to read input simulation result
data and where to output calculated progress coordinate probability
distribution data.

Both input and output files are *hdf5* format::

  -W, --WEST_H5FILE file
    Read simulation result data from file *file*. (**Default:** The
    *hdf5* file specified in the configuration file (default config file
    is *west.h5*))

  -o, --output file
    Store this tool's output in *file*. (**Default:** The *hdf5* file
    **pcpdist.h5**)

Iteration range options
~~~~~~~~~~~~~~~~~~~~~~~

Specify the range of iterations over which to construct the progress
coordinate probability distribution.::

  --first-iter n_iter
    Construct probability distribution starting with iteration *n_iter*
    (**Default:** 1)

  --last-iter n_iter
    Construct probability distribution's time evolution up to (and
    including) iteration *n_iter* (**Default:** Last completed
    iteration)

Probability distribution binning options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify the number of bins to use when constructing the progress
coordinate probability distribution. If using a multidimensional
progress coordinate, different binning schemes can be used for the
probability distribution for each progress coordinate.::

  -b binexpr
    *binexpr* specifies the number and formatting of the bins. Its
    format can be as follows:

        1. an integer, in which case all distributions have that many
        equal sized bins
        2. a python-style list of integers, of length corresponding to
        the number of dimensions of the progress coordinate, in which
        case each progress coordinate's probability distribution has the
        corresponding number of bins
        3. a python-style list of lists of scalars, where the list at
        each index corresponds to each dimension of the progress
        coordinate and specifies specific bin boundaries for that
        progress coordinate's probability distribution.

    (**Default:** 100 bins for all progress coordinates)

Examples
--------

Assuming simulation results are stored in *west.h5* (which is specified in the
configuration file named *west.cfg*), for a simulation with a 1-dimensional
progress coordinate:

Calculate a probability distribution histogram using all default options
(output file: *pdist.h5*; histogram binning: 100 equal sized bins; probability
distribution over the lowest reached progress coordinate to the largest; work
is parallelized over all available local cores using the 'processes' work
manager)::

  $WEST_ROOT/bin/w_pdist

Same as above, except using the serial work manager (which may be more
efficient for smaller datasets)::

  $WEST_ROOT/bin/w_pdist --serial
