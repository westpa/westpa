.. _w_fluxanl:

w_fluxanl
=========

``w_fluxanl`` calculates the probability flux of a weighted ensemble simulation
based on a pre-defined target state. Also calculates confidence interval of
average flux. Monte Carlo bootstrapping techniques are used to account for
autocorrelation between fluxes and/or errors that are not normally distributed.

Overview
--------

usage::

  $WEST_ROOT/bin/w_fluxanl [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                           [-W WEST_H5FILE] [-o OUTPUT]
                           [--first-iter N_ITER] [--last-iter N_ITER]
                           [-a ALPHA] [--autocorrel-alpha ACALPHA] [-N NSETS] [--evol] [--evol-step ESTEP]

Note: All command line arguments are optional for ``w_fluxanl``.

Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for more
information on the general options.

Input/output options
~~~~~~~~~~~~~~~~~~~~

These arguments allow the user to specify where to read input simulation result
data and where to output calculated progress coordinate probability
distribution data.

Both input and output files are *hdf5* format.::

  -W, --west-data file
    Read simulation result data from file *file*. (**Default:** The
    *hdf5* file specified in the configuration file)

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

Confidence interval and bootstrapping options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify alpha values of constructed confidence intervals.::

  -a alpha
    Calculate a (1 - *alpha*) confidence interval for the mean flux
    (**Default:** 0.05)

  --autocorrel-alpha ACalpha
    Identify autocorrelation of fluxes at *ACalpha* significance level.
    Note: Specifying an *ACalpha* level that is too small may result in
    failure to find autocorrelation in noisy flux signals (**Default:**
    Same level as *alpha*)

  -N n_sets, --nsets n_sets
    Use *n_sets* samples for bootstrapping (**Default:** Chosen based
    on *alpha*)

  --evol
    Calculate the time evolution of flux confidence intervals
    (**Warning:** computationally expensive calculation)

  --evol-step estep
    (if ``'--evol'`` specified) Calculate the time evolution of flux
    confidence intervals for every *estep* iterations (**Default:** 1)

Examples
--------

Calculate the time evolution flux every 5 iterations::

  $WEST_ROOT/bin/w_fluxanl --evol --evol-step 5

Calculate mean flux confidence intervals at 0.01 signicance level and
calculate autocorrelations at 0.05 significance::

  $WEST_ROOT/bin/w_fluxanl --alpha 0.01 --autocorrel-alpha 0.05

Calculate the mean flux confidence intervals using a custom bootstrap
sample size of 500::

  $WEST_ROOT/bin/w_fluxanl --n-sets 500
