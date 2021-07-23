.. _plothist_instant:

plothist_instant
================

usage::

 plothist instant [-h] [-o PLOT_OUTPUT] [--hdf5-output HDF5_OUTPUT] [--plot-contour]
                        [--title TITLE] [--linear | --energy | --zero-energy E | --log10]
                        [--range RANGE] [--postprocess-function POSTPROCESS_FUNCTION]
                        [--text-output TEXT_OUTPUT] [--iter N_ITER]
                        input [DIMENSION] [ADDTLDIM]

Plot a probability distribution for a single WE iteration. The probability
distribution must have been previously extracted with ``w_pdist`` (or, at
least, must be compatible with the output format of ``w_pdist``; see
``w_pdist --help`` for more information).

optional arguments::

  -h, --help            show this help message and exit

input options::

  input                 HDF5 file containing histogram data
  DIMENSION             Plot for the given DIMENSION, specified as INT[:[LB,UB]:LABEL], where INT is a
                        zero-based integer identifying the dimension in the histogram, LB and UB are
                        lower and upper bounds for plotting, and LABEL is the label for the plot axis.
                        (Default: dimension 0, full range.)
  ADDTLDIM              For instantaneous/average plots, plot along the given additional dimension,
                        producing a color map.
  --iter N_ITER         Plot distribution for iteration N_ITER (default: last completed iteration).

output options::

  -o PLOT_OUTPUT, --output PLOT_OUTPUT, --plot-output PLOT_OUTPUT
                        Store plot as PLOT_OUTPUT. This may be set to an empty string (e.g. --plot-
                        output='') to suppress plotting entirely. The output format is determined by
                        filename extension (and thus defaults to PDF). Default: "hist.pdf".
  --hdf5-output HDF5_OUTPUT
                        Store plot data in the HDF5 file HDF5_OUTPUT.
  --plot-contour        Determines whether or not to superimpose a contour plot over the heatmap for 2D
                        objects.
  --text-output TEXT_OUTPUT
                        Store plot data in a text format at TEXT_OUTPUT. This option is only valid for
                        1-D histograms. (Default: no text output.)

plot options::

  --title TITLE         Include TITLE as the top-of-graph title
  --linear              Plot the histogram on a linear scale.
  --energy              Plot the histogram on an inverted natural log scale, corresponding to (free)
                        energy (default).
  --zero-energy E       Set the zero of energy to E, which may be a scalar, "min" or "max"
  --log10               Plot the histogram on a base-10 log scale.
  --range RANGE         Plot histogram ordinates over the given RANGE, specified as "LB,UB", where LB
                        and UB are the lower and upper bounds, respectively. For 1-D plots, this is the
                        Y axis. For 2-D plots, this is the colorbar axis. (Default: full range.)
  --postprocess-function POSTPROCESS_FUNCTION
                        Names a function (as in module.function) that will be called just prior to
                        saving the plot. The function will be called as ``postprocess(hist, midpoints,
                        binbounds)`` where ``hist`` is the histogram that was plotted, ``midpoints`` is
                        the bin midpoints for each dimension, and ``binbounds`` is the bin boundaries
                        for each dimension for 2-D plots, or None otherwise. The plot must be modified
                        in place using the pyplot stateful interface.
