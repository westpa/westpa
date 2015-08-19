.. _plothist:

plothist
========

Use the ``plothist`` tool to plot the results of :ref:`w_pdist`. This tool uses
an *hdf5* file as its input (i.e. the output of another analysis tool), and
outputs a *pdf* image.

The ``plothist`` tool operates in one of three (mutually exclusive)
plotting modes:

-  **``evolution``**: Plots the relevant data as a time evolution over
   specified number of simulation iterations
-  **``average``**: Plots the relevant data as a time average over a
   specified number of iterations
-  **``instant``**: Plots the relevant data for a single specified
   iteration

Overview
--------

The basic usage, independent of plotting mode, is as follows:

usage:

| ``$WEST_ROOT/bin/plothist [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]``
| ``                        {instant,average,evolution} input ...``

Note that the user must specify a plotting mode (i.e. '``instant``\ ',
'``average``\ ', or '``evolution``\ ') and an input file, ``input``.

Therefore, this tool is always called as:

``$WEST_ROOT/bin/plothist mode input_file [``\ *``other``
``options``*\ ``]``

'``instant``\ ' mode
~~~~~~~~~~~~~~~~~~~~

usage:

| ``$WEST_ROOT/bin/plothist instant [-h] input [-o PLOT_OUTPUT]``
| ``                                [--hdf5-output HDF5_OUTPUT] [--text-output TEXT_OUTPUT]``
| ``                                [--title TITLE] [--range RANGE] [--linear | --energy | --log10]``
| ``                                [--iter N_ITER] ``
| ``                                [DIMENSION] [ADDTLDIM]``

'``average``\ ' mode
~~~~~~~~~~~~~~~~~~~~

usage:

| ``$WEST_ROOT/bin/plothist average [-h] input [-o PLOT_OUTPUT]``
| ``                                [--hdf5-output HDF5_OUTPUT] [--text-output TEXT_OUTPUT]``
| ``                                [--title TITLE] [--range RANGE] [--linear | --energy | --log10]``
| ``                                [--first-iter N_ITER] [--last-iter N_ITER]                           ``
| ``                                [DIMENSION] [ADDTLDIM]``

'``evolution``\ ' mode
~~~~~~~~~~~~~~~~~~~~~~

usage:

| ``$WEST_ROOT/bin/plothist evolution [-h] input [-o PLOT_OUTPUT]``
| ``                                  [--hdf5-output HDF5_OUTPUT]``
| ``                                  [--title TITLE] [--range RANGE] [--linear | --energy | --log10]``
| ``                                  [--first-iter N_ITER] [--last-iter N_ITER]``
| ``                                  [--step-iter STEP]                                   ``
| ``                                  [DIMENSION]``

Command-Line Options
--------------------

See the :ref:`command-line tool index <command_line_tool_index>` for more
information on the general options.

Unless specified (as a **Note** in the command-line option description), the
command-line options below are shared for all three plotting modes

Input/output options
~~~~~~~~~~~~~~~~~~~~

No matter the mode, an input *hdf5* file must be specified. There are
three possible outputs that are mode or user-specified: A text file, an
*hdf5* file, and a pdf image.

Specifying input file
^^^^^^^^^^^^^^^^^^^^^

***``input``***
    Specify the input *hdf5* file ''``input``. This is the output file
    from a previous analysis tool (e.g. 'pcpdist.h5')

Output plot pdf file
^^^^^^^^^^^^^^^^^^^^

**``-o ''plot_output'', --plot_output ''plot_output''``**
    Specify the name of the pdf plot image output (**Default:**
    'hist.pdf').
    **Note:** You can suppress plotting entirely by specifying an empty string
    as *plot\_output* (i.e. ``-o ''`` or ``--plot_output ''``)

Additional output options
^^^^^^^^^^^^^^^^^^^^^^^^^

Note: ``plothist`` provides additional, optional arguments to output the
data points used to construct the plot:

**``--hdf5-output ''hdf5_output''``**
    Output plot data *hdf5* file ``'hdf5_output'`` (**Default:** No
    *hdf5* output file)

**``--text-output ''text_output''``**
    Output plot data as a text file named ``'text_output'``
    (**Default:** No text output file)
    **Note:** This option is only available for 1 dimensional histogram
    plots (that is, ``'average'`` and ``'instant'`` modes only)

Plotting options
~~~~~~~~~~~~~~~~

The following options allow the user to specify a plot title, the type
of plot (i.e. energy or probability distribution), whether to apply a
log transformation to the data, and the range of data values to include.

**``--title ''title'' ``**
    Optionally specify a title, *``title``*, for the plot (**Default:**
    No title)

**``--range ''<nowiki>'</nowiki>LB, UB<nowiki>'</nowiki>''``**
    Optionally specify the data range to be plotted as "``LB, UB``\ "
    (e.g. ``' --range "-1, 10" '`` - note that the quotation marks are
    necessary if specifying a negative bound). For 1 dimensional
    histograms, the range affects the y axis. For 2 dimensional plots
    (e.g. evolution plot with 1 dimensional progress coordinate), it
    corresponds to the range of the color bar

Mutually exclusive plotting options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following three options determine how the plotted data is
represented (**Default:** ``'--energy'``)

**``--energy ``**
    Plots the probability distribution on an inverted natural log scale
    (i.e. -ln[P(x)] ), corresponding to the free energy (**Default**)

**``--linear ``**
    Plots the probability distribution function as a linear scale

**``--log10 ``**
    Plots the (base-10) logarithm of the probability distribution

Iteration selection options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on plotting mode, you can select either a range or a single
iteration to plot.

**``'instant'``** mode only:

**``--iter ''n_iter'' ``**
    Plot the distribution for iteration ``''n_iter''`` (**Default:**
    Last completed iteration)

**``'average'``** and **``'evolution'``** modes only:

**``--first-iter ''first_iter'' ``**
    Begin averaging or plotting at iteration *``first_iter``*
    (**Default:** 1)

**``--last-iter ''last_iter'' ``**
    Average or plot up to and including *``last_iter``* (**Default:**
    Last completed iteration)

**``'evolution'``** mode only:

**``--iter_step ''n_step'' ``**
    Average every *``n_step``* iterations together when plotting in
    ``'evolution'`` mode (**Default:** 1 - i.e. plot each iteration)

Specifying progress coordinate dimension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For progress coordinates with dimensions greater than 1, you can specify
the dimension of the progress coordinate to use, the of progress
coordinate values to include, and the progress coordinate axis label
with a single positional argument:

**``dimension ``**
    Specify ``'dimension'`` as '``int[:[LB,UB]:label]``\ ', where
    '``int``\ ' specifies the dimension (starting at 0), and,
    optionally, '``LB,UB``\ ' specifies the lower and upper range
    bounds, and/or '``label``\ ' specifies the axis label (**Default:**
    ``int`` = 0, full range, default label is 'dimension ``int``'; e.g
    'dimension 0')

For ``'average'`` and ``'instant'`` modes, you can plot two dimensions
at once using a color map if this positional argument is specified:

**``addtl_dimension ``**
    Specify the other dimension to include as ``'addtl_dimension'``

Examples
--------

These examples assume the input file is created using w\_pcpdist and is
named 'pcpdist.h5'

Basic plotting
~~~~~~~~~~~~~~

Plot the energy ( -ln(P(x)) ) for the last iteration

``$WEST_ROOT/bin/plothist instant pcpdist.h5``

Plot the evolution of the log10 of the probability distribution over all
iterations

``$WEST_ROOT/bin/plothist evolution pcpdist.h5 --log10 ``

Plot the average linear probability distribution over all iterations

``$WEST_ROOT/bin/plothist average pcpdist.h5 --linear``

Specifying progress coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plot the average probability distribution as the energy, label the
x-axis 'pcoord', over the entire range of the progress coordinate

``$WEST_ROOT/bin/plothist average pcpdist.h5 0::pcoord``

Same as above, but only plot the energies for with progress coordinate
between 0 and 10

``$WEST_ROOT/bin/plothist average pcpdist.h5 '0:0,10:pcoord'``

(Note: the quotes are needed if specifying a range that includes a
negative bound)

(For a simulation that uses at least 2 progress coordinates) plot the
probability distribution for the 5th iteration, representing the first
two progress coordinates as a heatmap

``$WEST_ROOT/bin/plothist instant pcpdist.h5 0 1 --iter 5 --linear``
