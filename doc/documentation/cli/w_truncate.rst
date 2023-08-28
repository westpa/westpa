w_truncate
==========

``w_truncate`` removes all iterations after a certain point

Overview
--------

Usage::

  w_truncate [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                  [-n N_ITER] [-W WEST_H5FILE]

Remove all iterations after a certain point in a

Command-Line Options
--------------------

See the `command-line tool index <command_line_tool_index>` for more
information on the general options.

Iteration Options
~~~~~~~~~~~~~~~~~

::

  -n N_ITER, --iter N_ITER
    Truncate this iteration and those following.

  -W WEST_H5FILE, --west-data WEST_H5FILE
    PATH of H5 file to truncate. By default, it will read from the RCFILE (e.g., west.cfg).
    This option will have override whatever's provided in the RCFILE.

Examples
--------

Running the following will remove iteration 50 and all iterations after 50 from multi.h5.

::

    w_truncate -n 50 -W multi.h5


westpa.cli.core.w\_truncate module
----------------------------------

.. automodule:: westpa.cli.core.w_truncate
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
