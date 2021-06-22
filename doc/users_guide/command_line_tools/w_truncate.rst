.. _w_truncate:

w_truncate
==========

NOTE: w_truncate only deletes iteration groups from the HDF5 data store. It is recommended that any
iteration data saved to the file system (e.g. in the traj_segs directory) is deleted or moved for the
corresponding iterations.

usage::

 w_truncate [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [-n N_ITER]

Remove all iterations after a certain point in a WESTPA simulation.

optional arguments::

  -h, --help            show this help message and exit
  -n N_ITER, --iter N_ITER
                        Truncate this iteration and those following.

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

