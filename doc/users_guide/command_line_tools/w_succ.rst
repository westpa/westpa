.. _w_succ:

w_succ
======

usage::

 w_succ [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [-A H5FILE] [-W WEST_H5FILE]
              [-o OUTPUT_FILE]

List segments which successfully reach a target state.

optional arguments::

  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Store output in OUTPUT_FILE (default: write to standard output).

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

general analysis options::

  -A H5FILE, --analysis-file H5FILE
                        Store intermediate and final results in H5FILE (default: analysis.h5).

WEST input data options::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).