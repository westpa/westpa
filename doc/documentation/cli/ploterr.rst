ploterr
=======

usage::

 ploterr [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                {help,d.kinetics,d.probs,rw.probs,rw.kinetics,generic} ...

Plots error ranges for weighted ensemble datasets.

Command-line options
--------------------

optional arguments::

  -h, --help            show this help message and exit

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

supported input formats::

  {help,d.kinetics,d.probs,rw.probs,rw.kinetics,generic}
    help                print help for this command or individual subcommands
    d.kinetics          output of w_direct kinetics
    d.probs             output of w_direct probs
    rw.probs            output of w_reweight probs
    rw.kinetics         output of w_reweight kinetics
    generic             arbitrary HDF5 file and dataset

westpa.cli.tools.ploterr module
-------------------------------

.. automodule:: westpa.cli.tools.ploterr
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
