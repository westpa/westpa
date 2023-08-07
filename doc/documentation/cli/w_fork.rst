w_fork
======

usage::

 w_fork [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [-i INPUT_H5FILE]
              [-I N_ITER] [-o OUTPUT_H5FILE] [--istate-map ISTATE_MAP] [--no-headers]

Prepare a new weighted ensemble simulation from an existing one at a particular point. A new HDF5 file
is generated. In the case of executable propagation, it is the user's responsibility to prepare the new
simulation directory appropriately, particularly making the old simulation's restart data from the
appropriate iteration available as the new simulations initial state data; a mapping of old simulation
segment to new simulation initial states is created, both in the new HDF5 file and as a flat text file,
to aid in this. Target states and basis states for the new simulation are taken from those in the
original simulation.

optional arguments::

  -h, --help            show this help message and exit
  -i INPUT_H5FILE, --input INPUT_H5FILE
                        Create simulation from the given INPUT_H5FILE (default: read from configuration
                        file.
  -I N_ITER, --iteration N_ITER
                        Take initial distribution for new simulation from iteration N_ITER (default:
                        last complete iteration).
  -o OUTPUT_H5FILE, --output OUTPUT_H5FILE
                        Save new simulation HDF5 file as OUTPUT (default: forked.h5).
  --istate-map ISTATE_MAP
                        Write text file describing mapping of existing segments to new initial states
                        in ISTATE_MAP (default: istate_map.txt).
  --no-headers          Do not write header to ISTATE_MAP

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit


westpa.cli.tools.w\_fork module
-------------------------------

.. automodule:: westpa.cli.core.w_fork
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
