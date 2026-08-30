w_mdacrawl
==========

The ``w_mdacrawl`` tool integrates MDAnalysis with WESTPA. It allows you to run any standard or custom MDAnalysis tool (like RMSD, RMSF, or RDF) directly on a WESTPA HDF5 file.

usage::

  w_mdacrawl [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version] [--column COLUMN] [--save DATASET_NAME] [--select SELECTION] [-j N] [--overwrite] west_h5 analysis

positional arguments::

  west_h5               Path to the WESTPA west.h5 file
  analysis              Dotted path to the MDAnalysis analysis class (e.g., MDAnalysis.analysis.rms.RMSD) or a user-defined module.ClassName

options::

  -h, --help            show this help message and exit
  --column COLUMN       Extract a specific column index from 2D results (e.g., 2 for RMSD values)
  --save DATASET_NAME   Save the results back into west.h5 under iterations/iter_XXXXXXXX/auxdata/DATASET_NAME
  --select SELECTION    MDAnalysis atom selection string (e.g., "name CA", "protein")
  -j N, --n-workers N   Number of parallel workers for multiprocessing backend (default: 1 -> serial)
  --overwrite           Overwrite the dataset if it already exists in the auxdata dataset

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit


westpa.cli.tools.w\_mdacrawl module
-----------------------------------

.. automodule:: westpa.cli.tools.w_mdacrawl
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members: