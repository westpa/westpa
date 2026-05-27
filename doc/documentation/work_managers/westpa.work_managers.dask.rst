westpa.work\_managers.dask package
==================================

Overview
--------

The Dask work manager is not installed by default and is only avilable on WESTPA v2022.16+. 

Install Dask dependencies with WESTPA with pip using the ``[dask]`` flag (e.g. ``python -m pip install westpa[dask]``). Or you can install the dependencies 
using conda with the following command: ``conda install -c conda-forge dask distributed``.

Usage
-----

Run on the default local cluster::

  $ w_run --work-manager dask

Run on a local cluster with 6 workers::

  $ w_run --work_manager dask --n-workers 6

Start up a cluster and point WESTPA to the scheduler address::

  $ dask scheduler
  INFO -   Scheduler at:     tcp://10.0.0.164:8786
  $ dask worker 10.0.0.164:8786
  $ w_run --work-manager dask --dask-scheduler-address 10.0.0.164:8786


Start up a cluster and point WESTPA to the scheduler file::

  $ dask scheduler --scheduler-file /path/to/scheduler.json
  $ dask worker --scheduler-file /path/to/scheduler.json
  $ w_run --work-manager dask --dask-scheduler-file /path/to/scheduler.json


By default, clients/schedulers/cluster/workers that are passed to WESTPA are considered self-managed and will not be shut down when ``DaskWorkManager`` shuts down. However, all workers will be restarted to ensure a clean state. Users can pass the ``--dask-shutdown-completely`` flag to forcefully shutdown the scheduler/workers. This complete shutdown behavior is default when the dask client/cluster is started by WESTPA.


Users may also modify the dask workers WESTPA manages with the ``--dask-memory-limit <str>`` (equivalent to ``--memory-limit`` for ``dask worker``) and ``--dask-threads-per-worker <int>`` (equivalent to ``--nthreads`` for ``dask worker``) flags. The general ``--n-workers <int>`` flag is equivalent to the ``--nworkers`` flag in dask.


More examples of dask work manager usage can be found in the following location: https://github.com/westpa/westpa_tutorials/tree/main/additional_tutorials.


westpa.work\_managers.dask module
---------------------------------

.. automodule:: westpa.work_managers.dask
   :members:
   :undoc-members:
   :show-inheritance:

