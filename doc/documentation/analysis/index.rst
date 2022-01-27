westpa.analysis package
=======================

This subpackage provides an API to facilitate the analysis of WESTPA
simulation data. Its core abstraction is the ``Run`` class.
A ``Run`` instance provides a read-only view of a WEST HDF5 ("west.h5") file.

**API reference:** `<https://westpa.readthedocs.io/en/latest/documentation/analysis/>`_

How To
------

Open a run::

    >>> from westpa.analysis import Run
    >>> run = Run.open('west.h5')
    >>> run
    <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>

Iterate over iterations and walkers::

    >>> for iteration in run:
    ...     for walker in iteration:
    ...         pass
    ...


Access a particular iteration::

    >>> iteration = run.iteration(10)
    >>> iteration
    Iteration(10, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))

Access a particular walker::

    >>> walker = iteration.walker(4)
    >>> walker
    Walker(4, Iteration(10, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))


Get the weight and progress coordinate values of a walker::

    >>> walker.weight
    9.876543209876543e-06
    >>> walker.pcoords
    array([[3.1283207],
           [3.073721 ],
           [2.959221 ],
           [2.6756208],
           [2.7888207]], dtype=float32)


Get the parent and children of a walker::

    >>> walker.parent
    Walker(2, Iteration(9, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    >>> for child in walker.children:
    ...     print(child)
    ...
    Walker(0, Iteration(11, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(1, Iteration(11, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(2, Iteration(11, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(3, Iteration(11, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(4, Iteration(11, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))

Trace the ancestry of a walker::

    >>> trace = walker.trace()
    >>> trace
    Trace(Walker(4, Iteration(10, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>)))
    >>> for walker in trace:
    ...     print(walker)
    ...
    Walker(1, Iteration(1, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(4, Iteration(2, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(5, Iteration(3, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(6, Iteration(4, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(9, Iteration(5, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(8, Iteration(6, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(8, Iteration(7, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(13, Iteration(8, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(2, Iteration(9, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))
    Walker(4, Iteration(10, <WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>))

Close a run (and its underlying HDF5 file)::

    >>> run.close()
    >>> run
    <Closed WESTPA Run at 0x7fcaf8f0d5b0>
    >>> run.h5file
    <Closed HDF5 file>


Retrieving Trajectories
-----------------------

Built-in Reader
^^^^^^^^^^^^^^^

MD trajectory data stored in an identical manner as in the
`Basic NaCl tutorial <https://github.com/westpa/westpa_tutorials/tree/main/basic_nacl>`_
may be retrieved using the built-in ``BasicMDTrajectory`` reader with its
default settings::

    >>> from westpa.analysis import BasicMDTrajectory
    >>> trajectory = BasicMDTrajectory()

Here ``trajectory`` is a callable object that takes either a ``Walker`` or
a ``Trace`` instance as input and returns an
`MDTraj Trajectory <https://mdtraj.org/1.9.5/api/generated/mdtraj.Trajectory.html>`_::

    >>> traj = trajectory(walker)
    >>> traj
    <mdtraj.Trajectory with 5 frames, 33001 atoms, 6625 residues, and unitcells at 0x7fcae484ad00>
    >>> traj = trajectory(trace)
    >>> traj
    <mdtraj.Trajectory with 41 frames, 33001 atoms, 6625 residues, and unitcells at 0x7fcae487c790>

Minor variations of the "basic" trajectory storage protocol (e.g., use of
different file formats) can be handled by changing the parameters of the
``BasicMDTrajectory`` reader. For example, suppose that instead of storing
the coordinate and topology data for trajectory segments in separate
files ("seg.dcd" and "bstate.pdb"), we store them together in a
`MDTraj HDF5 <https://mdtraj.org/1.9.5/hdf5_format.html>`_ trajectory file
("seg.h5"). This change can be accommodated by explicitly setting the
``traj_ext`` and ``top`` parameters of the trajectory reader::

    >>> trajectory = BasicMDTrajectory(traj_ext='.h5', top=None)

Trajectories that are saved with the HDF5 Framework can use ``HDF5MDTrajectory`` reader instead.


Custom Readers
^^^^^^^^^^^^^^

For users requiring greater flexibility, custom trajectory readers can be
implemented using the ``westpa.analysis.Trajectory`` class. Implementing
a custom reader requires two ingredients:

#. A function for retrieving individual trajectory segments. The function
   must take a ``Walker`` instance as its first argument and return a sequence
   (e.g., a list, NumPy array, or MDTraj Trajectory) representing the
   trajectory of the walker. Moreover, it must accept a Boolean keyword
   argument ``include_initpoint``, which specifies whether the returned
   trajectory includes its initial point.
#. A function for concatenating trajectory segments. A default implementation
   is provided by the ``concatenate()`` function in the
   ``westpa.analysis.trajectories`` module.

westpa.analysis.core module
---------------------------

.. automodule:: westpa.analysis.core
   :members:

westpa.analysis.trajectories module
-----------------------------------

.. automodule:: westpa.analysis.trajectories
   :members:

westpa.analysis.statistics module
---------------------------------

.. automodule:: westpa.analysis.statistics
   :members:
