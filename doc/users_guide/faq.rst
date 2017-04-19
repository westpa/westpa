Frequently asked questions (FAQ)
=================================

Simulation
-----------

- How can I cleanly shutdown a simulation (without corrupting the h5 
  file)? 

It is generally safe to shutdown a WESTPA simulation by simply canceling
the job through your queue management. However, to ensure data integrity
in the h5 file, you should wait until the WESTPA log indicates that an
iteration has begun or is occurring; canceling a job too quickly after
submission can result in the absolute corruption of the h5 file and
should be avoided.

- Storage of Large Files

During a normal WESTPA run, many small files are created and it is
convenient to tar these into a larger file (one tarball per iteration,
for instance). It is generally best to do this 'offline'. An important
aspect to consider is that some disk systems, such as LUSTRE, will
suffer impaired performance if very large files are created. On
Stampede, for instance, any file larger than 200 GB must be 'striped'
properly (such that its individual bits are spread across numerous
disks).

Within the user guide for such systems, there is generally a section on
how to handle large files. Some computers have special versions of tar
which stripe appropriately; others do not (such as Stampede). For those
that do not, it may be necessary to contact the sysadmin, and/or create
a directory where you can place your tarball with a different stripe
level than the default.

- H5py Inflate() Failed error

While running or analyzing a simulation, you may run into an error such
as ``IOError: Can't write data (Inflate() failed)``. These errors may be
related to an open bug in H5py. However, the following tips may help you
to find a workaround.

WESTPA may present you with such an error when unable to read or write a
data set. In the case that a simulation gives this error when you
attempt to run it, it may be helpful to check if a data set may be read
or written to using an interactive Python session. Restarting the
simulation may require deleting and remaking the data set. Also, this
error may be related to compression and other storage options. Thus, it
may be helpful to disable compression and chunked storage. Note that
existing datasets will retain compression and other options given to
them at the time of their creation, so it may be necessary to truncate
an iteration (for example, using ``w_truncate``) in order for changes to
take effect.

This error may also occur during repeated opening (e.g., 1000s of times)
of an HDF5 data set. Thus, this error may occur while running analysis
scripts. In this case, it may be helpful to cache data sets in physical
memory (RAM) as ``numpy`` arrays when they are read, so that the script
loads the dataset a minimal number of times.

- Dynamics Packages

WESTPA was designed to work cleanly with any dynamics package available
(using the executable propagator); however, many of the tips and tricks
available on the web or the user manual for these packages make the
(reasonable) assumption that you will be running a set of brute force
trajectories. As such, some of their guidelines for handling periodic
boundary conditions may not be applicable.

- How can I restart a WESTPA simulation?

In general restarting a westpa simulation will restart an incomplete 
iteration, retaining data from segments that have completed and 
re-running segments that were incomplete (or never started).

In case that the iteration data got corrupted or you want to go
back to an specific iteration and change something, you need to 
delete all the trajectory secgments and other files related to that 
iteration and run w_truncate on that iteration. This will delete westpa's 
information about the nth iteration, which includes which segments have 
run and which have not. Then restarting your westpa simulation will 
restart that iteration afresh.


GROMACS
--------

- Periodic Boundary Conditions

While many of the built in tools now handle periodic boundary conditions
cleanly (such as g\_dist) with relatively little user interaction,
others, such as g\_rms, do not. If your simulation analysis protocol
requires you to run such a tool, you must correct for the periodic
boundary conditions before running it. While there are guidelines
available to help you correct for whatever conditions your system may
have
`here <http://www.gromacs.org/Documentation/Terminology/Periodic_Boundary_Conditions>`__,
there is an implicit assumption that you have one long running
trajectory.

It will be necessary, within your executable propagator (usually
runseg.sh) to run trjconv (typically, two or three times, depending on
your needs: once to remove the periodic boundary conditions, then to
make molecules whole, then to remove any jumps). If no extra input is
supplied (the -s flag in GROMACS 4.X), GROMACS uses the first frame of
your segment trajectory as a reference state to remove jumps. If your
segment's parent ended the previous iteration having jumped across the
box barrier, trjconv will erroneously assume this is the correct state
and 'correct' any jump back across the barrier. **This can result in
unusually high RMSD values for one segment for one or more iterations,**
and can show as discontinuities on the probability distribution. It is
important to note that a lack of discontinuities does not imply a lack
of imaging problems.

To fix this, simply pass in the last frame of the imaged parent
trajectory and use that as the reference structure for trjconv. This
will ensure that trjconv is aware if your segment has crossed the
barrier at time 0 and will make the appropriate corrections.

Development
-------------

- I'm trying to profile a parallel script using the --profile
   option of bin/west. I get a PicklingError. What gives?

When executing a script using --profile, the following error may crop
up:

::

    PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed

The cProfile module used by the --profile option modifies function
definitions such that they are no longer pickleable, meaning that they
cannot be passed through the work manager to other processes. If you
absolutely must profile a parallel script, use the threads work manager.

