# `westpa.analysis`

This module provides an API for analyzing weighted-ensemble simulation data
generated using the WESTPA framework.
 
## Getting Started

```py
from westpa.analysis import Run
```
A `Run` represents a completed WESTPA simulation. 
Runs may be constructed from WESTPA HDF5 simulation data files.
```py
run = Run('west.h5')
```

### How To

Iterate over iterations and walkers:
```py
for iteration in run:
    for walker in iteration:
        ...
```

Retrieve a particular walker:
```py
walker = run.iteration(10).walker(4)
```

Get the weight and progress coordinate values of a walker:
```py
weight, pcoords = walker.weight, walker.pcoords
```

Get the parent and children of a walker:
```py
parent, children = walker.parent, walker.children
```

Trace the ancestry of a walker back to an initial state:
```py
trace = walker.trace()
```

### Retrieving Trajectories

The `@Trajectory` decorator transforms a function for retrieving the
trajectory of single walker into a callable that can be used to retrieve 
both walker and trace trajectories. The following code, which uses the
[MDTraj](https://www.mdtraj.org/1.9.5/index.html) library, demonstrates its
use:

```py
import mdtraj as md
from westpa.analysis import Trajectory

@Trajectory
def trajectory(walker, atom_indices=None):
    filename = f'traj_segs/{walker.iteration.number:06d}/{walker.index:06d}/seg.h5'
    return md.load(filename, atom_indices=atom_indices)
```
The decorated function can be used to get the trajectory of either a walker:
```py
seg = trajectory(walker)
```
or a trace:
```py
traj = trajectory(walker.trace())
```
To (potentially) speed things up, the trajectory segments of a trace
can be retrieved asynchronously using a pool of threads. This behavior can 
be controlled by setting the `use_threads` and `max_workers` properties
of the trajectory's `segment_collector` attribute. For example:
```py
trajectory.segment_collector.use_threads = True
trajectory.segment_collector.max_workers = 8
```
Even when using multiple threads, loading numerous trajectory segments from 
disk (or worse, downloading them from a remote host) can take time. 
To monitor progress, a progress bar can be displayed:
```py
trajectory.segment_collector.show_progress = True
```

