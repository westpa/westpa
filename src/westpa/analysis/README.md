# `westpa.analysis`

This subpackage provides an API to facilitate the analysis of WESTPA
simulation data.
 
## Getting Started

The core abstraction of the `westpa.analysis` package is the `Run` class. A `Run` instance provides a read-only view of a WEST HDF5 ("west.h5") file.
To create a  run, use the `open()` class method:
```python_repl
>>> from westpa.analysis import Run
>>> run = Run.open('west.h5')
>>> run
<WESTPA Run with 500 iterations at 0x7fcaf8f0d5b0>
```

### How To

Iterate over iterations and walkers:
```python
for iteration in run:
    for walker in iteration:
        ...
```

Retrieve a particular walker:
```python
walker = run.iteration(10).walker(4)
```

Get the weight and progress coordinate values of a walker:
```python
weight, pcoords = walker.weight, walker.pcoords
```

Get the parent and children of a walker:
```python
parent, children = walker.parent, walker.children
```

Trace the ancestry of a walker:
```python
trace = walker.trace()
```

### Retrieving Trajectories
 
The `Trajectory` class provides the basic functionality for retrieving 
trajectory data. 
```python
from westpa.analysis import Trajectory
```
Constructing a trajectory requires two ingredients:
1. A function for retrieving individual trajectory segments. The function must
  take a `Walker` object as its first argment and return a sequence (e.g., a 
  list, NumPy array, or 
  [MDTraj Trajectory](https://mdtraj.org/1.9.5/api/generated/mdtraj.Trajectory.html)) 
  representing the trajectory of the walker. Moreover, it
  must accept a Boolean keyword argment `include_initpoint`, which specifies 
  whether the returned trajectory segment includes its initial point.
2. A function for concatenating trajectory segments. The default implementation
  is the `concatenate()` function in the
  `westpa.analysis.trajectories` module. This implementation includes special 
  handling for NumPy arrays and MDTraj Trajectories.

```py
import mdtraj
import os


def get_segment(walker, include_initpoint=True, sim_root='.'):
    """Return the trajectory segment of a walker."""
    seg_dir = os.path.join(
        sim_root,
        'traj_segs',
        format(walker.iteration.number, '06d'),
        format(walker.index, '06d'),
    )
    
    filename = os.path.join(seg_dir, 'seg.dcd')
    top = os.path.join(seg_dir, 'bstate.pdb')
    traj = mdtraj.load(filename, top=top)
    
    if include_initpoint:
        filename = os.path.join(seg_dir, 'parent.xml')
        parent = mdtraj.load(filename, top=traj.top)
        return parent.join(traj, check_topology=False)
    else:
        return traj
```

```python
trajectory = Trajectory(get_segment)
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
To monitor progress, progress bars can be displayed:
```py
trajectory.segment_collector.show_progress = True
```

