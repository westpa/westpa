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

### Retrieving Trajectory Segments

The `trajectory_segment()` decorator transforms a function for loading the
trajectory of a particular walker into a trajectory descriptor attached to
both `Walker` and `Trace` instances. The following code (which uses the
[MDTraj](https://www.mdtraj.org/1.9.5/index.html) library) demonstrates its
use:

```py
import mdtraj

from westpa.analysis import trajectory_segment

@trajectory_segment
def traj(walker):
    filename = f'traj_segs/{walker.iteration.number}/{walker.index}/seg.h5'
    return mdtraj.load(filename)
```

The trajectory segment associated with a given walker can then be accessed via the `traj` attribute:
```py
segment = walker.traj
```
The trajectory of a trace (i.e., the concatentation of a sequence of 
trajectory segments) can be accessed similarly:
```py
trace = walker.trace()
trace_traj = trace.traj
```

