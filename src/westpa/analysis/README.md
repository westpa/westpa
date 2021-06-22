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

Iterate over simulation iterations and segments:
```py
for iteration in run:
    for segment in iteration:
        ...
```

Retrieve a particular segment:
```py
segment = run.iteration(10).segment(4)
```

Get the weight and progress coordinate values of a segment:
```py
weight, pcoords = segment.weight, segment.pcoords
```

Get the parent and children of a segment:
```py
parent, children = segment.parent, segment.children
```

Trace the ancestry of a segment back to its origin:
```py
trace = segment.trace()
```

### Loading Segment Trajectories

```py
import mdtraj

from westpa.analysis import segment_trajectory

@segment_trajectory
def traj(segment):
    filename = f'traj_segs/{segment.iteration.number}/{segment.index}/seg.h5'
    return mdtraj.load(filename)
```

Then the trajectory of a segment can be accessed via the `traj`
property:
```py
segment_traj = segment.traj
```
The trajectory of a trace can be accessed similarly:
```py
trace = segment.trace()
trace_traj = trace.traj
```
