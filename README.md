# WESTPAnalysis

This package provides an API for analyzing weighted-ensemble simulation data generated using the [WESTPA](http://westpa.github.io/westpa/) framework.
 
## Getting Started

```py
from westpanalysis import Run
```
A `Run` represents a completed WESTPA simulation. 
Runs may be "loaded" from WESTPA HDF5 simulation data files.
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

Trace a segment back to its origin:
```py
trace = segment.trace()
```

### Loading Segment Trajectories

To enable loading segment trajectories, set the `segment_traj_loader` attribute
of the `Run`. An example using [MDTraj](https://www.mdtraj.org) might look like this:
```py
import mdtraj

def load_segment_traj(iteration_number, segment_index):
    return mdtraj.load(
        f'traj_segs/{iteration_number:06d}/{segment_index:06d}/seg.h5')

run.segment_traj_loader = load_segment_traj
```

Then the trajectory of a segment can be loaded using
```py
segment_traj = segment.trajectory()
```

and the trajectory of its trace can be loaded using
```py
trace_traj = segment.trace().trajectory()
```

