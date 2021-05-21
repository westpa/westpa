# WESTPAnalysis

This package provides an API for analyzing weighted-ensemble simulation data generated using the [WESTPA](http://westpa.github.io/westpa/) framework.

*This package is currently in pre-alpha phase.*
 
## Getting Started

```py
from westpanalysis import WESTDataset
```
A `WESTDataset` provides a view over a WESTPA HDF5 simulation data file.
```py
dataset = WESTDataset('west.h5')
```

### How To

Iterate over simulation iterations and segments:
```py
for iteration in dataset.iterations:
    for segment in iteration.segments:
        ...
```

Retrieve a particular segment:
```py
segment = dataset.iteration(10).segment(4)
```

Get the weight and progress coordinate values of a segment:
```py
weight, pcoords = segment.weight, segment.pcoords
```

Get the parent and children of a segment:
```py
parent, children = segment.parent, segment.children
```

Get the trace (ancestral line) of a segment:
```py
trace = segment.trace()
```

### Loading Segment Trajectories

To enable loading segment trajectories, set the `seg_traj_loader` attribute
of the dataset. An example using [MDTraj](https://www.mdtraj.org) might look like this:
```py
import mdtraj

def load_segment_traj(iteration_number, segment_index):
    return mdtraj.load(
        f'traj_segs/{iteration_number:06d}/{segment_index:06d}.h5')

dataset.segment_traj_loader = load_segment_traj
```

Then the trajectory of a segment can be loaded using
```py
segment_traj = segment.trajectory()
```

and the trajectory of its trace can be loaded using
```py
trace_traj = segment.trace_trajectory()
```

