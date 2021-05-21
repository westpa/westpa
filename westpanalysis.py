from westpa.core.h5io import WESTPAH5File


class Segment:

    def __init__(self, index, iteration):
        self.index = index
        self.iteration = iteration

    @property
    def dataset(self):
        return self.iteration.dataset

    @property
    def info(self):
        return self.iteration.segment_info[self.index]

    @property
    def weight(self):
        return self.info['weight']

    @property
    def pcoords(self):
        return self.iteration.pcoord[self.index]

    @property
    def num_snapshots(self):
        return self.pcoords.shape[0]

    def __len__(self):
        return self.num_snapshots

    @property
    def parent(self):
        if not self.iteration.prev:
            return None
        return Segment(self.info['parent_id'], self.iteration.prev)

    @property
    def children(self):
        if not self.iteration.next:
            return []
        return [segment for segment in self.iteration.next.segments
                if segment.parent == self]

    def trace(self):
        return trace(self)

    def trajectory(self):
        if self.dataset.segment_traj_loader:
            return self.dataset.segment_traj_loader(
                self.iteration.number, self.index)
        return None

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class Iteration:

    def __init__(self, number, dataset):
        self.number = number
        self.dataset = dataset

    @property
    def group(self):
        return self.dataset.iteration_group(self.number)

    @property
    def segment_info(self):
        return self.group['seg_index']

    @property
    def pcoord(self):
        return self.group['pcoord']

    @property
    def num_segments(self):
        return self.segment_info.shape[0]

    @property
    def segment_weights(self):
        return self.segment_info['weight']

    @property
    def segments(self):
        return (Segment(index, self) for index in range(self.num_segments))

    def segment(self, index):
        valid_range = range(self.num_segments)
        if index not in valid_range:
            raise ValueError(f'segment index must be in {valid_range}')
        return Segment(index, self)

    @property
    def prev(self):
        if self.number == 1:
            return None
        return self.dataset.iteration(self.number - 1)

    @property
    def next(self):
        if self.number == self.dataset.num_iterations:
            return None
        return self.dataset.iteration(self.number + 1)

    def __eq__(self, other):
        return self.number == other.number and self.dataset is other.dataset

    def __repr__(self):
        return f'{self.__class__.__name__}({self.number}, {self.dataset})'


class WESTDataset:

    def __init__(self, h5filename, segment_traj_loader=None):
        self.h5file = WESTPAH5File(h5filename, 'r')
        self.segment_traj_loader = segment_traj_loader

    @property
    def segment_traj_loader(self):
        return self._segment_traj_loader

    @segment_traj_loader.setter
    def segment_traj_loader(self, value):
        if value is not None:
            _ = value(1, 0)
        self._segment_traj_loader = value

    @property
    def num_iterations(self):
        return len(self.h5file['iterations'])

    @property
    def iterations(self):
        return (Iteration(n, self) for n in range(1, self.num_iterations + 1))

    def iteration(self, number):
        valid_range = range(1, self.num_iterations + 1)
        if number not in valid_range:
            raise ValueError(f'iteration number must be in {valid_range}')
        return Iteration(number, self)

    def iteration_group(self, number):
        return self.h5file.get_iter_group(number)

    def __len__(self):
        return self.num_iterations

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.h5file.filename}')"


def trace(segment):
    segments = [segment]
    segment = segment.parent
    while segment:
        segments.append(segment)
        segment = segment.parent
    return segments[::-1]
