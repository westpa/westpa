import h5py


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
        return Segment(self.info['parent_id'],
                       self.iteration.prev)

    @property
    def children(self):
        if not self.iteration.next:
            return []
        return [segment for segment in self.iteration.next.segments
                if segment.parent == self]

    def trace(self):
        return trace(self)

    def __eq__(self, other):
        return self.index == other.index and self.iteration == other.iteration

    def __repr__(self):
        return f'{self.__class__.__name__}({self.index}, {self.iteration})'


class Iteration:

    def __init__(self, number,  dataset):
        self.number = number
        self.dataset = dataset

    @property
    def segment_info(self):
        return self.dataset.segment_info(self.number)

    @property
    def pcoord(self):
        return self.dataset.pcoord(self.number)

    @property
    def num_segments(self):
        return self.segment_info.shape[0]

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

    @property
    def segments(self):
        return (Segment(index, self) for index in range(self.num_segments))

    def segment(self, index):
        valid_range = range(self.num_segments)
        if index not in valid_range:
            raise ValueError(f'segment index must be in {valid_range}')
        return Segment(index, self)

    def __eq__(self, other):
        return self.number == other.number and self.dataset is other.dataset

    def __repr__(self):
        return f'{self.__class__.__name__}({self.number}, {self.dataset})'


class WEDataset:

    def __init__(self, h5filename):
        self.h5file = h5py.File(h5filename)

    def _iter_group(self, iteration_number):
        return self.h5file['iterations'][f'iter_{iteration_number:08d}']

    def segment_info(self, iteration_number):
        return self._iter_group(iteration_number)['seg_index']

    def pcoord(self, iteration_number):
        return self._iter_group(iteration_number)['pcoord']

    @property
    def num_iterations(self):
        return len(self.h5file['iterations'])

    @property
    def iterations(self):
        return (Iteration(iteration_number, self)
                for iteration_number in range(1, self.num_iterations + 1))

    def iteration(self, number):
        valid_range = range(1, self.num_iterations + 1)
        if number not in valid_range:
            raise ValueError(f'iteration number must be in {valid_range}')
        return Iteration(number, self)

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
