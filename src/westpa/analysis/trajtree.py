from collections.abc import Sequence
from itertools import chain
from operator import attrgetter
from typing import NamedTuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pygraphviz as pgv
from matplotlib.ticker import MaxNLocator

from ..core._data_manager import DataManager  # noqa


class SegmentPointer(NamedTuple):
    n_iter: int
    seg_id: int


class TrajectoryTree:
    """Network representation of WE trajectory data.

    Parameters
    ----------
    datafile : str or io.BytesIO
        HDF5 file containing simulation output. Either a pathname (e.g.,
        ``'west.h5'``) or an in-memory stream may be provided.
    load_pcoords : bool, default True
        Whether to load progress coordinates when retrieving trajectory segments.
    load_auxdata : bool, default False
        Whether to load auxiliary data when retrieving trajectory segments.

    Attributes
    ----------
    datafile : str or io.BufferedIOBase
    load_pcoords : bool
    load_auxdata : bool
    n_iters : int

    Methods
    -------
    segment_count
    get_segment
    get_segments
    get_parent_ids
    parent
    children
    trace
    close

    Examples
    --------

    >>> import westpa
    >>> trajtree = westpa.TrajectoryTree('west.h5')

    Retrieve a single segment:

    >>> trajtree.get_segment(1, 0)
    <Segment n_iter=1, seg_id=0, weight=0.2, parent_id=-1 at 0x1694df380>

    Retrieve multiple segments from a given iteration:

    >>> trajtree.get_segments(5, [13, 9, 17])
    [<Segment n_iter=5, seg_id=13, weight=0.02, parent_id=4 at 0x1632a36e0>,
     <Segment n_iter=5, seg_id=9, weight=0.01, parent_id=7 at 0x1632a3680>,
     <Segment n_iter=5, seg_id=17, weight=0.02, parent_id=5 at 0x1632a35f0>]

    Get the parent of a segment:

    >>> segment = trajtree.get_segment(5, 9)
    >>> trajtree.parent(segment)
    <Segment n_iter=4, seg_id=7, weight=0.02, parent_id=14 at 0x179e445c0>

    Trace the lineage of a segment:

    >>> traj = trajtree.trace(segment)
    >>> traj
    <Trajectory with 5 segments at 0x15c8b0500>
    >>> list(traj)
    [<Segment n_iter=1, seg_id=1, weight=0.2, parent_id=-2 at 0x16bcb6b70>,
     <Segment n_iter=2, seg_id=13, weight=0.04, parent_id=1 at 0x16bcb68a0>,
     <Segment n_iter=3, seg_id=14, weight=0.08, parent_id=13 at 0x16bcb6930>,
     <Segment n_iter=4, seg_id=7, weight=0.02, parent_id=14 at 0x16bcb6270>,
     <Segment n_iter=5, seg_id=9, weight=0.01, parent_id=7 at 0x16bb8be30>]

    Close the HDF5 file:

    >>> trajtree.close()

    Use as a context manager to close the HDF5 file automatically:

    >>> with westpa.TrajectoryTree('west.h5') as trajtree:
    ...     segment = trajtree.get_segment(5, 9)
    ...
    >>> segment
    <Segment n_iter=5, seg_id=9, weight=0.01, parent_id=7 at 0x161da3b60>

    """

    def __init__(
        self,
        datafile,
        load_pcoords=True,
        load_auxdata=False,
    ):
        self.data_manager = DataManager(datafile)
        self.data_manager.open_backing(mode='r')

        self._load_pcoords = None
        self._load_auxdata = None

        self.load_pcoords = load_pcoords
        self.load_auxdata = load_auxdata

        graph = nx.DiGraph()

        for n_iter in range(1, self.n_iters + 1):
            seg_ids = range(self.segment_count(n_iter))

            nodes = [SegmentPointer(n_iter, seg_id) for seg_id in seg_ids]
            weights = self.data_manager.get_weights(n_iter, seg_ids)
            for node, weight in zip(nodes, weights):
                graph.add_node(node, subset=n_iter, weight=weight)

            if n_iter == 1:
                continue

            parent_ids = self.data_manager.get_all_parent_ids(n_iter).tolist()
            for node, parent_id in zip(nodes, parent_ids):
                if parent_id < 0:
                    continue
                graph.add_edge(SegmentPointer(n_iter - 1, parent_id), node)

        self._graph = graph

    def close(self):
        """Close the HDF5 file."""
        self.data_manager.close_backing()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    @property
    def datafile(self):
        """HDF5 data file."""
        return self.data_manager.we_h5filename

    @property
    def load_pcoords(self):
        """Whether to load progress coordinates."""
        return self._load_pcoords

    @load_pcoords.setter
    def load_pcoords(self, value):
        if not isinstance(value, bool):
            raise TypeError("'load_pcoords' must be True or False")
        self._load_pcoords = value

    @property
    def load_auxdata(self):
        """Whether to load auxiliary data."""
        return self._load_auxdata

    @load_auxdata.setter
    def load_auxdata(self, value):
        if not isinstance(value, bool):
            raise TypeError("'load_auxdata' must be True or False")
        self._load_auxdata = value

    @property
    def n_iters(self):
        """Number of iterations (layers) in the trajectory tree."""
        return self.data_manager.current_iteration - 1

    def _get_iter_group(self, n_iter):
        if not isinstance(n_iter, int):
            raise TypeError("'n_iter' must be an integer")
        if n_iter not in range(1, iter_stop := self.data_manager.current_iteration):
            raise ValueError(f"'n_iter' must be in range(1, {iter_stop})")
        return self.data_manager.get_iter_group(n_iter)

    def segment_count(self, n_iter=None):
        """Return the total number of segments in a given iteration or in the
        entire trajectory tree.

        Parameters
        ----------
        n_iter : int, optional
            Iteration number. If provided, the total number of segments in the
            given iteration will be returned. If None (the default), the total
            number of segments in the trajectory tree will be returned.

        Returns
        -------
        n_segments : int
            Total number of segments in the given iteration if `n_iter` was
            provided; else the total number of segments in the trajectory tree.

        """
        if n_iter is None:
            return self._graph.number_of_nodes()
        else:
            return self._get_iter_group(n_iter)['seg_index'].shape[0]

    def get_segments(self, n_iter, seg_ids=None):
        """Retrieve multiple segments from a given iteration.

        Parameters
        ----------
        n_iter : int
            Iteration number.
        seg_ids : list of int, optional
            Segment indices. If not provided, all the segments from iteration
            `n_iter` will be returned.

        Returns
        -------
        segments : list of Segment
            Selected segments.

        """
        if seg_ids is not None:
            if not all(isinstance(seg_id, int) for seg_id in seg_ids):
                raise TypeError("'seg_ids' must be a list of integers")
            n_segments = self.segment_count(n_iter)
            for seg_id in seg_ids:
                if seg_id not in range(n_segments):
                    raise ValueError(f'segment index {seg_id} out of range for iteration {n_iter}')

        segments = self.data_manager.get_segments(n_iter, seg_ids, self.load_pcoords, self.load_auxdata)

        if seg_ids is not None:
            # reorder segments to match order of indices in `seg_ids`
            segments_by_id = {s.seg_id: s for s in segments}
            segments = [segments_by_id[seg_id] for seg_id in seg_ids]

        return segments

    def get_segment(self, n_iter, seg_id):
        """Retrieve a single segment.

        Parameters
        ----------
        n_iter : int
            Iteration number.
        seg_id : int
            Segment index.

        Returns
        -------
        segment : Segment
            Selected segment.

        """
        return self.get_segments(n_iter, seg_ids=[seg_id])[0]

    def get_parent_ids(self, n_iter, seg_ids=None):
        """Return the parent indices of all or selected segments from a given iteration.

        Parameters
        ----------
        n_iter : int
            Iteration number.
        seg_ids : list of int, optional
            Indices of the segments for which to retrieve parent indices.
            If not provided, the parent indices of all the segments (in order)
            will be returned.

        Returns
        -------
        parent_ids : list of int
            Parent indices of the selected segments.

        """
        return self.data_manager.get_parent_ids(n_iter, seg_ids)

    def parent(self, segment, n=1):
        """Return the `n`-th level parent of a segment.

        Parameters
        ----------
        segment : Segment
            Segment to find the parent of.
        n : int, default 1
            Number of iterations by which the parent is removed from `segment`
            (1 = direct parent, 2 = grandparent, etc.).

        Returns
        -------
        parent : Segment or None
            `n`-th level parent of the given segment, or None if the segment
            has no `n`-th level parent.

        """
        if not isinstance(n, int):
            raise TypeError("'n' must be an integer")
        if n < 1:
            raise ValueError("'n' must be greater than or equal to 1")

        node = SegmentPointer(segment.n_iter, segment.seg_id)

        for _ in range(n):
            predecessor = next(self._graph.predecessors(node), None)
            if predecessor is None:
                return None
            node = predecessor

        return self.get_segment(node.n_iter, node.seg_id)

    def children(self, segment, n=1):
        """Return the `n`-th level children of a segment.

        Parameters
        ----------
        segment : Segment
            Segment to find the children of.
        n : int, default 1
            Number of iterations by which the children are removed from
            `segment` (1 = direct children, 2 = grandchildren, etc.).

        Returns
        -------
        children : list of Segment
            `n`-th level children of the given segment (empty if the
            segment has no `n`-th level children).

        """
        if not isinstance(n, int):
            raise TypeError("'n' must be an integer")
        if n < 1:
            raise ValueError("'n' must be greater than or equal to 1")

        nodes = [SegmentPointer(segment.n_iter, segment.seg_id)]

        for _ in range(n):
            successors = list(chain.from_iterable(map(self._graph.successors, nodes)))
            if not successors:
                return []
            nodes = successors

        return self.get_segments(nodes[0].n_iter, [node.seg_id for node in nodes])

    def trace(self, segment, maxlen=None):
        """Trace the lineage of a segment.

        Parameters
        ----------
        segment : Segment
            Segment to trace.
        maxlen : int, optional
            Maximum number of segments in the returned trajectory trace.

        Returns
        -------
        traj : Trajectory
            Trajectory leading up to and including `segment`.

        """
        if maxlen is not None:
            if not isinstance(maxlen, int):
                raise TypeError("'maxlen' must be an integer")
            if maxlen < 1:
                raise ValueError("'maxlen' must be positive")

        segments = [segment]
        while parent := self.parent(segment):
            if maxlen is not None and len(segments) == maxlen:
                break
            segments.append(parent)
            segment = parent

        return Trajectory(reversed(segments))

    def to_networkx(self, first_iter=1, last_iter=None, stride=1, roots=None, copy=True):
        """Return a NetworkX representation of the trajectory tree.

        Parameters
        ----------
        first_iter : int, optional
        last_iter : int, optional
        stride : int, optional
        roots : tuple or iterable or tuple, optional
        copy : bool, optional

        Returns
        -------
        graph : networkx.DiGraph
            Directed graph in which nodes represent trajectory segments.

        """
        last_iter = last_iter or self.n_iters
        iter_range = range(first_iter, last_iter + 1, stride)

        nodes = (node for node in self._graph if node.n_iter in iter_range)
        graph = self._graph.subgraph(nodes)

        return graph.copy() if copy else graph

    def view(self, **kwargs):
        """Display an interactive view of the trajectory tree.

        Parameters
        ----------
        **kwargs
            View options. See the :class:`TrajectoryTreeView` class documentation
            for details.

        Returns
        -------
        view : TrajectoryTreeView
            Interface for querying and modifying the interactive view.

        """
        return TrajectoryTreeView(self, **kwargs)

    def __repr__(self):
        return f'<{type(self).__name__} with {self.n_iters} iterations at {hex(id(self))}>'


def _to_pygraphviz(trajtree, first_iter=1, last_iter=None):
    graph = trajtree.to_networkx(first_iter=first_iter, last_iter=last_iter, copy=False)

    agraph = pgv.AGraph(directed=True, ordering='out')
    subgraphs = {}
    for u, d in graph.nodes(data=True):
        n_iter = d['subset']
        if n_iter not in subgraphs:
            subgraphs[n_iter] = agraph.add_subgraph(rank='same')
        subgraphs[n_iter].add_node(u, weight=d['weight'])
    for u, v in graph.edges():
        agraph.add_edge(u, v)

    return agraph


class TrajectoryTreeView:
    """Interactive view of a trajectory tree.

    Parameters
    ----------
    trajtree : TrajectoryTree
    ax : matplotlib.axes.Axes, optional
    first_iter : int, optional
    last_iter : int, optional
    stride : int, optional
    roots : tuple or iterable of tuple, optional
    x_func : callable, optional
    x_label : str, optional
    cmap : str, optional
    default_edge_color : str, optional
    highlight_color : str, optional

    """

    def __init__(
        self,
        trajtree,
        ax=None,
        first_iter=1,
        last_iter=None,
        stride=1,
        roots=None,
        x_func=None,
        x_label=None,
        cmap='viridis_r',
        default_edge_color='gray',
        highlight_color='magenta',
    ):
        self._trajtree = trajtree
        self._ax = ax
        self._first_iter = first_iter
        self._last_iter = last_iter or self._trajtree.n_iters
        self._stride = stride
        self._roots = roots
        self._x_func = x_func
        self._x_label = x_label
        self._cmap = cmap
        self._default_edge_color = default_edge_color
        self._highlight_color = highlight_color

        self._fig = None

        self._graph = None
        self._nodes = None
        self._x_pos = None
        self._highlighted_subgraph = None
        self._edge_artists = None
        self._node_artist = None
        self._highlighted_node_artist = None
        self._text = None

        self._draw()

    def _clear(self):
        for artist in self._edge_artists.values():
            artist.remove()
        self._node_artist.remove()
        self._highlighted_node_artist.remove()

        self._graph = None
        self._nodes = None
        self._x_pos = None
        self._highlighted_subgraph = None
        self._edge_artists = None
        self._node_artist = None
        self._highlighted_node_artist = None
        self._text = None

    def _draw(self):
        if (ax := self._ax) is None:
            fig, ax = plt.subplots()
            ax.set_ylabel('WE Iteration')
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            fig = ax.get_figure()
            if self._graph is not None:
                self._clear()

        graph = self.trajtree.to_networkx(self.first_iter, self.last_iter, copy=False)

        x_pos = {}
        if self.x_func is not None:
            for n_iter in range(self.first_iter, self.last_iter + 1):
                for segment in self.trajtree.get_segments(n_iter):
                    u = SegmentPointer(segment.n_iter, segment.seg_id)
                    x_pos[u] = self.x_func(segment)
            ax.set_xlabel(self.x_label or 'x')
        else:
            agraph = _to_pygraphviz(self.trajtree, self.first_iter, self.last_iter)
            agraph.graph_attr['rankdir'] = 'BT'
            agraph.layout(prog='dot')  # see https://graphviz.org/docs/layouts/dot/
            for u, au in zip(graph, agraph):
                x_pos[u] = float(au.attr['pos'].split(',')[0])
            ax.set_xticks([])

        edge_artists = {}
        for u, v in graph.edges():
            artist, *_ = ax.plot(
                [x_pos[u], x_pos[v]],
                [u.n_iter, v.n_iter],
                c=self.default_edge_color,
                linewidth=0.75,
                zorder=1,
            )
            edge_artists[u, v] = artist

        log_weights = np.log([float(d['weight']) for u, d in graph.nodes(data=True)])

        norm = mpl.colors.Normalize(vmin=log_weights.min(), vmax=log_weights.max())
        cmap = mpl.colormaps[self.cmap]

        y = list(map(attrgetter('n_iter'), graph.nodes()))
        c = list(map(cmap, map(norm, log_weights)))
        node_artist = ax.scatter([x_pos[u] for u in graph], y, c=c, s=3, picker=True)
        highlighted_node_artist = ax.scatter(x=[], y=[], c=self.highlight_color, s=3, alpha=0.8)

        text_ax = fig.add_axes((0.1, 0.88, 0.8, 0.1))
        text_ax.axis('off')
        text = text_ax.text(0.5, 0.25, '', ha='center')
        text.set_visible(False)

        fig.canvas.mpl_connect('pick_event', self._node_pick)
        fig.canvas.mpl_connect('motion_notify_event', self._node_hover)

        self._ax = ax
        self._graph = graph
        self._nodes = list(graph)
        self._x_pos = x_pos
        self._highlighted_subgraph = None
        self._edge_artists = edge_artists
        self._node_artist = node_artist
        self._highlighted_node_artist = highlighted_node_artist
        self._text = text

    @property
    def trajtree(self):
        return self._trajtree

    @property
    def first_iter(self):
        return self._first_iter

    @property
    def last_iter(self):
        return self._last_iter

    @property
    def stride(self):
        return self._stride

    @property
    def roots(self):
        return self._roots

    @property
    def x_func(self):
        return self._x_func

    @property
    def x_label(self):
        return self._x_label

    @property
    def cmap(self):
        return self._cmap

    @property
    def default_edge_color(self):
        return self._default_edge_color

    @property
    def highlight_color(self):
        return self._highlight_color

    def clear_highlights(self):
        """Clear highlighted nodes and edges."""
        if self._highlighted_subgraph is not None:
            for edge in self._highlighted_subgraph.edges():
                self._edge_artists[edge].set_color(self._default_edge_color)
            self._highlighted_node_artist.set_offsets([float('inf'), float('inf')])
            self._highlighted_subgraph = None

    def _highlight(self, subgraph):
        if self._highlighted_subgraph is not None:
            self.clear_highlights()
        for edge in subgraph.edges():
            self._edge_artists[edge].set_color(self.highlight_color)
        self._highlighted_node_artist.set_offsets([[self._x_pos[u], u.n_iter] for u in subgraph.nodes()])
        self._highlighted_subgraph = subgraph

    def highlight_trace(self, n_iter, seg_id):
        """Highlight the trace of a given segment.

        Parameters
        ----------
        n_iter : int
        seg_id : int

        """
        u = SegmentPointer(n_iter, seg_id)
        subgraph = self._graph.subgraph(nx.ancestors(self._graph, u) | {u})
        self._highlight(subgraph)

    def highlight_subtree(self, n_iter, seg_id):
        """Highlight the subtree rooted at a given segment.

        Parameters
        ----------
        n_iter : int
        seg_id : int

        """
        u = SegmentPointer(n_iter, seg_id)
        subgraph = self._graph.subgraph(nx.descendants(self._graph, u) | {u})
        self._highlight(subgraph)

    def _node_pick(self, event):
        self.clear_highlights()
        n_iter, seg_id = self._nodes[event.ind[0]]
        if event.mouseevent.button == plt.MouseButton.LEFT:
            if event.mouseevent.dblclick:
                self.highlight_subtree(n_iter, seg_id)
            else:
                self.highlight_trace(n_iter, seg_id)

    def _node_hover(self, event):
        contains, details = self._node_artist.contains(event)
        if contains:
            u = self._nodes[details['ind'][0]]
            self._text.set_text(f'(n_iter, seg_id) = ({u.n_iter}, {u.seg_id})')
            self._text.set_visible(True)
        else:
            self._text.set_visible(False)

    def highlighted_segments(self):
        if self._highlighted_subgraph is not None:
            iter_range = range(self.first_iter, self.last_iter + 1)
            seg_ids_by_iter = {n_iter: [] for n_iter in iter_range}
            for u in self._highlighted_subgraph:
                seg_ids_by_iter[u.n_iter].append(u.seg_id)
            for n_iter in iter_range:
                seg_ids = sorted(seg_ids_by_iter[n_iter])
                yield from self._trajtree.get_segments(n_iter, seg_ids)
        else:
            yield from ()


class Trajectory(Sequence):
    """A contiguous sequence of trajectory segments.

    Parameters
    ----------
    segments : iterable of Segment
        Segments comprising the trajectory.

    Attributes
    ----------
    states : iterator of State
    initial_state : State
    final_state : State
    iter_range : range

    """

    def __init__(self, segments):
        self.segments = tuple(segments)

    def __len__(self):
        return len(self.segments)

    def __getitem__(self, index):
        return self.segments[index]

    @property
    def states(self):
        """Sequence of states visited by the trajectory."""
        for segment in self.segments:
            yield segment.initial_state
        yield self.segments[-1].final_state

    @property
    def initial_state(self):
        """Initial state of the trajectory."""
        return self.segments[0].initial_state

    @property
    def final_state(self):
        """Final state of the trajectory."""
        return self.segments[-1].final_state

    @property
    def iter_range(self):
        """Range of iterations spanned by the trajectory."""
        return range(self.segments[0].n_iter, self.segments[-1].n_iter + 1)

    def __repr__(self):
        return f'<{type(self).__name__} with {len(self)} segments at {hex(id(self))}>'
