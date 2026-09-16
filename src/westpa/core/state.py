import numpy as np


class State:
    """Represents a specific configuration of the model being simulated.

    Either `coord` or `file` must be specified. The required parameter and the
    significance of its value are determined by the `propagator <propagators.html>`_. If using a
    built-in propagator type, see its documentation for details on creating
    compatible initial states.

    Parameters
    ----------
    coord : 1-D array-like, optional
        Coordinate vector of the state. The array length (dimension of state
        space) and data type are determined by the propagator.
    file : str, optional
        Address of a file containing coordinate data for the state. The format
        of the address (e.g., absolute path or URL) and the file format are
        determined by the propagator.

    Attributes
    ----------
    coord : numpy.ndarray or None
    file : str or None

    Examples
    --------

    >>> import westpa

    Coordinates stored in memory:

    >>> westpa.State(coord=[0., 0.])
    State(coord=array([0., 0.]))

    Coordinates stored on disk:

    >>> westpa.State(file='/path/to/file.xyz'))
    State(file='/path/to/file.xyz')

    """

    def __init__(self, coord=None, file=None):
        if (coord is None) == (file is None):  # xnor
            raise ValueError("either 'coord' or 'file' must be specified")

        if coord is not None:
            coord = np.asarray(coord)
            if not np.issubdtype(coord.dtype, np.number):
                raise TypeError("scalar type of 'coord' must be numeric")
            if coord.ndim != 1:
                raise ValueError("'coord' must be a 1-D array")
            coord.flags.writeable = False

        self._coord = coord
        self._file = str(file) if file is not None else None

    @property
    def coord(self):
        """Coordinate vector (as a read-only array)."""
        return self._coord

    @property
    def file(self):
        """File containing coordinate data."""
        return self._file

    def __repr__(self):
        if self._coord is not None:
            return type(self).__name__ + '(coord=' + repr(self._coord) + ')'
        else:
            return type(self).__name__ + '(file=' + repr(self._file) + ')'
