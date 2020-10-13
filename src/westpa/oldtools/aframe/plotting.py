import logging

import numpy as np

from westpa.oldtools.aframe import AnalysisMixin

log = logging.getLogger(__name__)

try:
    import matplotlib
except ImportError:
    matplotlib = None
    pyplot = None
    cm_hovmol = None
    cm_hovmol_r = None
else:
    try:
        matplotlib.use('PDF')
        from matplotlib import pyplot
    except Exception as e:
        log.info('could not select matplotlib PDF backend: {}'.format(e))
        matplotlib = None
        pyplot = None

    cmap_data = np.array(
        [
            (124, 0, 24),
            (211, 0, 32),
            (244, 184, 0),
            (245, 235, 0),
            (129, 183, 2),
            (32, 128, 38),
            (21, 27, 87),
            (36, 62, 137),
            (178, 220, 245),
            (255, 255, 255),
        ],
        np.float32,
    )
    cmap_data /= 255.0
    cm_hovmol = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol', cmap_data)
    cm_hovmol_r = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol_r', np.flipud(cmap_data))

    _pdr_data = {
        'red': [
            (0.0, 1.0, 1.0),
            (0.25, 1.0, 1.0),
            (0.50, 1.0, 1.0),
            (0.75, 0.0, 0.0),
            (0.90, 0.50, 0.50),
            (0.95, 0.0, 0.0),
            (1.0, 1.0, 1.0),
        ],
        'green': [
            (0.0, 1.0, 1.0),
            (0.25, 0.0, 0.0),
            (0.50, 1.0, 1.0),
            (0.75, 1.0, 1.0),
            (0.90, 0.86, 0.86),
            (0.95, 0.0, 0.0),
            (1.0, 1.0, 1.0),
        ],
        'blue': [(0.0, 1.0, 1.0), (0.25, 0.0, 0.0), (0.75, 0.0, 0.0), (0.90, 1.0, 1.0), (0.95, 1.0, 1.0), (1.0, 1.0, 1.0)],
    }

    cm_pdr = matplotlib.colors.LinearSegmentedColormap('pdr', _pdr_data, 2048)
    cm_pdr_r = matplotlib.colors.LinearSegmentedColormap('pdr_r', matplotlib.cm.revcmap(_pdr_data), 2048)

    matplotlib.cm.register_cmap('pdr', cm_pdr)
    matplotlib.cm.register_cmap('pdr_r', cm_pdr_r)
    matplotlib.cm.register_cmap('hovmol', cm_hovmol)
    matplotlib.cm.register_cmap('hovmol_r', cm_hovmol_r)

    del cmap_data


class PlottingMixin(AnalysisMixin):
    def __init__(self):
        global matplotlib, pyplot

        super().__init__()

        self.matplotlib_avail = matplotlib is not None and pyplot is not None

    def require_matplotlib(self):
        global matplotlib
        if not self.matplotlib_avail:
            raise RuntimeError('matplotlib is not available')
        else:
            return matplotlib
