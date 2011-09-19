from __future__ import division, print_function; __metaclass__ = type

import logging
log = logging.getLogger(__name__)

import numpy
import wemd
from wemdtools.aframe import AnalysisMixin

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
        
    cmap_data = numpy.array( [ (124,   0,  24),
                               (211,   0,  32),
                               (244, 184,   0),
                               (245, 235,   0),
                               (129, 183,   2),
                               ( 32, 128,  38),
                               ( 21,  27,  87),
                               ( 36,  62, 137),
                               (178, 220, 245),
                               (255, 255, 255) ], numpy.float32)
    cmap_data /= 255.0
    cm_hovmol   = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol', cmap_data)
    cm_hovmol_r = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol_r', numpy.flipud(cmap_data))
        

class PlottingMixin(AnalysisMixin):
    def __init__(self):
        global matplotlib, pyplot
        
        super(PlottingMixin,self).__init__()
        
        self.matplotlib_avail = (matplotlib is not None and pyplot is not None)
        
    def require_matplotlib(self):
        global matplotlib
        if not self.matplotlib_avail:
            raise RuntimeError('matplotlib is not available')
        else:
            return matplotlib

