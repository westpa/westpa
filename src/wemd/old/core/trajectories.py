import numpy

class Trajectory(object):
    def __init__(self, segments, squeeze_data = True):
        self.segments = list(segments)
        self._squeeze_data = True
        
        self.total_cputime = numpy.sum(numpy.fromiter((segment.cputime for segment in segments),
                                                      numpy.float64))
        self.total_walltime = numpy.sum(numpy.fromiter((segment.walltime for segment in segments),
                                                       numpy.float64))
                        
        self._seg_ids = None
        self._pcoord = None
        self._weight = None
        self._data = None
    
    def __repr__(self):
        return '<Trajectory leaf_id=%d len=%d>' % (self.segments[-1].seg_id,
                                                   len(self.segments))
        
    def __hash__(self):
        return hash(self.segments[-1].seg_id)
    
    def __len__(self):
        return len(self.segments)
                                                   
    
    def _get_squeeze_data(self):
        return self._squeeze_data

    def _set_squeeze_data(self, val):
        if self._squeeze_data != val:
            self._squeeze_data = val
            self.clear_arrays()
        
    def populate_arrays(self):
        n_segs = len(self.segments)
        n_pcoord_steps = self.segments[0].pcoord.shape[0]
        
        if self._squeeze_data:
            pcoord_shape = ((n_pcoord_steps-1)*(n_segs-1)+n_pcoord_steps,)  \
                           + self.segments[0].pcoord.shape[1:]
        else:
            pcoord_shape = (n_pcoord_steps * n_segs,) \
                           + self.segments[0].pcoord.shape[1:]
                           
        weight_shape = pcoord_shape[0:1]
        
        seg_ids = numpy.fromiter((segment.seg_id for segment in self.segments),
                                 numpy.uint64)
        data = numpy.array([segment.data for segment in self.segments],
                              numpy.object_)
        
        pcoord = numpy.empty(pcoord_shape, self.segments[0].pcoord.dtype)
        weight = numpy.empty(weight_shape, numpy.float64)
        
        weight[0:n_pcoord_steps] = self.segments[0].weight
        pcoord[0:n_pcoord_steps] = self.segments[0].pcoord
        
        n_squeeze_data = int(self._squeeze_data)
        for (oiseg, segment) in enumerate(self.segments[1:]):
            iseg = oiseg+1
            
#            This indexing trick reduces to the following:
#            if squeeze_data:
#                lb = (n_pcoord_steps-1)*(iseg-1) + n_pcoord_steps
#                ub = (n_pcoord_steps-1)*iseg     + n_pcoord_steps
#                
#                # Skip first pcoord row, as we already have it in the last point
#                # from the last segment (but double check)
#                pco = 1
#                #assert (seg.pcoord[0] == pcoord[lb-1]).all()
#            else:
#                lb = (n_pcoord_steps) * (iseg-1)
#                ub = (n_pcoord_steps) * iseg
#                pco = 0

            lb = (n_pcoord_steps - n_squeeze_data)*(iseg-1) + n_squeeze_data*n_pcoord_steps
            ub = (n_pcoord_steps - n_squeeze_data)*iseg     + n_squeeze_data*n_pcoord_steps
            pco = n_squeeze_data
            
            weight[lb:ub] = segment.weight
            pcoord[lb:ub] = segment.pcoord[pco:]
        self._seg_ids = seg_ids
        self._data = data
        self._pcoord = pcoord
        self._weight = weight
                           
    def clear_arrays(self):
        self._seg_ids = None
        self._pcoord = None
        self._weight = None
        self._data = None
        
    def get_seg_ids(self):
        if self._seg_ids is None:
            self.populate_arrays()
        return self._seg_ids
    
    def get_pcoord(self):
        if self._pcoord is None:
            self.populate_arrays()
        return self._pcoord
    
    def get_weight(self):
        if self._weight is None:
            self.populate_arrays()
        return self._weight
    
    def get_data(self):
        if self._data is None:
            self.populate_arrays()
        return self._data
            
    seg_ids = property(get_seg_ids)
    pcoord = property(get_pcoord)
    weight = property(get_weight)
    data = property(get_data)
