import numpy

class Trajectory(object):
    def __init__(self, seg_ids, weight, pcoord,
                 endpoint_type,
                 cputime = None, walltime = None,
                 startdate = None, enddate = None,
                 data = None):
        self.seg_ids = numpy.array(seg_ids, numpy.uint64)
        self.weight = numpy.array(weight, numpy.float64)
        self.pcoord = numpy.array(pcoord)
        self.endpoint_type = endpoint_type
        self.cputime = cputime
        self.walltime = walltime
        self.startdate = startdate
        self.enddate = enddate
        if data is None:
            self.data = numpy.fromiter({} for i in xrange(0, len(self.seg_ids)))
        else:
            self.data = data
        
    def __repr__(self):
        brepr = super(Trajectory,self).__repr__()[:-1]
        
        return '%s: %d segments, leaf_id=%d>' % (brepr, len(self.seg_ids), self.seg_ids[-1])
    
    @classmethod
    def from_segment_list(cls, segments, squeeze_data = True):
        n_segs = len(segments)
        n_pcoord_steps = segments[0].pcoord.shape[0]
        
        if squeeze_data:
            pcoord_shape = ((n_pcoord_steps-1)*(n_segs-1) + n_pcoord_steps,) \
                           + segments[0].pcoord.shape[1:]
        else:
            pcoord_shape = (n_pcoord_steps * n_segs,) \
                           + segments[0].pcoord.shape[1:]
                           
        weight_shape = pcoord_shape[0:1]
        
        seg_ids = numpy.empty((n_segs,), numpy.uint64)
        pcoord = numpy.empty(pcoord_shape, segments[0].pcoord.dtype)
        weight = numpy.empty(weight_shape, numpy.float64)
        data = numpy.empty((n_segs,), numpy.object_)
        
        seg_ids[0] = segments[0].seg_id
        weight[0:n_pcoord_steps] = segments[0].weight
        pcoord[0:n_pcoord_steps] = segments[0].pcoord
        cputime = segments[0].cputime or 0.0
        walltime = segments[0].walltime or 0.0
        data[0] = segments[0].data
        
        for (oiseg, seg) in enumerate(segments[1:]):
            iseg = oiseg+1
            seg_ids[iseg] = seg.seg_id
            data[iseg] = seg.data
            
            if squeeze_data:
                lb = (n_pcoord_steps-1)*(iseg-1) + n_pcoord_steps
                ub = (n_pcoord_steps-1)*iseg     + n_pcoord_steps
                
                # Skip first pcoord row, as we already have it in the last point
                # from the last segment (but double check)
                pco = 1
                assert (seg.pcoord[0] == pcoord[lb-1]).all()
            else:
                lb = (n_pcoord_steps) * (iseg-1)
                ub = (n_pcoord_steps) * iseg
                pco = 0
                
            # Constant weight through the segment
            weight[lb:ub] = seg.weight
            
            pcoord[lb:ub] = seg.pcoord[pco:]
            
            cputime += segments[iseg].cputime
            walltime += segments[iseg].walltime
        
            
        return cls(seg_ids, weight, pcoord, 
                   endpoint_type = segments[-1].endpoint_type, 
                   cputime = cputime, walltime = walltime, 
                   startdate = segments[0].startdate, 
                   enddate = segments[0].enddate,
                   data = data)
