__metaclass__ = type
import numpy
from numpy import asarray

class Segment:
    SEG_STATUS_UNSET = None
    SEG_STATUS_PREPARED = 1
    SEG_STATUS_RUNNING  = 2
    SEG_STATUS_COMPLETE = 3
    SEG_STATUS_FAILED   = 4
    
    SEG_ENDPOINT_TYPE_UNKNOWN = 0
    SEG_ENDPOINT_TYPE_CONTINUATION = 1
    SEG_ENDPOINT_TYPE_MERGED = 2
    SEG_ENDPOINT_TYPE_RECYCLED = 3
    
    status_names = {}
    endpoint_type_names = {}
    
    def __hash__(self):
        return hash(self.seg_id)

    def __init__(self, seg_id = None, status = None, p_parent = None, parents = None,
                 endpoint_type = None, weight = None, pcoord = None, walltime = None,
                 cputime = None):
        self.seg_id = seg_id
        self.status = status
        self.p_parent = p_parent
        self.parents = set(parents) if parents else set()
        self.endpoint_type = endpoint_type
        self.weight = weight
        self.pcoord = asarray(pcoord) if pcoord is not None else None
        self.walltime = walltime
        self.cputime = cputime

    def __repr__(self):
        return '<%s(%s) seg_id=%s weight=%s>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.seg_id, self.weight)
            
    status_text = property((lambda s: s.status_names[s.status]))
    endpoint_type_text = property((lambda s: s.endpoint_type_names[s.endpoint_type]))
    
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_STATUS_')):
    _val = getattr(Segment, _attr)
    Segment.status_names[_val] = _attr[11:].lower()
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_ENDPOINT_TYPE_')):
    _val = getattr(Segment, _attr)
    Segment.endpoint_type_names[_val] = _attr[18:].lower()    

del _attr, _val


class Particle:
    def __init__(self):
        self.particle_id = None
        self.p_parent = None
        self.parents = None
        self.weight = None
        self.pcoord = None
        
