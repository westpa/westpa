__metaclass__ = type
import numpy
from math import isnan

class Segment:
    SEG_STATUS_UNSET = None
    SEG_STATUS_PREPARED = 1
    SEG_STATUS_RUNNING  = 2
    SEG_STATUS_COMPLETE = 3
    SEG_STATUS_FAILED   = 4
    
    SEG_ENDPOINT_TYPE_NOTSET = 0
    SEG_ENDPOINT_TYPE_CONTINUES = 1
    SEG_ENDPOINT_TYPE_MERGED = 2
    SEG_ENDPOINT_TYPE_RECYCLED = 3
    
    status_names = {}
    endpoint_type_names = {}
    
    def __init__(self, n_iter = None, seg_id = None, status = None, 
                 n_parents = None, p_parent_id = None, parent_ids = None,
                 endpoint_type = None, weight = None, pcoord = None, walltime = None, cputime = None,
                 data = None):
        # NaNs appear sometimes if a WEMD program is terminated unexpectedly; replace with zero
        walltime = 0.0 if walltime is None or isnan(walltime) else walltime
        cputime  = 0.0 if cputime  is None or isnan(cputime)  else cputime
        
        # the int() and float() calls are required so that new-style string formatting doesn't barf
        # assuming that the respective fields are actually strings, probably after implicitly 
        # calling __str__() on them.  Not sure if this is a numpy, h5py, or python problem
        self.n_iter = int(n_iter)  if n_iter is not None else None
        self.seg_id = long(seg_id) if seg_id is not None else None
        self.status = int(status)  if status is not None else None
        self.p_parent_id = long(p_parent_id) if p_parent_id is not None else None
        self.parent_ids = set(map(long,parent_ids)) if parent_ids else set()
        self.n_parents = int(n_parents) if n_parents else len(self.parent_ids)
        self.endpoint_type = int(endpoint_type) if endpoint_type is not None else None
        self.weight = float(weight) if weight is not None else None
        self.pcoord = numpy.asarray(pcoord) if pcoord is not None else None
        self.walltime = walltime
        self.cputime = cputime
        self.data = data if data else {}

    def __repr__(self):
        return '<%s(%s) n_iter=%r seg_id=%r weight=%r p_parent_id=%r parent_ids=%r>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.n_iter, self.seg_id, self.weight, self.p_parent_id, tuple(sorted(self.parent_ids)))
            
    status_text = property((lambda s: s.status_names[s.status]))
    endpoint_type_text = property((lambda s: s.endpoint_type_names[s.endpoint_type]))
        
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_STATUS_')):
    _val = getattr(Segment, _attr)
    Segment.status_names[_val] = _attr[11:].lower()
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_ENDPOINT_TYPE_')):
    _val = getattr(Segment, _attr)
    Segment.endpoint_type_names[_val] = _attr[18:].lower()    

del _attr, _val
