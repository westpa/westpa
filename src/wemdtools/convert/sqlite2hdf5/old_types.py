__metaclass__ = type

class OldSegment:
    SEG_STATUS_UNSET = None
    SEG_STATUS_PREPARED = 1
    SEG_STATUS_RUNNING  = 2
    SEG_STATUS_COMPLETE = 3
    SEG_STATUS_FAILED   = 4
    SEG_STATUS_DELETED  = 5
    
    SEG_ENDPOINT_TYPE_UNKNOWN = 0
    SEG_ENDPOINT_TYPE_CONTINUATION = 1
    SEG_ENDPOINT_TYPE_MERGED = 2
    SEG_ENDPOINT_TYPE_RECYCLED = 3
    
    status_names = {}
    endpoint_type_names = {}
    
    def __hash__(self):
        return hash(self.seg_id)
    
    def __init__(self, 
                 seg_id = None, n_iter = None, status = None,
                 p_parent_id = None, parent_ids = None,
                 p_parent = None, parents = None,
                 endpoint_type = None,
                 weight = None, pcoord = None,
                 data_ref = None,
                 walltime = None, cputime = None,
                 starttime = None, endtime = None,
                 data = None):
        
        self.seg_id = seg_id
        self.n_iter = n_iter
        self.status = status
        self.p_parent = p_parent
        self.parents = parents or set()
        self.endpoint_type = endpoint_type
        self.weight = weight
        self.pcoord = pcoord
        self.data_ref = data_ref
        self.walltime = walltime
        self.cputime = cputime
        self.starttime = starttime
        self.endtime = endtime
        self.data = data or {}
        
    def __repr__(self):
        return '<%s(%s) seg_id=%s p_parent_id=%s weight=%s>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.seg_id, self.p_parent_id, self.weight)
            
    status_text = property((lambda s: s.status_names[s.status]))
    endpoint_type_text = property((lambda s: s.endpoint_type_names[s.endpoint_type]))
    
for _attr in (attr for attr in dir(OldSegment) if attr.startswith('SEG_STATUS_')):
    _val = getattr(OldSegment, _attr)
    OldSegment.status_names[_val] = _attr[11:].lower()
for _attr in (attr for attr in dir(OldSegment) if attr.startswith('SEG_ENDPOINT_TYPE_')):
    _val = getattr(OldSegment, _attr)
    OldSegment.endpoint_type_names[_val] = _attr[18:].lower()    

del _attr, _val

class WESimIter:
    """
    Describes per-iteration information (summary or otherwise) for
    a WE simulation.
    """
    
    def __init__(self, n_iter = None, n_particles = None, norm = None,
                 cputime = None, walltime = None, 
                 data = None):
        self.n_iter = n_iter
        self.n_particles = n_particles
        self.norm = norm
        self.cputime = cputime
        self.walltime = walltime
        self.data = data or {}


