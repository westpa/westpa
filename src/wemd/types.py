__metaclass__ = type
import numpy

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
                 endpoint_type = None, weight = None, pcoord = None, walltime = None,
                 cputime = None):
        self.n_iter = n_iter
        self.seg_id = seg_id
        self.status = status
        self.p_parent_id = p_parent_id
        self.parent_ids = set(parent_ids) if parent_ids else set()
        self.n_parents = n_parents or len(self.parent_ids)
        self.endpoint_type = endpoint_type
        self.weight = weight
        self.pcoord = numpy.asarray(pcoord) if pcoord is not None else None
        self.walltime = walltime
        self.cputime = cputime

    def __repr__(self):
        return '<%s(%s) seg_id=%r weight=%r p_parent_id=%r parent_ids=%r>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.seg_id, self.weight, self.p_parent_id, tuple(sorted(self.parent_ids)))
            
    status_text = property((lambda s: s.status_names[s.status]))
    endpoint_type_text = property((lambda s: s.endpoint_type_names[s.endpoint_type]))
    
    def update_propagation_data(self, segment):
        '''Update this segment with propagation data from the given segment object. 
        Non-propagation data (seg_id, parentage, n_iter, etc.) is IGNORED.'''            
        self.status = segment.status
        self.walltime = segment.walltime
        self.cputime = segment.cputime
        self.pcoord[...] = segment.pcoord[...]
    
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_STATUS_')):
    _val = getattr(Segment, _attr)
    Segment.status_names[_val] = _attr[11:].lower()
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_ENDPOINT_TYPE_')):
    _val = getattr(Segment, _attr)
    Segment.endpoint_type_names[_val] = _attr[18:].lower()    

del _attr, _val


class Particle:
    GEN_CONTINUATION = 0
    GEN_SPLIT = 1
    GEN_MERGE = 2
    GEN_RECYCLE = 3
    
    def __init__(self, seg_id = None, weight = None, pcoord = None, 
                 p_parent_id = None, parent_ids = None, source_id = None, genesis = None):
        self.seg_id = seg_id
        self.p_parent_id = p_parent_id
        self.parent_ids = set(parent_ids) if parent_ids else set()
        self.weight = weight
        self.pcoord = numpy.asarray(pcoord)
        self.source_id = source_id
        self.genesis = genesis
                
    def __repr__(self):
        return '<%s(%s) seg_id=%r weight=%r p_parent_id=%r parent_ids=%r>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.seg_id, self.weight, self.p_parent_id, tuple(sorted(self.parent_ids)))
