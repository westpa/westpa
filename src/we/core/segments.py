__metaclass__ = type

class Segment:
    SEG_STATUS_UNSET = None
    SEG_STATUS_PREPARED = 1
    SEG_STATUS_RUNNING = 2
    SEG_STATUS_COMPLETE = 3
    
    status_names = {}

    def __init__(self, 
                 seg_id = None, we_iter = None, status = None,
                 p_parent_id = None, p_parent = None, parents = None,
                 weight = None, final_pcoord = None,
                 data_ref = None,
                 walltime = None, cputime = None,
                 supplementary_data = None):
        
        self.seg_id = seg_id            
        self.we_iter = we_iter
        self.status = status
        self.p_parent_id = p_parent_id
        self.p_parent = p_parent
        self.parents = parents
        self.weight = weight
        self.final_pcoord = final_pcoord
        self.data_ref = data_ref
        self.walltime = walltime
        self.cputime = cputime
        self.supplementary_data = supplementary_data
        
    def __repr__(self):
        return '<%s(%s) seg_id=%s p_parent_id=%s weight=%s>' \
               % (self.__class__.__name__, hex(id(self)),
                  self.seg_id, self.p_parent_id, self.weight)
            
    status_text = property((lambda s: s.status_names[s.status]))
    
for _attr in (attr for attr in dir(Segment) if attr.startswith('SEG_STATUS_')):
    _val = getattr(Segment, _attr)
    Segment.status_names[_val] = _attr[11:].lower()
del _attr, _val
