import wemd
from wt2.tool_classes.selected_segs import AllSegmentSelection
from _trajtree import _trajtree_base

class TrajTreeSet(_trajtree_base):    
    def __init__(self, segsel = None, data_manager = None):
        self.data_manager = data_manager or wemd.rc.get_data_manager()
        self.segsel = segsel or AllSegmentSelection(data_manager = self.data_manager)
         
        #self._build_table(self.segsel, self.data_manager)
        self._build_table()
        
    def _build_table(self):
        _trajtree_base._build_table(self, self.segsel, self.data_manager)    
    
    
#    def _build_table(self):
#        del self.trajtable, self.childtable, self.iter_offsets
#        
#        segsel = self.segsel
#        start_iter = segsel.start_iter
#        
#        seg_count = len(segsel)
#        trajtable  = self.trajtable = numpy.empty((seg_count,), dtype=_tt_dtype)
#        len_childtable = 0
#        childtable = self.childtable = numpy.empty((CT_CHUNKSIZE,), dtype=numpy.uint64)
#        iter_offsets = self.iteroffsets = {}
#        
#        data_manager = self.data_manager
#        tt_offset = 0
#        
#        last_iter_indices = {} # mapping of seg_id -> table row for previous iteration
#        
#        for n_iter in xrange(start_iter, segsel.stop_iter):
#            iter_offsets[n_iter] = tt_offset
#            this_iter_indices = {}
#            
#            seg_ids = segsel.from_iter(n_iter)
#            iter_group = data_manager.get_iter_group(n_iter)
#            seg_index = iter_group['seg_index']
#            weights = seg_index['weight']
#            all_parent_ids = data_manager.get_all_parent_ids(n_iter)
#            
#            for seg_id in seg_ids:
#                this_iter_indices[seg_id] = tt_offset
#                trajtable[tt_offset]['n_iter'] = n_iter
#                trajtable[tt_offset]['seg_id'] = seg_id
#                trajtable[tt_offset]['weight'] = weights[seg_id]
#                trajtable[tt_offset]['n_children'] = 0
#                trajtable[tt_offset]['children_offset'] = 0
#                
#                if n_iter == start_iter:
#                    parent_row = NO_PARENT
#                else:
#                    parent_id = all_parent_ids[seg_id]
#                    
#                    if parent_id < 0:
#                        parent_row = NO_PARENT
#                    else:                        
#                        # this will raise KeyError if segsel is missing segments along a trajectory
#                        parent_row = trajtable[tt_offset]['parent_offset'] = last_iter_indices[parent_id]
#                
#                trajtable[tt_offset]['parent_offset'] = parent_row
#                    
#                if parent_row != NO_PARENT:
#                    if trajtable[parent_row]['n_children'] == 0:
#                        # First child
#                        trajtable[parent_row]['n_children'] = 1
#                        trajtable[parent_row]['children_offset'] = len_childtable
#                    else:
#                        trajtable[parent_row]['n_children'] += 1
#                    
#                    if len_childtable == len(childtable):
#                        childtable.resize((len(childtable)+CT_CHUNKSIZE,))
#                        childtable[len_childtable] = tt_offset
#                        len_childtable += 1
#                
#                tt_offset += 1
#            del all_parent_ids, weights, seg_index, iter_group, seg_ids                     
#            last_iter_indices = this_iter_indices
#            
#        childtable.resize((len(childtable),))
                
    def __len__(self):
        return len(self.trajtable)
