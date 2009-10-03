from wemd.util.lazy_loader import lazy_load
import sa

__metaclass__ = type
    
class DataManagerBase:
    we_sim = lazy_load('_we_sim')
    
    def __init__(self, source):
        self.source = source
        self.config = {}
                
    def configure(self, config):
        self.config.update(config)

    def initialize(self):
        raise NotImplementedError
        
    def store_state_object(self, name, obj):
        raise NotImplementedError
        
    def load_state_object(self, name):
        raise NotImplementedError
    
    def __lazy_load__(self, name):
        return self.load_state_object(name[1:])
    
    def save_state(self):
        assert(self._we_sim.data_manager is self)
        self._we_sim.data_manager = None
        try:
            self.store_state_object('we_sim', self._we_sim)
        finally:
            self._we_sim.data_manager = self
        self.save_config()
    
    def restore_state(self):
        self.restore_config()
        self._we_sim = self.load_state_object('we_sim')
        self._we_sim.data_manager = self
    
    def save_config(self):
        raise NotImplementedError

    def restore_config(self):
        raise NotImplementedError
                                                      
    def is_propagation_complete(self):
        """Returns True if propagation of all segments for this iteration is
        complete"""
        raise NotImplementedError

    def create_segments(self, segments):
        raise NotImplementedError
    
    def get_segment(self, we_iter, seg_id, load_parent = False):
        raise NotImplementedError
    
    def get_segments(self, we_iter):
        raise NotImplementedError    
    
    def get_next_segment(self, pretend = False):
        raise NotImplementedError
    
    def record_lineage(self, segment, parent_particles):
        raise NotImplementedError
    
    def record_data_item(self, we_iter, name, obj):
        raise NotImplementedError
    
    def load_data_item(self, we_iter, name):
        raise NotImplementedError
        
    def load_data_sequence(self, name):
        raise NotImplementedError
    
    def get_trajectory(self, we_iter, seg_id):
        node = self.get_segment(we_iter, seg_id, load_parent=False)
        traj = [node]
        while node.p_parent_id is not None:
            node = self.get_segment(node.we_iter-1, node.p_parent_id, False)
            traj.append(node)
        return list(reversed(traj))

    def get_trajectories(self):
        we_iter = self.we_sim.current_iteration-1
        segs = self.get_segments(we_iter)
        trajs = [self.get_trajectory(we_iter, seg.seg_id) for seg in segs]
        return trajs
        
    def scrub_crashed_segments(self, we_iter):
        raise NotImplementedError

from sqlite import SQLiteDataManager

def make_data_manager(source):
    return SQLiteDataManager(source)
