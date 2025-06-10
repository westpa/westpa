import westpa
from westpa.core.yamlcfg import check_bool

class AssignPorts:
    def __init__(self, sim_manager, plugin_config):
        if not sim_manager.work_manager.is_master:
            return
        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system

        # Switch for the plugin    
        self.assign_ports = check_bool(plugin_config.get('assign_ports', False))

        if self.assign_ports:
            sim_manager.register_callback(sim_manager.pre_propagation, self.pre_propagation, 0)

    def pre_propagation(self):
        n_segs = len(self.sim_manager.segments)
        westpa.rc.pstatus('-westext.assign_port -----------------\n')
        westpa.rc.pstatus(f'Assigning ports to {n_segs} segments.\n')
        westpa.rc.pflush()

        # for seg in self.sim_manager.segments:

def assign_port(seg_id):
    """
    Assign a port to a segment based on its ID.
    """
    
    return seg_id