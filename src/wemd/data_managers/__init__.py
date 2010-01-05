
class DataManagerBase(object):
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        
    def prepare_backing(self, sim_config):
        raise NotImplementedError

    def create_we_sim_iter(self, we_sim_iter):
        raise NotImplementedError
    
    def update_we_sim_iter(self, we_sim_iter):
        raise NotImplementedError
    
    def get_we_sim_iter(self, we_sim_iter_id):
        raise NotImplementedError
    
    def create_segment(self, segment):
        raise NotImplementedError
    
    def update_segment(self, segment):
        raise NotImplementedError

def make_data_manager(runtime_config):
    storage_engine = runtime_config.get('data.storage_engine', 'hdf5')
    
    if storage_engine == 'hdf5':
        import hdf5_data_manager
        return hdf5_data_manager.HDF5DataManager(runtime_config)
    else:
        from wemd.core.errors import ConfigError
        raise ConfigError('invalid data storage engine %r' % storage_engine)
