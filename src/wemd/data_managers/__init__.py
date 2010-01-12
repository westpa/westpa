
class DataManagerBase(object):
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config
        
    def prepare_backing(self, sim_config):
        raise NotImplementedError

    def create_we_sim_iter(self, we_sim_iter):
        raise NotImplementedError
    
    def update_we_sim_iter(self, we_sim_iter):
        raise NotImplementedError
    
    def get_we_sim_iter(self, n_iter):
        raise NotImplementedError
    
    def create_segments(self, we_sim_iter, segments):
        raise NotImplementedError
    
    def update_segments(self, we_sim_iter, segments, **kwargs):
        raise NotImplementedError
    
    def num_incomplete_segments(self, we_iter):
        raise NotImplementedError
    
    def get_segment(self, we_iter, seg_id, **kwargs):
        raise NotImplementedError
    
    def get_segments(self, we_iter, **kwargs):
        raise NotImplementedError
    
    def get_prepared_segments(self, we_iter, **kwargs):
        raise NotImplementedError
        
def make_data_manager(runtime_config):
    storage_engine = runtime_config.get('data.storage_engine', 'hdf5')
    
    if storage_engine == 'sql':
        import sa_data_manager
        return sa_data_manager.SQLAlchemyDataManager(runtime_config)        
    else:
        from wemd.core.errors import ConfigError
        raise ConfigError('invalid data storage engine %r' % storage_engine)
