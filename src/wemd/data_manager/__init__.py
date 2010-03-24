import schema, core, mappingtable, versioning

def make_data_manager(runtime_config):
    return core.SQLAlchemyDataManager(runtime_config)        
