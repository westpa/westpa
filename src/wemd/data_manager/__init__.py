import hdf5
def make_data_manager(runtime_config):
    return hdf5.HDF5DataManager(runtime_config)
