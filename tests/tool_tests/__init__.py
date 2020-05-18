from westpa import rc, h5io

data_manager = rc.get_data_manager()

##Store west.h5 file in RAM for testing
west_file_name = 'west.h5'
west_file = h5io.WESTPAH5File(west_file_name, driver='core', backing_store=False)

data_manager.we_h5file = west_file
data_manager.we_h5file_version = int(west_file['/'].attrs['west_file_format_version'])
