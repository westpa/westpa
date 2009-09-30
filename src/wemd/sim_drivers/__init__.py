from wemd.core import ConfigError
from fixed_bins import FixedBinWESim

def make_sim_driver(data_manager, config):
    bin_type = data_manager.config['bins.type'] = config['bins.type']
    if bin_type == 'fixed':
        return FixedBinWESim(data_manager)
    else:
        raise ConfigError('invalid bin type (%s) specified' % bin_type)


        