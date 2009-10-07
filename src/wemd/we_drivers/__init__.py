from wemd.core import ConfigError
from fixed_bins import FixedBinWEDriver

def make_we_driver(sim_config):
    bin_type = sim_config['bins.type']
    if bin_type == 'fixed':
        driver = FixedBinWEDriver()
    else:
        raise ConfigError('invalid bin type (%s) specified' % bin_type)
    
    return driver


        