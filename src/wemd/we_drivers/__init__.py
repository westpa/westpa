import logging
log = logging.getLogger(__name__)

def get_we_driver(driver_name):
    if driver_name == 'fixed':
        log.info('using fixed-bin WE driver')
        from fixed_bins import FixedBinWEDriver
        driver = FixedBinWEDriver
    elif driver_name == 'fixed_bin_particles':
        log.info('using fixed-bin-particles WE driver')
        log.warn('fixed_bin_particles WE driver is deprecated')
        from fixed_bins_particles import FixedBinParticlesWEDriver
        driver = FixedBinParticlesWEDriver
    else:
        from wemd.core import ConfigError
        raise ConfigError('invalid bin type (%s) specified' % bin_type)
    
    return driver
