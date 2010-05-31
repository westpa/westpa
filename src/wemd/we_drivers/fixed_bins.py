import re

import logging
log = logging.getLogger('wemd.we_drivers.fixed_bins')

import numpy
from we_driver import WEDriver
from wemd.core.binarrays import Bin, BinArray
from wemd.core import ConfigError

class FixedBinWEDriver(WEDriver):
    def __init__(self):
        super(FixedBinWEDriver,self).__init__()
        
    def sim_init(self, sim_config, sim_config_src):
        super(FixedBinWEDriver,self).sim_init(sim_config, sim_config_src)
        bintype = sim_config['bins.type']
        assert bintype in ('fixed', 'fixed_bin_particles')

        ndim = sim_config['bins.ndim'] = sim_config_src.get_int('bins.ndim')
        bin_limits = []
        
        if ndim == 1:
            boundary_entries = [key for key in sim_config_src if key.startswith('bins.boundaries')]
            if not boundary_entries:
                raise ConfigError('no bin boundaries provided')
            elif len(boundary_entries) > 1:
                raise ConfigError('more than one bin boundary set provided')
            else:
                boundary_entry = boundary_entries[0]
            
            if boundary_entry in ('bins.boundaries', 'bins.boundaries_0'):
                boundaries = [float(bound) for bound in sim_config_src.get_list(boundary_entry)]
            elif boundary_entry in ('bins.boundaries_expr', 'bins.boundaries_expr_0'):
                boundaries = eval(sim_config_src[boundary_entry])
            else:
                raise ConfigError('invalid bin boundary specification')
                
            bin_limits.append(numpy.array(boundaries))
        else:
            reIsBoundary = re.compile('bins\.boundaries(_expr)?_(\d+)')
            bin_limits = [None] * ndim
            
            boundary_entries = [entry for entry in sim_config_src if entry.startswith('bins.boundaries')]
            for boundary_entry in boundary_entries:
                m = reIsBoundary.match(boundary_entry)
                if not m:
                    raise ConfigError('invalid bin boundary specification')
                else:
                    idim = int(m.group(2))
                    if m.group(1):
                        boundaries = eval(sim_config_src[boundary_entry])
                    else:
                        boundaries = [float(lim) for lim in sim_config_src.get_list(boundary_entry)]
                    bin_limits[idim] = numpy.array(boundaries)
                    
            if None in bin_limits:
                raise ConfigError('missing bin boundaries for at least one dimension')
        
        sim_config['bins.boundaries'] = self.bin_boundaries = bin_limits
        sim_config['bins.particles_per_bin'] = self.particles_per_bin = sim_config_src.get_int('bins.particles_per_bin')
    
    def make_bins(self):
        return BinArray(boundaries = self.bin_boundaries,
                        ideal_num = self.particles_per_bin)
        # non-time-dependent bin info stored here
