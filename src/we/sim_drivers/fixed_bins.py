import numpy
from we.core.we_sim import WESim
from we.core.binarrays import Bin, BinArray
from we.core import ConfigError

class FixedBinWESim(WESim):
    def configure_bin_params(self, config):
        bintype = config['bins.type']
        assert(bintype == 'fixed')

        ndim = config.get_int('bins.ndim')
        bin_limits = []
        
        if ndim == 1:
            boundary_entries = [key for key in config if key.startswith('bins.boundaries')]
            if not boundary_entries:
                raise ConfigError('no bin boundaries provided')
            elif len(boundary_entries) > 1:
                raise ConfigError('more than one bin boundary set provided')
            else:
                boundary_entry = boundary_entries[0]
            
            if boundary_entry in ('bins.boundaries', 'bins.boundaries_0'):
                boundaries = [float(bound) for bound in config.get_list(boundary_entry)]
            elif boundary_entry in ('bins.boundaries_expr', 'bins.boundaries_expr_0'):
                boundaries = eval(config[boundary_entry])
            else:
                raise ConfigError('invalid bin boundary specification')
                
            bin_limits.append(numpy.array(boundaries))
        else:
            reIsBoundary = re.compile('bins\.boundaries(_expr)?_(\d+)')
            bin_limits = [None] * ndim
            
            boundary_entries = [entry for entry in config if entry.startswith('bins.boundaries')]
            for boundary_entry in boundary_entries:
                m = reIsBoundary.match(boundary_entry)
                if not m:
                    raise ConfigError('invalid bin boundary specification')
                else:
                    ndim = int(m.group(2))
                    if m.group(1):
                        boundaries = eval(config[boundary_entry])
                    else:
                        boundaries = [float(lim) for lim in config.get_list(boundary_entry)]
                    boundary_entries[ndim] = numpy.array(boundaries)
                    
            if None in bin_limits:
                raise ConfigError('missing bin boundaries for at least one dimension')
        
        self.data_manager.config['bins.boundaries'] =  bin_limits
        self.data_manager.config['bins.particles_per_bin'] \
            = config.get_int('bins.particles_per_bin')
        self.data_manager.config['bins.split_threshold'] \
            = config.get_float('bins.split_threshold', 2.0)
        self.data_manager.config['bins.merge_threshold_min'] \
            = config.get_float('bins.merge_threshold_min', 0.5)
        self.data_manager.config['bins.merge_threshold_max'] \
            = config.get_float('bins.merge_threshold_max', 1.5)        
        
        self.bins = self.make_bins()
    
    def make_bins(self):
        config = self.data_manager.config
        return BinArray(boundaries = config['bins.boundaries'],
                        ideal_num = config['bins.particles_per_bin'],
                        split_threshold = config['bins.split_threshold'],
                        merge_threshold_min = config['bins.merge_threshold_min'],
                        merge_threshold_max = config['bins.merge_threshold_max'])
        # non-time-dependent bin info stored here
