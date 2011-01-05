import re

import logging
log = logging.getLogger(__name__)

import numpy
from we_driver import WEDriver
from wemd.core.binarrays import Bin, BinArray
from wemd.core import ConfigError

class FixedBinWEDriver(WEDriver):
    def __init__(self, sim_manager):
        super(FixedBinWEDriver,self).__init__(sim_manager)
        
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

        #read in the recycling regions
        #[ [lower bound, uppder bound], ...]
        target_pcoords = []
        target_entries = [key for key in sim_config_src if key.startswith('bins.target_pcoord')]
        for target_entry in target_entries:
            target_pcoords.append( eval(sim_config_src.get(target_entry)) )

        #check to make sure they are correct
        for i in xrange(0, len(target_pcoords)):

            if len(target_pcoords[i]) != ndim:
                raise ConfigError("Invalid target pcoord specified, check dimensions %r" % (target_pcoords[i]))

            for idim in xrange(0, len(target_pcoords[i])):
                if len(target_pcoords[i][idim]) != 2:
                    raise ConfigError("Invalid target pcoord specified, [ [lower bound, upper bound], ...] %r" %(target_pcoords[i]))

        if not target_pcoords:
            raise ConfigError("No recycling targets specified")
            
        #read in the source regions
        #source_pcoords['Region Name'] is a dictionary with the keys 'pcoord', 'region', and 'weight'
        source_pcoords = {}
        source_entries = [key for key in sim_config_src if key.startswith('bins.source_pcoord_')]
        reIsSourcePcoord = re.compile('bins\.source_pcoord_(.+)_(.+)')

        for source_entry in source_entries:
            m = reIsSourcePcoord.match(source_entry)
            if not m:
                raise ConfigError('Invalid source pcoord specified')
            else:
                if m.group(2) not in source_pcoords.keys():
                    source_pcoords[m.group(2)] = {} 

                if m.group(1) == 'weight':
                    source_pcoords[m.group(2)][m.group(1)] = sim_config_src.get_float(source_entry)
                elif m.group(1) == 'pcoord':
                    pcoord_vals = [float(x) for x in sim_config_src.get_list(source_entry)]
                    source_pcoords[m.group(2)][m.group(1)] = numpy.array( pcoord_vals, dtype=numpy.float64 )
                elif m.group(1) == 'region':
                    source_pcoords[m.group(2)][m.group(1)] = eval(sim_config_src.get(source_entry))
                else:
                    raise ConfigError('Invalid source pcoord specified')     
                               
        #check to make sure they are correct
        for source_name in source_pcoords:

            for key in ('weight', 'pcoord', 'region'):
                if key not in source_pcoords[source_name]:
                    raise ConfigError("Invalid source pcoord specified -- missing key %s -- source_pcoord_%s" % (key, source_name))

            if (source_pcoords[source_name]['weight'] < 0) or (source_pcoords[source_name]['weight'] > 1):
                raise ConfigError("Invalid initial weight -- source_pcoord_%s" % (source_name))
            
            if len(source_pcoords[source_name]['pcoord']) != ndim:
                raise ConfigError("Invalid source pcoord specified, check pcoord dimensions -- source_pcoord_%s" % (source_name))
            
            if len(source_pcoords[source_name]['region']) != ndim:
                raise ConfigError("Invalid source region specified, check region dimensions -- source_pcoord_%s" % (source_name))
            
            for idim in xrange(0, len(source_pcoords[source_name]['region'])):
                if len(source_pcoords[source_name]['region'][idim]) != 2:
                    raise ConfigError("Invalid source pcoord region specified "+\
                                      "[ [lower bound, upper bound], ...] -- source_pcoord_%s" %(source_name))
        
        sim_config['bin.target_pcoords'] = self.target_pcoords = numpy.array(target_pcoords,dtype=numpy.float64)
        sim_config['bin.source_pcoords'] = self.source_pcoords = source_pcoords
        sim_config['bins.boundaries'] = self.bin_boundaries = bin_limits
        sim_config['bins.particles_per_bin'] = self.particles_per_bin = sim_config_src.get_int('bins.particles_per_bin')       
    
    def make_bins(self):
        return BinArray(boundaries = self.bin_boundaries,
                        ideal_num = self.particles_per_bin)
        # non-time-dependent bin info stored here
