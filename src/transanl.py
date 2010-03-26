from __future__ import division

import os, sys
from optparse import OptionParser
import numpy, h5py
from wemd.analysis.transitions import TransitionEventAccumulator, OneDimRegionSet
from wemd.util.config_dict import ConfigDict, ConfigError

transcfg = ConfigDict()

parser = OptionParser(usage='%prog [ANLCONF]',
                      description='Perform transition analysis')
parser.add_option('-c', '--chunk-size', type='int', dest='chunk_size',
                  help='retrieve CHUNK_SIZE coordinates from disk at a time')
(opts, args) = parser.parse_args()
if args:
    transcfg.read_config_file(args[0])
else:
    transcfg.read_config_file('analysis.cfg')

transcfg.require_all(['regions.edges', 'regions.names',
                      'data.source', 'data.source_node'])

region_names = transcfg.get_list('regions.names')
region_edges = transcfg.get_list('regions.edges', type=float)
if len(region_edges) != len(region_names) + 1:
    sys.stderr.write('region names and boundaries do not match')
    sys.exit(1)

region_boundaries = []
for (irr, rname) in enumerate(region_names):
    region_boundaries.append((region_edges[irr], region_edges[irr+1]))

regions = OneDimRegionSet(region_names, region_boundaries)

try:
    translog = transcfg.get_file_object('output.transition_log', mode='w')
except KeyError:
    translog = None
    

source_filename = transcfg['data.source']
source_node = transcfg['data.source_node']

h5file = h5py.File(source_filename)
data = h5file[source_node]
sys.stdout.write('source data is %s\n' % 'x'.join([str(x) for x in data.shape]))

trans_finder = TransitionEventAccumulator(regions, data_overlaps=False,
                                          transition_log = translog)

chunk_size = opts.chunk_size
if not chunk_size:
    if data.chunks is not None and data.chunks[0] >= 10:
        chunk_size = data.chunks[0]
    else:
        chunk_size = 1000
weights = numpy.ones((chunk_size,), numpy.float64)

sys.stdout.write('using chunks of %d\n' % chunk_size)
if len(data.shape) > 1 and transcfg.get_bool('data.contains_time', True):
    contains_time = True
    timestep = None
else:
    contains_time = False
    timestep = transcfg.get_float('data.timestep', 1.0)
    sys.stdout.write('using dt=%g\n' % timestep)
    
sys.stdout.write('identifying transitions...\n')
for istart in xrange(0, data.shape[0], chunk_size):
    if istart+chunk_size > data.shape[0]:
        iend = data.shape[0]
    else:
        iend = istart+chunk_size
    
    if contains_time:
        if timestep is None:
            data_chunk = data[istart:iend,:]
            timestep = data_chunk[1,0] - data_chunk[0,0]
            sys.stdout.write('(using timestep = %g)\n' % timestep)
            pc_chunk = data_chunk[:,1:]
            del data_chunk
        else:
            pc_chunk = data[istart:iend,1:]
    else:
        pc_chunk = data[istart:iend]
    
    sys.stdout.write('%d/%d (%.1f%%)\n' % (istart, data.shape[0], 100*istart/data.shape[0]))
    trans_finder.timestep = timestep
    trans_finder.identify_transitions(pc_chunk,weights)
    #sys.stdout.write('%s\n\n' % trans_finder.event_counts)
    
sys.stdout.write('event count (row->column, states %s)\n' % ', '.join(regions.names))
sys.stdout.write('%s\n' % trans_finder.event_counts)
for ((ir1, ir2), ed_list) in trans_finder.eds.iteritems():
    region1_name = regions.names[ir1]
    region2_name = regions.names[ir2]
    if len(ed_list) == 0:
        sys.stdout.write('No %s->%s transitions observed\n'
                         % (region1_name, region2_name))
        continue
    
    ed_array = numpy.array(ed_list, numpy.float64)
    ed_array[:,0] *= timestep
    
    sys.stdout.write('\nStatistics for %s->%s:\n'
                     % (region1_name, region2_name))
    sys.stdout.write('Number of events: %d\n' % ed_array.shape[0])
    sys.stdout.write('ED average:       %g\n' % ed_array[:,0].mean())
    ed_stdev = ed_array[:,0].std()
    ed_sem = ed_stdev / ed_array.shape[0] ** 0.5
    sys.stdout.write('ED st. dev.:      %g\n' % ed_stdev)
    sys.stdout.write('ED S.E.M.:        %g\n' % ed_sem)
    sys.stdout.write('ED min:           %g\n' % ed_array[:,0].min())
    sys.stdout.write('ED median:        %g\n' % ed_array[ed_array.shape[0]/2,0])
    sys.stdout.write('ED max:           %g\n' % ed_array[:,0].max())
    
    ed_file = open('ed_%s_%s.txt' % (region1_name, region2_name), 'wt')
    for irow in xrange(0, ed_array.shape[0]):
        ed_file.write('%20.16g    %20.16g\n' % tuple(ed_array[irow,0:2]))
    ed_file.close()
    
    
    fpt_list = trans_finder.fpts.get((ir1, ir2), None)
    if not fpt_list:
        sys.stdout.write('No %s->%s FPTs observed\n'
                         % (region1_name, region2_name))
    else:
        sys.stdout.write('%d %s->%s FPTs observed\n'
                         % (len(fpt_list), region1_name, region2_name))
        fpt_array = numpy.array(fpt_list, numpy.float64)
        fpt_array[:,0] *= timestep
        
        fpt_mean = fpt_array[:,0].mean()
        fpt_stdev = fpt_array[:,0].std()
        fpt_sem = fpt_stdev / fpt_array.shape[0] ** 0.5

        sys.stdout.write('FPT average:      %g\n' % fpt_mean)
        sys.stdout.write('FPT st. dev.:     %g\n' % fpt_stdev)
        sys.stdout.write('FPT S.E.M.:       %g\n' % fpt_sem)
        sys.stdout.write('FPT min:          %g\n' % fpt_array[:,0].min())
        sys.stdout.write('FPT median:       %g\n' % fpt_array[fpt_array.shape[0]/2,0])
        sys.stdout.write('FPT max:          %g\n' % fpt_array[:,0].max())
        
        rate_mean = 1.0/fpt_mean
        rate_stdev = fpt_stdev / fpt_mean**2
        rate_sem = rate_stdev / fpt_array.shape[0] ** 0.5
        
        sys.stdout.write('rate average:     %g\n' % rate_mean)
        sys.stdout.write('rate st. dev.:    %g\n' % rate_stdev)
        sys.stdout.write('rate S.E.M.:      %g\n' % rate_sem)
        
        fpt_file = open('fpt_%s_%s.txt' % (region1_name, region2_name), 'wt')
        for irow in xrange(0, fpt_array.shape[0]):
            fpt_file.write('%20.16g    %20.16g\n' % tuple(fpt_array[irow,0:2]))
        fpt_file.close()
