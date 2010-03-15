import os, sys
from optparse import OptionParser
import numpy, h5py
from wemd.analysis.transitions import OneDimTransitionEventFinder
from wemd.util.config_dict import ConfigDict, ConfigError

transcfg = ConfigDict()

parser = OptionParser(usage='%prog [ANLCONF]',
                      description='Perform transition analysis')
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

try:
    translog = transcfg.get_file_object('output.transition_log', mode='w')
except KeyError:
    translog = None
    
regions = []
for (irr, rname) in enumerate(region_names):
    regions.append((rname, (region_edges[irr], region_edges[irr+1])))

source_filename = transcfg['data.source']
source_node = transcfg['data.source_node']

h5file = h5py.File(source_filename)
data = h5file[source_node][...]
sys.stdout.write('source data is %s\n' % 'x'.join([str(x) for x in data.shape]))

if data.ndim > 1 and transcfg.get_bool('data.contains_time', True):
    pcoord = data[:, 1:]
    sys.stdout.write('using time data from source file\n')
    timestep = data[1,0] - data[0,0]
else:
    pcoord = data
    timestep = transcfg.get_float('data.timestep', 1.0)
sys.stdout.write('using dt=%g\n' % timestep)

event_durations = {}
fpts = {}
for irr1 in xrange(0, len(regions)):
    for irr2 in xrange(0, len(regions)):
        if abs(irr1-irr2) > 1:
            event_durations[irr1,irr2] = numpy.empty((0,2), numpy.float64)
            fpts[irr1, irr2] = numpy.empty((0,2), numpy.float64)

trans_finder = OneDimTransitionEventFinder(regions,
                                           pcoord,
                                           dt = timestep,
                                           transition_log = translog)

trans_finder.identify_regions()
trans_finder.identify_transitions()

for ((region1, region2), ed_array) in trans_finder.event_durations.iteritems():
    region1_name = regions[region1][0]
    region2_name = regions[region2][0]
    if ed_array.shape[0] == 0:
        sys.stdout.write('No %s->%s transitions observed\n')
        continue
    
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
    
    try:
        fpt_array = trans_finder.fpts[region1, region2]
    except KeyError:
        sys.stdout.write('No %s->%s FPTs observed\n'
                         % (region1_name, region2_name))
    else:
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
