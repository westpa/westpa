from __future__ import print_function, division
import os, sys, argparse, math
import numpy
import wemd
from wemd.pcoords import RectilinearRegionSet

region_set = RectilinearRegionSet([[0, 0.4, 1.0, float('inf')]])
all_bins = region_set.get_all_bins()
nbins = len(all_bins)
all_indices = set(xrange(0,nbins))
labels = [bin.label for bin in all_bins]

last_crossing   = numpy.zeros((nbins,nbins), numpy.uintp)
last_entry      = numpy.zeros((nbins,), numpy.uintp)
last_exit       = numpy.zeros((nbins,), numpy.uintp)
last_completion = numpy.zeros((nbins,nbins), numpy.uintp)

n_crossings     = numpy.zeros((nbins,nbins), numpy.uintp)
n_completions   = numpy.zeros((nbins,nbins), numpy.uintp)

ed_weight_acc   = numpy.zeros((nbins,nbins), numpy.float64)
ed_acc          = numpy.zeros((nbins,nbins), numpy.float64)
ed_sq_acc       = numpy.zeros((nbins,nbins), numpy.float64)
fpt_weight_acc  = numpy.zeros((nbins,nbins), numpy.float64)
fpt_acc         = numpy.zeros((nbins,nbins), numpy.float64)
fpt_sq_acc      = numpy.zeros((nbins,nbins), numpy.float64)

dist = numpy.load('dist.npy', 'r')

tfile = open('transitions.txt', 'wt')
edfile = open('durations.txt', 'wt')
fptfile = open('fpts.txt', 'wt')
dwellfile = open('dwells.txt', 'wt')

last_index = region_set.map_to_indices([dist[0]])[0]
CHUNK=100000
for ii in xrange(1, len(dist), CHUNK):
    print('{:d}/{:d}'.format(ii, len(dist)))
    dist_chunk = numpy.expand_dims(dist[ii:ii+CHUNK],1)
    indices_chunk = region_set.map_to_indices(dist_chunk)
    
    for ij in xrange(0, CHUNK, 1):
        i = ii+ij
        weight = 1
        current_index = indices_chunk[ij]
    
        if current_index != last_index:
            clabel = '{:d}->{:d}'.format(long(last_index), long(current_index))
            tfile.write('{time:20d}    {label:20s}    {weight:20.14e}\n'.format(time=i, label=clabel, weight=weight))
            n_crossings[last_index,current_index] += 1
                
            #not_current = list(sorted(all_indices - {current_index}))
            #not_last    = list(sorted(all_indices - {last_index}))
                    
            # event duration X->Y: (time of this entry into Y) - (time of last exit from X)
            for istart in all_indices:
                if istart == current_index: continue    # renewal
                if istart == last_index: continue       # just a crossing
                tlabel = '{:d}->{:d}'.format(long(istart), long(current_index))
                
                # initial region has been visited since last completion of initial->current;
                # this indicates a transition from initial->current
                if last_entry[istart] > last_completion[istart, current_index]:
                    # We can always compute event durations for nonadjacent states
                    ed = i - last_exit[istart]
                    fed = float(ed)
                    ed_weight_acc[istart,current_index] += weight
                    ed_acc[istart,current_index] += fed*weight
                    ed_sq_acc[istart,current_index] += fed*fed*weight
                    edfile.write('{time:20d}    {label:20s}    {ed:20d}    {weight:20.14e}\n'
                                 .format(time=i, label=tlabel, ed=long(ed), weight=weight))
                    
                    # If we have seen a current->initial transition, then we can compute an FPT
                    # for initial->current
                    if last_completion[current_index, istart] > 0:
                        fpt = i - last_completion[current_index, istart]
                        ffpt = float(fpt)
                        fpt_weight_acc[istart,current_index] += weight
                        fpt_acc[istart,current_index] += ffpt*weight
                        fpt_sq_acc[istart,current_index] += ffpt*ffpt*weight
                        fptfile.write('{time:20d}    {label:20s}    {fpt:20d}    {weight:20.14e}\n'
                                     .format(time=i, label=tlabel, fpt=long(fpt), weight=weight))
                        
                    last_completion[istart,current_index] = i
                    n_completions[istart,current_index] += 1
                
            last_exit[last_index] = i
            last_entry[current_index] = i
            last_crossing[last_index, current_index] = i
                
        last_index = current_index

avg_ed = ed_acc / ed_weight_acc
stdev_ed = (ed_sq_acc/ed_weight_acc - avg_ed*avg_ed)**0.5
avg_fpt = fpt_acc / fpt_weight_acc
stdev_fpt = (fpt_sq_acc/fpt_weight_acc - avg_fpt*avg_fpt)**0.5

print('Number of crossings:')
print(n_crossings)

print('Number of completed transitions:')
print(n_completions)

print('Average event durations:')
print(avg_ed)
print('Standard deviation of event durations:')
print(stdev_ed)

print('Average first passage times:')
print(avg_fpt)
print('Standard deviation of first passage times:')
print(stdev_fpt)
