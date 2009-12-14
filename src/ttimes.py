import os, sys
from optparse import OptionParser
from configparse import SafeConfigParser
import numpy
from wemd.util.binfiles import UniformTimeData
from itertools import izip

#region_edges = [0.0, 0.4, 0.65, 0.85, 1.0, float('inf')]
#region_names = ['B', 'T3', 'T2', 'T1', 'A']
region_edges = [0.0, 0.4, 1.0, float('inf')]
region_names = ['B', 'T', 'A']

tb_by_type = {(0,2): [], (2,0): []}
fpt_by_type = {(0,2): [], (2,0): []}

assert len(region_edges) == len(region_names) + 1
nreg = len(region_names)

dist_utd = UniformTimeData()
dist = dist_utd.mmap_array_from(open('dist.dat', 'rb'))
dt = dist_utd.dt

#td = numpy.load('dist.pkl')
#dist = numpy.zeros(td.shape, numpy.float)
#dist[:,0] = td[:,1] * 0.1
#dist[:,1] += 1.0
#dt = td[1,0] - td[0,0]
#print dist[0:4]

completion_times = numpy.zeros((nreg,nreg), numpy.int64) - 1
event_count = numpy.zeros((nreg,nreg), numpy.int64)

for ir in xrange(0, nreg):
    q = dist[0,0]
    if region_edges[ir] <= dist[0,0] < region_edges[ir+1]:
        last_iregion = ir
        break
else:
    raise ValueError('particle is outside specified coordinate space')
    
for it in xrange(1, dist.shape[0]):
    q = dist[it,0]
    for ir in xrange(0, nreg):
        lb = region_edges[ir]
        ub = region_edges[ir+1]
        
        if lb <= q < ub:
            iregion = ir
            break
    else:
        raise ValueError('particle is outside specified coordinate space')
    
    if iregion != last_iregion:        
        # Crossing event...do all the real work here
    
        completion_times[last_iregion, iregion] = it
        event_count[last_iregion, iregion] += 1
                
        #print "%d (%s->%s)" % (it, region_names[last_iregion], region_names[iregion])
        
        for irr in xrange(1, nreg):
            if iregion > last_iregion:
                trans_iregion = last_iregion - irr
                tb_iregion = trans_iregion + 1
            else:
                trans_iregion = last_iregion + irr
                tb_iregion = trans_iregion - 1
                
            if 0 <= trans_iregion < nreg and trans_iregion != iregion:
                if completion_times[iregion, last_iregion] > completion_times[trans_iregion, last_iregion]:
                    # Fluctuation between regions
                    #print "  cannot be a %s->%s transition" % (region_names[trans_iregion], region_names[iregion])
                    pass
                else:
                    try:
                        type_tbs = tb_by_type[trans_iregion, iregion]
                        type_fpts = fpt_by_type[trans_iregion, iregion]
                    except KeyError:
                        pass
                    else:
                        #print "  %s->%s transition complete" % (region_names[trans_iregion], region_names[iregion])
                
                        # First passage time:  time since last transition ending in starting state
                        if completion_times[iregion, trans_iregion] >= 0:
                            fpt = it-completion_times[iregion, trans_iregion]
                            type_fpts.append(fpt) 
                            #print "  %s->%s first passage time: %d" \
                            #      % (region_names[trans_iregion], region_names[iregion], it-completion_times[iregion, trans_iregion])
                    
                        # Event duration: time since last transition originating from starting state
                        # trans_iregion -> trans_iregion +/- 1
                        if completion_times[trans_iregion, tb_iregion] >= 0:
                            tb = it-completion_times[trans_iregion, tb_iregion]
                            type_tbs.append(tb)
                            #print "  %s->%s event duration: %d" \
                            #      % (region_names[trans_iregion], region_names[iregion], 
                            #         it-completion_times[trans_iregion, tb_iregion])
                    # Update completion times matrix
                    completion_times[trans_iregion, iregion] = it
                    event_count[trans_iregion, iregion] += 1
        last_iregion = iregion
    #if it == 1000000: break
del dist
print "event counts:"
print event_count

tb_AB = numpy.array(tb_by_type[2,0], numpy.float)
tb_AB = tb_AB * dt
fpt_AB = numpy.array(fpt_by_type[2,0], numpy.float)
fpt_AB = fpt_AB * dt
print "average A->B duration:", tb_AB.mean()
print "stdev   A->B duration:", tb_AB.std()
print "SEM     A->B duration:", tb_AB.std()/(tb_AB.size**0.5)
print
print "average A->B FPT:", fpt_AB.mean()
print "stdev   A->B FPT:", fpt_AB.std()
print "SEM     A->B FPT:", fpt_AB.std()/(fpt_AB.size**0.5)
