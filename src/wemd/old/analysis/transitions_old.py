import numpy
from itertools import izip

import logging
log = logging.getLogger(__name__)

__docformat__ = "restructuredtext en"

class TransitionEventFinder(object):
    """A class for identifying transition events within progress coordinate
    space. Subclasses are responsible for providing a method for determining
    which region the particle is in for each (time, progress coordinate)
    pair.  Currently limited to analysis of transitions within sub-regions
    which can be mapped to a continuous one-dimensional space.
    
    :IVariables:
      regions : sequence
        A sequence of ``(region name, region boundary)`` tuples. In most cases,
        an index into this sequence is used to identify regions.
      pcoords : array
        The 2-D array of progress coordinate data used by the transition
        search routine. The first dimension indexes time.
      t0 : float
        The initial time of ``pcoords``. Defaults to 0.0.
      dt : float
        The time interval of ``pcoords``. Defaults to 1.0.
      time_threshold : float
        Only transitions ending after this time will be recorded. 
        Defaults to ``t0``. 
      weights : array
        The trajectory weight as a function of time.
      traj_id : object (printable) or ``None``
        When logging individual transitions and ``traj_id`` is not ``None``, 
        each identified transition will be labeled with this identifier when 
        written to the transition log.
      transition_log : object (file-like) or ``None``
        When not ``None``, individual transitions are written to this file.
      pcoord_regions : array
        A 1-D array of indices into the ``regions`` sequence identifying what
        region the progress coordinate indicates as a function of time.
        Available after calling `identify_regions()`
      event_counts : array
        A 2-D array of indices into the ``regions`` sequence.
        ``event_counts[a,b]`` gives the number of ``a -> b`` transitions found.
        Available after calling `identify_transitions()`.
      event_durations : dict
        A mapping of a region index to an array of ``(duration, weight)``. 
        Available after calling `identify_transitions()`.
      fpts : dict
        A mapping of a region index to an array of ``(fpt, weight)``.
        Available after calling `identify_transitions()`.
        
    """
    def __init__(self, 
                 regions,
                 pcoords,
                 t0 = 0.0,
                 dt = 1.0,
                 weights = None, 
                 traj_id = 0,
                 transition_log = None,
                 save_durations = None,
                 save_fpts = None):
        """Set up a transition finder.
        
        :Parameters:
          regions : sequence
            A sequence of ``(region name, region boundary)`` tuples.
          pcoords : array
            The 2-D array of progress coordinate data used by the transition
            search routine. The first dimension indexes time.
          t0 : float
            The initial time of ``pcoords``. Defaults to 0.0.
          dt : float
            The time interval of ``pcoords``. Defaults to 1.0.
          time_threshold : float
            Only transitions ending after this time will be recorded. 
            Defaults to `t0`. 
          weights : array
            The trajectory weight as a function of time.
          traj_id : object (printable) or ``None``
            When logging individual transitions and ``traj_id`` is not ``None``, 
            each identified transition will be labeled with this identifier 
            when written to the transition log.
          transition_log : object (file-like) or ``None``
            When not ``None``, individual transitions are written to this file.
          save_durations : dict
            A set of ``(int,int)`` tuples indicating the transition 
            types for which to store event duration data. If not provided
            (or ``None`` is provided), event duration data is stored for every 
            possible transition type. If an empty dictionary is provided,
            no event duration time data will be stored. 
          save_fpts : dict
            A set of ``(int,int)`` tuples indicating the transition 
            types for which to store first passage time data. If not provided
            (or `None` is provided), FPT data is stored for every 
            possible transition type. If an empty dictionary is provided,
            no FPT will be stored. 
              
        """
        self.regions = regions
        self.pcoords = pcoords
        self.t0 = t0
        self.dt = dt
        self.pcoord_regions = None
        self.event_counts = None
        self.event_durations = {}
        self.fpts = {}
        
        if weights is None:
            self.weights = numpy.ones((pcoords.shape[0],), dtype=numpy.float64)
        else:
            self.weights = weights
        
        self.traj_id = traj_id
        self.transition_log = transition_log
        
        self.save_durations = set()
        if save_durations is None:
            for irr1 in xrange(0, len(self.regions)):
                for irr2 in xrange(0, len(self.regions)):
                    if abs(irr1-irr2) > 1:
                        self.save_durations.add((irr1,irr2))
        else:
            self.save_durations = save_durations
        
        self.save_fpts = set()    
        if save_fpts is None:
            for irr1 in xrange(0, len(self.regions)):
                for irr2 in xrange(0, len(self.regions)):
                    if abs(irr1-irr2) > 1:
                        self.save_fpts.add((irr1,irr2))
        else:
            self.save_fpts = save_fpts
    
    def identify_regions(self):
        raise NotImplementedError
    
    def identify_transitions(self):
        nreg = len(self.regions)
        pcoords = self.pcoords
        weights = self.weights
        regions = self.regions
        
        if self.pcoord_regions is None:
            self.identify_regions()
        pcoord_regions = self.pcoord_regions
        
        completion_times = numpy.zeros((nreg,nreg), numpy.int64) - 1
        event_counts = self.event_counts = numpy.zeros((nreg,nreg), 
                                                       numpy.uint64)
        transition_log = self.transition_log
        t0 = self.t0
        dt = self.dt
        traj_id = self.traj_id
        
        event_durations = {}
        fpts = {}
                
        for (irr1,irr2) in self.save_durations:
            event_durations[irr1,irr2] = []
        for (irr1,irr2) in self.save_fpts:
            fpts[irr1,irr2] = []
            
        for it in xrange(1, pcoords.shape[0]):
            last_iregion = pcoord_regions[it-1]
            iregion = pcoord_regions[it]
                        
            if iregion != last_iregion:
                # A crossing has happened
                completion_times[last_iregion, iregion] = it
                event_counts[last_iregion, iregion] += 1
                
                t = self.t0 + self.dt * it
                
                if transition_log:
                    self.transition_log.write('%s    %g    %g    %s->%s\n'
                                              % (self.traj_id or '', 
                                                 t, weights[it],
                                                 regions[last_iregion][0],
                                                 regions[iregion][0]))
                
                for irr in xrange(1, nreg):
                    # Index adjustment.  What fun.
                    if iregion > last_iregion:
                        # Forward motion; a transition began in the region
                        # "before" last_region.
                        trans_iregion = last_iregion - irr
                        tb_iregion = trans_iregion + 1
                    else:
                        # Backward motion; a transition began in the region
                        # "after" last_region
                        trans_iregion = last_iregion + irr
                        tb_iregion = trans_iregion - 1
                    if 0 <= trans_iregion < nreg and trans_iregion != iregion:
                        if completion_times[iregion, last_iregion] > completion_times[trans_iregion, last_iregion]:
                            # Fluctuation between regions
                            pass
                        else:                            
                            # First passage time: time since last transition
                            # ending in starting state
                            
                            # Test to see if the destination state has been
                            # visited at least once
                            if completion_times[iregion, trans_iregion] >= 0:
                                fpt = it-completion_times[iregion, trans_iregion]
                                log.debug('found %s->%s transition (FPT %d)' \
                                       % (self.regions[trans_iregion][0], 
                                          self.regions[iregion][0],
                                          fpt))
                                

                                log.debug('recording %s->%s FPT'
                                          % (self.regions[trans_iregion][0],
                                             self.regions[iregion][0]))
                                try:
                                    fpts[int(trans_iregion), 
                                         int(iregion)].append((fpt, weights[it],t0+it*dt))
                                except KeyError:
                                    pass
                                    
                            # Event duration: time since last transition
                            # originating from a neighboring state
                            
                            # If the neighboring state has been visited...
                            if completion_times[trans_iregion, tb_iregion] >= 0:
                                tb = it - completion_times[trans_iregion, tb_iregion]-1
                                log.debug('found %s->%s transition (duration %d)' \
                                       % (self.regions[trans_iregion][0], 
                                          self.regions[iregion][0],
                                          tb))
                                
                                log.debug('recording %s->%s transition time'
                                          % (self.regions[trans_iregion][0],
                                             self.regions[iregion][0]))
                                try:
                                    event_durations[int(trans_iregion),
                                                    int(iregion)].append((tb,weights[it],t0+it*dt))
                                except KeyError:
                                    pass
                                
                            # Update off-diagonal elements of the completion
                            # times matrix
                            completion_times[trans_iregion, iregion] = it
                            event_counts[trans_iregion, iregion] += 1
                last_iregion = iregion
                            
        for (k,v) in fpts.iteritems():
            if len(v):
                self.fpts[k] = numpy.array(v, numpy.float64)
                self.fpts[k][:,0] *= self.dt
        for (k,v) in event_durations.iteritems():
            if len(v):
                self.event_durations[k] = numpy.array(v, numpy.float64)
                self.event_durations[k][:,0] *= self.dt
                    
class OneDimTransitionEventFinder(TransitionEventFinder):
    def identify_regions(self):
        pcoords = self.pcoords
        pcoord_regions = self.pcoord_regions = numpy.empty((pcoords.shape[0],),
                                                           int)
        regions = self.regions
        
        for it in xrange(0, pcoords.shape[0]):
            q = pcoords[it,0]
            for irr in xrange(0, len(regions)):
                lb, ub = regions[irr][1]
                if lb <= q < ub:
                    pcoord_regions[it] = irr
                    break

class AltOneDimTransitionEventFinder(TransitionEventFinder):
    def identify_regions(self):
        return
    
    def identify_transitions(self):
        nreg = len(self.regions)
        pcoords = self.pcoords
        weights = self.weights
        regions = self.regions 
                
        completion_times = numpy.zeros((nreg,nreg), numpy.int64) - 1
        event_counts = self.event_counts = numpy.zeros((nreg,nreg), 
                                                       numpy.uint64)
        transition_log = self.transition_log
        t0 = self.t0
        dt = self.dt
        traj_id = self.traj_id
        
        save_durations = self.save_durations
        save_fpts = self.save_fpts
        event_durations = {}
        fpts = {}
        
        for (irr1,irr2) in save_durations:
            event_durations[irr1,irr2] = []
        for (irr1,irr2) in save_fpts:
            fpts[irr1,irr2] = []
        
        last_forward_trans_entry = [0, 1]
        last_backward_trans_entry = [0, 1]
        last_forward_trans_exit = [0, 1]
        last_backward_trans_exit = [0, 1]
        last_forward_completion = [0, 1]
        last_backward_completion = [0, 1]
        
        q = pcoords[0,0]
        for irr in xrange(0, nreg):
            lb, ub = regions[irr][1]
            if lb <= q < ub:
                last_iregion = irr
                break
        
        for it in xrange(1, pcoords.shape[0]):
            q = pcoords[it, 0]
            for irr in xrange(0, nreg):
                lb, ub = regions[irr][1]
                if lb <= q < ub:
                    iregion = irr
                    break
                
            if iregion != last_iregion:
                # A crossing has occurred
                
                t = self.t0 + self.dt * it
                
                if transition_log:
                    transition_log.write('%-12s    %21.16g    %21.16g    %s->%s\n'
                                          % (self.traj_id or '', 
                                             t, weights[it],
                                             regions[last_iregion][0],
                                             regions[iregion][0]))
                                
                if iregion == nreg-1: 
                    # "forward" crossing into endpoint
                    if last_forward_completion[0] < last_forward_trans_entry[0]:
                        if (0, iregion) in save_durations and last_forward_trans_entry[0] > 0:
                            event_durations[0, iregion].append((it - last_forward_trans_entry[0] - 1, weights[it]))
                        if (0, iregion) in save_fpts and last_backward_completion[0] > 0:
                            fpts[0, iregion].append((it-last_backward_completion[0], weights[it]))
                        last_forward_completion = (it, weights[it])
                elif iregion == 0: 
                    # "backward" crossing into endpoint
                    if last_backward_completion[0] < last_backward_trans_entry[0]:
                        if (nreg-1, 0) in save_durations and last_backward_trans_entry[0] > 0:
                            event_durations[nreg-1, 0].append((it - last_backward_trans_entry[0] - 1, weights[it]))
                        if (nreg-1, 0) in save_fpts and last_forward_completion[0] > 0:
                            fpts[nreg-1,0].append((it - last_forward_completion[0], weights[it]))

                        last_backward_completion = (it, weights[it])
                elif last_iregion == 0:
                    # "forward" crossing into transition region
                    last_forward_trans_entry = (it, weights[it])
                elif last_iregion == nreg-1:
                    # "backward" crossing into transition region
                    last_backward_trans_entry = (it, weights[it])
                
            last_iregion = iregion
            
        for (k,v) in fpts.iteritems():
            if len(v):
                self.fpts[k] = numpy.array(v, numpy.float64)
                self.fpts[k][:,0] *= self.dt + self.t0
        for (k,v) in event_durations.iteritems():
            if len(v):
                self.event_durations[k] = numpy.array(v, numpy.float64)
                self.event_durations[k][:,0] *= self.dt + self.t0
        
