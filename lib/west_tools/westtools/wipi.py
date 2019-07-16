# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import os, sys
import scipy.sparse as sp

from westtools import Plotter
import itertools

# A useful dataclass used as a wrapper for w_ipa to facilitate
# ease-of-use in ipython/jupyter notebooks/sessions.
# It basically just wraps up numpy arrays and dicts.

class WIPIDataset(object):
    def __init__(self, raw, key):
        self.__dict__ = {}
        self.raw = raw
        self.name = key
    def __repr__(self):
        if isinstance(self.__dict__['raw'], dict):
            return repr(self.__dir__())
        else:
            return repr(self.raw)
    def __getitem__(self, value):
        if not isinstance(value, str):
            return self.__dict__['raw'][value]
        if value in list(self.__dict__['raw'].keys()):
            return self.__dict__['raw'][value]
        elif value in list(self.__dict__.keys()):
            return self.__dict__[value]
    def __setitem__(self, key, value):
        self.__dict__[key] = value
    def __getattr__(self, value):
        # Check if it's an attribute of the underlying datatype.
        # If not, just use the getitem function.
        if value in dir(self.__dict__['raw']):
            return getattr(self.__dict__['raw'], value)
        else:
            return self.__getitem__(value)
    def __setattr__(self, key, value):
        self.__dict__[key] = value
    def __dir__(self):
        dict_keys = list(self.__dict__.keys())
        remove = ['raw', 'name', '__dict__', 'plotter']
        for i in remove:
            try:
                dict_keys.remove(str(i))
            except:
                pass
        # We don't enforce that this is a dictionary.
        if isinstance(self.__dict__['raw'], dict):
            return sorted(set(list(self.raw.keys()) + dict_keys))
        else:
            return sorted(set(dict_keys))
    def keys(self):
        print(self.__dir__())
    # We want to override the basic math functions, now, so... this is only valid for numpy sets.
    def __add__(self, other):
        return self.__dict__['raw'] + other
    def __radd__(self, other):
        return other + self.__dict__['raw']
    def __sub__(self, other):
        return self.__dict__['raw'] - other
    def __rsub__(self, other):
        return other - self.__dict__['raw']
    def __mul__(self, other):
        return self.__dict__['raw'] * other
    def __rmul__(self, other):
        return other * self.__dict__['raw']
    def __div__(self, other):
        return self.__dict__['raw'] / other
    def __floordiv__(self, other):
        return self.__dict__['raw'] // other
    def __rdiv__(self, other):
        return other / self.__dict__['raw']
    def __mod__(self, other):
        return self.__dict__['raw'] % other
    def __pow__(self, other):
        return self.__dict__['raw'] ** other
    def __lshift__(self, other):
        return self.__dict__['raw'] << other
    def __rshift__(self, other):
        return self.__dict__['raw'] >> other
    def __and__(self, other):
        return self.__dict__['raw'] & other
    def __eq__(self, other):
        return self.__dict__['raw'] == other
    def __ne__(self, other):
        return self.__dict__['raw'] != other
    def __lt__(self, other):
        return self.__dict__['raw'] < other
    def __gt__(self, other):
        return self.__dict__['raw'] > other
    def __le__(self, other):
        return self.__dict__['raw'] <= other
    def __ge__(self, other):
        return self.__dict__['raw'] >= other
    def __xor__(self, other):
        return self.__dict__['raw'] ^ other
    def __or__(self, other):
        return self.__dict__['raw'] | other
    #def __iadd__(self, other):
    #    return self.__dict__['raw'] += other
    #def __isub__(self, other):
    #    return self.__dict__['raw'] -= other
    #def __imul__(self, other):
    #    return self.__dict__['raw'] *= other
    #def __idiv__(self, other):
    #    return self.__dict__['raw'] /= other
    #def __ifloordiv__(self, other):
    #    return self.__dict__['raw'] //= other
    #def __imod__(self, other):
    #    return self.__dict__['raw'] %= other
    #def __ipow__(self, other):
    #    return self.__dict__['raw'] **= other
    #def __ilshift__(self, other):
    #    return self.__dict__['raw'] <<= other
    #def __irshift__(self, other):
    #    return self.__dict__['raw'] >>= other
    #def __iand__(self, other):
    #    return self.__dict__['raw'] &= other
    #def __ixor__(self, other):
    #    return self.__dict__['raw'] ^= other
    #def __ior__(self, other):
    #    return self.__dict__['raw'] |= other

# Similar to the above, but slightly expanded to contain information from analysis files.
class KineticsIteration(object):
    def __init__(self, kin_h5file, index, assign, iteration=-1):
        self.__dict__ = {}
        self.h5file = kin_h5file
        # Keys:
        _2D_h5keys = [ 'conditional_flux_evolution', 'rate_evolution' ]
        _1D_h5keys = [ 'state_pop_evolution', 'color_prob_evolution', 'target_flux_evolution' ]
        for key in _2D_h5keys:
            try:
                self.__dict__[key] = self.__2D_with_error__(key, index, assign)
            except:
                self.__dict__[key] = None
        for key in _1D_h5keys:
            try:
                self.__dict__[key] = self.__1D_with_error__(key, index, assign)
            except:
                self.__dict__[key] = None
        try:
            self.__dict__['total_fluxes'] = WIPIDataset(raw=np.array(self.h5file['total_fluxes']), key='total_fluxes')
            # We'll have to update this to make things better...
            #self.__dict__['total_fluxes'].plotter = Plotter(self.h5file['total_fluxes'][...], 'Total Fluxes', iteration=iteration, interface='text')
            #self.__dict__['total_fluxes'].plot = self.__dict__['total_fluxes'].plotter.plot
        except:
            pass

    def __repr__(self):
        return repr(self.__dir__())
    def __getitem__(self, value):
        if value in list(self.__dict__.keys()):
            return self.__dict__[value]
    def __setitem__(self, key, value):
        self.__dict__[key] = value
    def __getattr__(self, value):
        if value in list(self.__dict__.keys()):
            return self.__dict__[value]
    def __setattr__(self, key, value):
        self.__dict__[key] = value
    def __dir__(self):
        dict_keys = list(self.__dict__.keys())
        # We don't want to show the plotter class; just the plot function
        remove = [ 'h5file', '__dict__']
        for i in remove:
            try:
                dict_keys.remove(str(i))
            except:
                pass
        return sorted(set(dict_keys))
    def keys(self):
        print(self.__dir__())

    # We seriously need to rename this.
    # It's similar to the global WIPDataset, but has some nice pretty print functions.
    class __custom_dataset__(object):
        # This is just allow it to be indexed via properties.
        # Not a huge thing, but whatever.
        def __init__(self, raw, assign, key):
            self.__dict__ = {}
            self.raw = raw
            self.name = key
            self.assign = assign
            self.nstates = assign.attrs['nstates']
            self.dim = len(raw.shape)
        def __repr__(self):
            return repr(self.__dir__())
        def __getitem__(self, value):
            if value in self.__dict__['raw'].dtype.names:
                return self.__dict__['raw'][value]
            elif value in list(self.__dict__.keys()):
                return self.__dict__[value]
        def __setitem__(self, key, value):
            self.__dict__[key] = value
        def __getattr__(self, value):
            if value in self.__dict__['raw'].dtype.names:
                return self.__dict__['raw'][value]
            elif value in list(self.__dict__.keys()):
                return self.__dict__[value]
        def __setattr__(self, key, value):
            self.__dict__[key] = value
        def __dir__(self):
            dict_keys = list(self.__dict__.keys())
            # We don't want to show the plotter class; just the plot function
            remove = ['assign', 'dim', 'nstates', 'plotter', '__dict__']
            for i in remove:
                try:
                    dict_keys.remove(str(i))
                except:
                    pass
            return sorted(set(list(self.raw.dtype.names) + dict_keys))
        def keys(self):
            print(self.__dir__())
        def _repr_pretty_(self, p, cycle):
            if self.dim == 1:
                return self._1D_repr_pretty_(p, cycle)
            if self.dim == 2:
                return self._2D_repr_pretty_(p, cycle)
        def _1D_repr_pretty_(self, p, cycle):
           # We're just using this as a way to print things in a pretty way.  They can still be indexed appropriately.
           # Stolen shamelessly from westtools/kinetics_tool.py
            maxlabellen = max(list(map(len,self.assign['state_labels'])))
            p.text('')
            p.text('{name} data:\n'.format(name=self.name))
            for istate in range(self.nstates):
                p.text('{:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1\n'
                    .format(self.assign['state_labels'][istate],
                    self.raw['expected'][istate],
                    self.raw['ci_lbound'][istate],
                    self.raw['ci_ubound'][istate],
                    maxlabellen=maxlabellen))
            p.text('To access data, index via the following names:\n')
            p.text(str(self.__dir__()))
            return " "
        def _2D_repr_pretty_(self, p, cycle):
            # We're just using this as a way to print things in a pretty way.  They can still be indexed appropriately.
            # Stolen shamelessly from westtools/kinetics_tool.py
            maxlabellen = max(list(map(len,self.assign['state_labels'])))
            p.text('')
            p.text('{name} data:\n'.format(name=self.name))
            for istate in range(self.nstates):
                for jstate in range(self.nstates):
                    if istate == jstate: continue
                    p.text('{:{maxlabellen}s} -> {:{maxlabellen}s}: mean={:21.15e} CI=({:21.15e}, {:21.15e}) * tau^-1\n'
                        .format(self.assign['state_labels'][istate], self.assign['state_labels'][jstate],
                        self.raw['expected'][istate, jstate],
                        self.raw['ci_lbound'][istate, jstate],
                        self.raw['ci_ubound'][istate, jstate],
                        maxlabellen=maxlabellen))
            p.text('To access data, index via the following names:\n')
            p.text(str(self.__dir__()))
            return " "


    def __2D_with_error__(self, h5key, index, assign):
        # Check the start and stop, calculate the block size, and index appropriately.
        # While we could try and automatically generate this above, it's a little more consistent to try it here.
        # This should show the first block for which the current iteration has contributed data.
        self.step_iter = (self.h5file[h5key]['iter_stop'][0] - self.h5file[h5key]['iter_start'][0])[1,0]
        value = ((index-self.h5file.attrs['iter_start']) // self.step_iter)
        if value < 0:
            value = 0
        raw = self.h5file[h5key][value, :, :]
        error = (raw['ci_ubound'] - raw['ci_lbound']) / (2*raw['expected'])
        expected = raw['expected']
        raw = self.__custom_dataset__(raw, assign, h5key)
        raw.error = error
        raw.plotter = Plotter(self.h5file, h5key, iteration=value, interface='text')
        raw.plot = raw.plotter.plot
        return raw
    def __1D_with_error__(self, h5key, index, assign):
        self.step_iter = (self.h5file[h5key]['iter_stop'][0] - self.h5file[h5key]['iter_start'][0])[1]
        value = ((index-self.h5file.attrs['iter_start']) // self.step_iter)
        if value < 0:
            value = 0
        raw = self.h5file[h5key][value, :]
        error = (raw['ci_ubound'] - raw['ci_lbound']) / (2*raw['expected'])
        expected = raw['expected']
        raw = self.__custom_dataset__(raw, assign, h5key)
        raw.error = error
        raw.plotter = Plotter(self.h5file, h5key, iteration=value, interface='text')
        raw.plot = raw.plotter.plot
        return raw


class __get_data_for_iteration__(object):
    '''
    All interesting data from an iteration (current/past).  Whenever you change the scheme or iteration,
    this dictionary is automatically updated.  For the current iteration, it's keyed to the current seg_id.
    For the past iteration, it's keyed to the seg_id in the CURRENT iteration such that:

        w.current[X] & w.past[X]

    returns information about seg_id X in the current iteration and information on seg_ID X's PARENT in the
    preceding iteration.

    Can be indexed via a seg_id, or like a dictionary with the following keys:

        kinavg, weights, pcoord, auxdata (optional), parents, summary, seg_id, walkers, states, bins

    kinavg, states, and bins refer to the output from w_kinavg and w_assign for this iteration
    and analysis scheme.  They are NOT dynamics bins, but the bins defined in west.cfg.  
    
    Has the following properties:

        .minweight, .maxweight

    which return all properties of the segment that matches those criteria in the selected iteration.

    If you change the analysis scheme, so, too, will the important values.
    '''

    def __init__(self, parent, value, seg_ids = None):
        '''
        Initializes and sets the correct data.
        '''
        # We've classed this so that we can override some of the normal functions and allow indexing via seg_id
        self.__dict__ = {}
        # Is this function thread safe?
        iter_group = parent.data_reader.get_iter_group(value)
        #iter_group = parent.west['iterations/iter_{num:08d}'.format(num=value)]
        self.parent = parent
        current = {}
        current['iteration'] = value
        try: 
            if seg_ids == None:
                seg_ids = range(0, iter_group['seg_index']['weight'].shape[0])
        except:
            pass
        # Just make these easier to access.
        current['weights'] = iter_group['seg_index']['weight'][seg_ids]
        current['pcoord'] = iter_group['pcoord'][...][seg_ids, :, :]
        try:
            current['auxdata'] = {}
            for key in list(iter_group['auxdata'].keys()):
                current['auxdata'][key] = iter_group['auxdata'][key][...][seg_ids, :]
        except:
            pass
        current['parents'] = iter_group['seg_index']['parent_id'][seg_ids]
        current['summary'] = parent.data_reader.data_manager.get_iter_summary(int(value))
        current['seg_id'] = np.array(list(range(0, iter_group['seg_index'].shape[0])))[seg_ids]
        current['walkers'] = current['summary']['n_particles']
        current['states'] = parent.assign['trajlabels'][value-1, :current['walkers'], :][seg_ids]
        current['bins'] = parent.assign['assignments'][value-1, :current['walkers'], :][seg_ids]
        # Calculates the bin population for this iteration.
        nbins = parent.assign['state_map'].shape[0]
        # We have to take the 'unknown' state into account
        nstates = parent.assign['state_labels'].shape[0] + 1
        # Temporarily disabled while I sort out the fact that we shouldn't be using data from w_assign for state populations.
        #current['plot'] = Plotter(parent.direct, parent.reweight, parent.iteration, parent.assign['bin_labels'], parent.assign['state_labels'], current['populations'].states, current['populations'].bins, parent.interface)
        # Now we'll load up the results of the kinetics analysis.
        current['direct'] = KineticsIteration(parent.direct, value, parent.assign, value)
        evolution_datasets = [ 'rate_evolution', 'conditional_flux_evolution', 'state_pop_evolution', 'color_prob_evolution' , 'total_fluxes', 'target_flux_evolution']
        # We want to load these up as... oh, who knows, I suppose?
        try:
            current['reweight'] = KineticsIteration(parent.reweight, value, parent.assign, value)
            # We'll make this not a sparse matrix...
            matrix = parent.reweight['iterations/iter_{:08d}'.format(value)]
            # Assume color.
            current['instant_matrix'] = sp.coo_matrix((matrix['flux'][...], (matrix['rows'][...], matrix['cols'][...])), shape=((nbins-1)*2, (nbins-1)*2)).todense()
            reweighting = True
        except:
          # This analysis hasn't been enabled, so we'll simply return the default error message.
            current['reweight'] = parent.reweight['rate_evolution']
            current['instant_matrix'] = parent.reweight['bin_populations']
            current['matrix'] = parent.reweight['bin_populations']
            reweighting = False
        # Check if the analysis has been enabled.  If yes, make them specify dataset dictionaries.  If not, return the thing.
        if reweighting:
            for key in evolution_datasets:
                current[key] = WIPIDataset(raw={ 'direct': current['direct'][key], 'reweight': current['reweight'][key] }, key='a')
        else:
            for key in evolution_datasets:
                current[key] = WIPIDataset(raw={ 'direct': current['direct'][key] }, key='direct')

        self.raw = current
    def __repr__(self):
        '''
        Returns the dictionary containing the iteration's values.
        '''
        return repr(self.__dir__())

    def keys(self):
        '''
        Returns the keys function of the internal dictionary.
        '''
        return list(self.__dict__['raw'].keys())

    def __setitem__(self, key, value):
        self.__dict__[key] = value
    def __getattr__(self, value):
        if value in list(self.__dict__['raw'].keys()):
            return self.__dict__['raw'][value]
        elif value in list(self.__dict__.keys()):
            return self.__dict__[value]
    def __setattr__(self, key, value):
        self.__dict__[key] = value
    def __dir__(self):
        dict_keys = list(self.__dict__.keys())
        dict_keys += ['maxweight', 'minweight', 'walkers', 'aggregate_walkers', 'successful_trajectories']
        remove = ['__dict__']
        for i in remove:
            try:
                dict_keys.remove(str(i))
            except:
                pass
        return sorted(set(list(self.__dict__['raw'].keys()) + dict_keys))

    @property
    def maxweight(self):
        '''
        Returns information about the segment which has the largest weight for this iteration.
        '''
        # Is there a faster or cleaner way to do this?  Ah, maybe.
        walker = np.where(self.raw['weights'] == np.max(self.raw['weights']))[0][0]
        return self.__getitem__(walker)

    @property
    def minweight(self):
        '''
        Returns information about the segment which has the smallest weight for this iteration.
        '''
        walker = np.where(self.raw['weights'] == np.min(self.raw['weights']))[0][0]
        return self.__getitem__(walker)

    @property
    def successful_trajectories(self):
        '''
        Returns which trajectories are successful.
        '''
        #walker = np.where(self.raw['weights'] == np.min(self.raw['weights']))[0][0]
        # Find where we have a transition....
        state_changes = np.where(self.raw['states'][:,:-1] != self.raw['states'][:,1:])
        walkers = state_changes[0]
        # The index of the state change.
        new_states = state_changes[1] + 1
        old_states = state_changes[1]
        walker = {}
        for z, (i, j) in enumerate(zip(old_states, new_states)):
            #if self.raw['states'][walkers[z], i] == istate and self.raw['states'][walkers[z], j] == jstate:
            istate = self.raw['states'][walkers[z], i] 
            jstate = self.raw['states'][walkers[z], j]
            #print(z,i,j, istate, jstate)
            try:
                walker[istate,jstate].append(walkers[z])
            except:
                walker[istate,jstate] = [walkers[z]]

        walker = WIPIDataset(raw=walker, key=None)
        return walker

    @property
    def walkers(self):
        '''
        The number of walkers active in the current iteration.
        '''
        # Returns number of walkers for iteration X.  Assumes current iteration, but can go with different one.
        # Make this just... yeah, put this elsewhere.
        return self.parent.west['summary']['n_particles'][self.iteration-1]

    @property
    def aggregate_walkers(self):
        return self.parent.west['summary']['n_particles'][:self.iteration].sum()


    def __getitem__(self, value):
        '''
        Responsible for handling whether this is treated like a dictionary of data sets, or an array of walker data.
        '''
        # Check to see if we're indexing via any of the active string types.  We should probably break it down via string or int, instead of 'what exists and what doesn't', but it works for now.
        active_items = ['kinavg', 'statepops', 'weights', 'pcoord', 'auxdata', 'parents', 'summary', 'seg_id', 'walkers', 'states', 'bins', 'populations', 'plot', 'instant_matrix', 'kinrw', 'matrix', 'rwstatepops']
        #if value in active_items:
        if isinstance(value, str):
            # This should handle everything.  Otherwise...
            try:
                return self.raw[value]
            except:
                print('{} is not a valid data structure.'.format(value))
        elif isinstance(value, int) or isinstance(value, np.int64):
            # Otherwise, we assume they're trying to index for a seg_id.
            if value < self.walkers:
                current = {}
                current['plotter'] = {}
                for i in ['pcoord']:
                    current[i] = WIPIDataset(raw=self.raw[i][value,:,:], key=i)
                    current[i].plotter = Plotter(self.raw[i][value,:,:], i, iteration=self.iteration, interface='text')
                    current[i].plot = current[i].plotter.plot

                current['states'] = self.raw['states'][value, :]
                current['bins'] = self.raw['bins'][value, :]
                current['parents'] = self.raw['parents'][value]
                current['seg_id'] = self.raw['seg_id'][value]
                current['weights'] = self.raw['weights'][value]
                try:
                    current['auxdata'] = {}
                    for key in list(self.raw['auxdata'].keys()):
                        current['auxdata'][key] = self.raw['auxdata'][key][value]
                except:
                    pass
                current = WIPIDataset(current, 'Segment {} in Iter {}'.format(value, self.iteration))
                return current
            else:
                print('INVALID SEG_ID {}.  SEG_ID should be less than {}.'.format(value, self.walkers))

# This handles the 'schemes', and all assorted data.
class WIPIScheme(object):
    def __init__(self, scheme, name, parent, settings):
        self.__dict__ = {}
        self.raw = scheme
        #self.name = parent._schemename
        self.__analysis_schemes__ = scheme
        self.iteration = parent.iteration
        self.__dict__['name'] = None
        self.__settings = settings
        # Are these necessary?  We'll try to edit these out.
        self.parent = parent
        self.data_reader = parent.data_reader
    def __setattr__(self, key, value):
        self.__dict__[key] = value
    def __repr__(self):
        return self.__str__()
    def __str__(self):
        # Right now, this returns w.scheme, NOT necessarily what we're pulling from...
        # So you can rely on this, but it's confusing.
        if self.name!= None:
            # Set it to None, then return the original value.
            rtn_string = self.name
            self.name = None
            return rtn_string
        else:
            return str(self.scheme)
    def __getitem__(self, value):
        if not isinstance(value, str):
            for ischeme, schemename in enumerate(self.__dict__['raw'].keys()):
                if ischeme == value:
                    value = schemename
        # Check for some weird Ipython stuff.
        if '_ipython' in value:
            return self
        self.name = None
        if value in list(self.__dict__['raw'].keys()):
            # If we have it in there...
            self.name = value
            return self
        elif value in list(self.__dict__.keys()):
            self.name = value
            return self
        elif value in self.__dir__():
            self.name = value
            return self

    def __getattr__(self, value):
        return self.__getitem__(value)

    def __dir__(self):
        dict_keys = ['assign', 'direct', 'state_labels', 'bin_labels', 'west', 'reweight', 'current', 'past', 'iteration']
        if self.name != None:
            return sorted(set(dict_keys))
        else:
            return sorted(set(self.__analysis_schemes__.keys()))

    @property
    def scheme(self):
        self.name = None
        return self.parent._schemename

    @property
    def list_schemes(self):
        '''
        Lists what schemes are configured in west.cfg file.
        Schemes should be structured as follows, in west.cfg:

        west:
          system:
            analysis:
              directory: analysis
              analysis_schemes:
                scheme.1:
                  enabled: True
                  states:
                    - label: unbound
                      coords: [[7.0]]
                    - label: bound
                      coords: [[2.7]]
                  bins:
                    - type: RectilinearBinMapper
                      boundaries: [[0.0, 2.80, 7, 10000]]
        '''
        print("The following schemes are available:")
        print("")
        for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
            print('{}. Scheme: {}'.format(ischeme, scheme))
        print("")
        print("Set via name, or via the index listed.")
        print("")
        print("Current scheme: {}".format(self.scheme))

    @property
    def iteration(self):
        return self.parent.iteration

    @property
    def assign(self):
        return self.__analysis_schemes__[str(self.name)]['assign']

    @property
    def direct(self):
        '''
        The output from w_direct.py from the current scheme.
        '''
        return self.__analysis_schemes__[str(self.name)]['direct']

    @property
    def state_labels(self):
        print("State labels and definitions!")
        for istate, state in enumerate(self.assign['state_labels']):
            print('{}: {}'.format(istate, state))
        print('{}: {}'.format(istate+1, 'Unknown'))

    @property
    def bin_labels(self):
        print("Bin definitions! ")
        for istate, state in enumerate(self.assign['bin_labels']):
            print('{}: {}'.format(istate, state))

    @property
    def west(self):
        return self.data_reader.data_manager.we_h5file

    @property
    def reweight(self):
        # Need to fix this...
        if self.__settings['analysis_schemes'][str(self.name)]['postanalysis'] == True:
            return self.__analysis_schemes__[str(self.name)]['reweight']
        else:
            value = "This sort of analysis has not been enabled."
            current = { 'bin_prob_evolution': value, 'color_prob_evolution': value, 'conditional_flux_evolution': value, 'rate_evolution': value, 'state_labels': value, 'state_prob_evolution': value }
            current.update({ 'bin_populations': value, 'iterations': value })
            return current

    @property
    def current(self):
        '''
        The current iteration.  See help for __get_data_for_iteration__
        '''
        return __get_data_for_iteration__(value=self.iteration, parent=self)

    @property
    def past(self):
        '''
        The previous iteration.  See help for __get_data_for_iteration__
        '''
        if self.iteration > 1:
            return __get_data_for_iteration__(value=self.iteration - 1, seg_ids=self.current['parents'], parent=self)
        else:
            print("The current iteration is 1; there is no past.")
