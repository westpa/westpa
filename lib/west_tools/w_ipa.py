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

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=ImportWarning)
import numpy as np
import h5py

# Must be run with the WEST wrapper.
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object
import westpa
import os, sys
import w_assign, w_direct, w_reweight
warnings.filterwarnings('ignore')
import scipy.sparse as sp
import hashlib
import json
#sys.tracebacklimit = 5

from westtools import (WESTSubcommand, WESTParallelTool, WESTDataReader, WESTDSSynthesizer, BinMappingComponent, 
                       ProgressIndicatorComponent, IterRangeSelection, Plotter)

from westtools import WIPIDataset, __get_data_for_iteration__, WIPIScheme

class WIPI(WESTParallelTool):
    '''
        Welcome to w_ipa (WESTPA Interactive Python Analysis)!
        From here, you can run traces, look at weights, progress coordinates, etc.
        This is considered a 'stateful' tool; that is, the data you are pulling is always pulled
        from the current analysis scheme and iteration.
        By default, the first analysis scheme in west.cfg is used, and you are set at iteration 1.

        ALL PROPERTIES ARE ACCESSED VIA w or west
        To see the current iteration, try:

            w.iteration
            OR
            west.iteration

        to set it, simply plug in a new value.

            w.iteration = 100

        To change/list the current analysis schemes:

            w.list_schemes
            w.scheme = OUTPUT FROM w.list_schemes

        To see the states and bins defined in the current analysis scheme:

            w.states
            w.bin_labels

        All information about the current iteration is available in an object called 'current':

            w.current
            walkers, summary, states, seg_id, weights, parents, kinavg, pcoord, bins, populations, and auxdata, if it exists.

        In addition, the function w.trace(seg_id) will run a trace over a seg_id in the current iteration and return a dictionary
        containing all pertinent information about that seg_id's history.  It's best to store this, as the trace can be expensive.

        Run help on any function or property for more information!

        Happy analyzing!
                
    '''

    def __init__(self):
        super(WIPI,self).__init__()
        self.data_reader = WESTDataReader()
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager
        self.progress = ProgressIndicatorComponent()

        self._iter = 1
        self.config_required = True
        self.version = "1.0B"
        # Set to matplotlib if you want that.  But why would you?
        # Well, whatever, we'll just set it to that for now.
        self.interface = 'matplotlib'
        self._scheme = None
        global iteration

    def add_args(self, parser):
        self.progress.add_args(parser)
        self.data_reader.add_args(parser)
        rgroup = parser.add_argument_group('runtime options')
        rgroup.add_argument('--analysis-only', '-ao', dest='analysis_mode', action='store_true',
                             help='''Use this flag to run the analysis and return to the terminal.''')
        rgroup.add_argument('--reanalyze', '-ra', dest='reanalyze', action='store_true',
                             help='''Use this flag to delete the existing files and reanalyze.''')
        rgroup.add_argument('--ignore-hash', '-ih', dest='ignore_hash', action='store_true',
                             help='''Ignore hash and don't regenerate files.''')
        rgroup.add_argument('--debug', '-d', dest='debug_mode', action='store_true',
                             help='''Debug output largely intended for development.''')
        rgroup.add_argument('--terminal', '-t', dest='plotting', action='store_true',
                             help='''Plot output in terminal.''')
        # There is almost certainly a better way to handle this, but we'll sort that later.
        import argparse
        rgroup.add_argument('--f', '-f', dest='extra', default='blah',
                             help=argparse.SUPPRESS)
        
        parser.set_defaults(compression=True)

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.niters = self.data_reader.current_iteration - 1
        self.__config = westpa.rc.config
        self.__settings = self.__config['west']['analysis']
        for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
            if (self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None):
                self.scheme = scheme
        self.data_args = args
        self.analysis_mode = args.analysis_mode
        self.reanalyze = args.reanalyze
        self.ignore_hash = args.ignore_hash
        self.debug_mode = args.debug_mode
        if args.plotting:
            self.interface = 'text'

    def hash_args(self, args, extra=None, path=None):
        '''Create unique hash stamp to determine if arguments/file is different from before.'''
        '''Combine with iteration to know whether or not file needs updating.'''
        # Why are we not loading this functionality into the individual tools?
        # While it may certainly be useful to store arguments (and we may well do that),
        # it's rather complex and nasty to deal with pickling and hashing arguments through
        # the various namespaces.
        # In addition, it's unlikely that the functionality is desired at the individual tool level,
        # since we'll always just rewrite a file when we call the function.
        #return hashlib.md5(pickle.dumps([args, extra])).hexdigest()
        # We don't care about the path, so we'll remove it.
        # Probably a better way to do this, but who cares.
        cargs = list(args)
        for iarg, arg in enumerate(cargs):
            if path in arg:
                cargs[iarg] = arg.replace(path,'').replace('/', '')
            if arg == '--disable-averages':
                cargs.remove('--disable-averages')
        to_hash = cargs + [extra]
        #print(args)
        #print(to_hash)
        #print(str(to_hash).encode('base64'))
        if self.debug_mode:
            for iarg, arg in enumerate(to_hash):
                if not isinstance(arg, list):
                    print('arg {num:02d} -- {arg:<20}'.format(num=iarg, arg=arg))
                else:
                    for il, l in enumerate(arg):
                        print('arg {num:02d} -- {arg:<20}'.format(num=il+iarg, arg=l))
            #print('args: {}'.format(to_hash))
        # This SHOULD produce the same output, maybe?  That would be nice, anyway.
        # But we'll need to test it more.
        return hashlib.md5(str(to_hash).encode('base64')).hexdigest()

    def stamp_hash(self, h5file_name, new_hash):
        '''Loads a file, stamps it, and returns the opened file in read only'''
        h5file = h5io.WESTPAH5File(h5file_name, 'r+')
        h5file.attrs['arg_hash'] = new_hash
        h5file.close()
        h5file = h5io.WESTPAH5File(h5file_name, 'r')
        return h5file

    def analysis_structure(self):
        '''
        Run automatically on startup.  Parses through the configuration file, and loads up all the data files from the different 
        analysis schematics.  If they don't exist, it creates them automatically by hooking in to existing analysis routines 
        and going from there.  

        It does this by calling in the make_parser_and_process function for w_{assign,reweight,direct} using a custom built list
        of args.  The user can specify everything in the configuration file that would have been specified on the command line.

        For instance, were one to call w_direct as follows:

            w_direct --evolution cumulative --step-iter 1 --disable-correl

        the west.cfg would look as follows:

        west:
          analysis:
            w_direct:
              evolution: cumulative
              step_iter: 1
              extra: ['disable-correl']

        Alternatively, if one wishes to use the same options for both w_direct and w_reweight, the key 'w_direct' can be replaced
        with 'kinetics'.
        '''
        # Make sure everything exists.
        try:
            os.mkdir(self.__settings['directory'])
        except:
            pass
        # Now, check to see whether they exist, and then load them.
        self.__analysis_schemes__ = {}
        # We really need to implement some sort of default behavior if an analysis scheme isn't set.
        # Right now, we just crash.  That isn't really graceful.
        for scheme in self.__settings['analysis_schemes']:
            if self.__settings['analysis_schemes'][scheme]['enabled']:
                if self.work_manager.running == False:
                    self.work_manager.startup()
                path = os.path.join(os.getcwd(), self.__settings['directory'], scheme)
                #if 'postanalysis' in self.__settings['analysis_schemes'][scheme] and 'postanalysis' in self.__settings['postanalysis']:
                # Should clean this up.  But it uses the default global setting if a by-scheme one isn't set.
                if 'postanalysis' in self.__settings:
                    if 'postanalysis' in self.__settings['analysis_schemes'][scheme]:
                        pass
                    else:
                        self.__settings['analysis_schemes'][scheme]['postanalysis'] = self.__settings['postanalysis']
                try:
                    os.mkdir(path)
                except:
                    pass
                self.__analysis_schemes__[scheme] = {}
                try:
                    if self.__settings['analysis_schemes'][scheme]['postanalysis'] == True or self.__settings['postanalysis'] == True:
                        analysis_files = ['assign', 'direct', 'reweight']
                    else:
                        analysis_files = ['assign', 'direct']
                except:
                    analysis_files = ['assign', 'direct']
                    self.__settings['analysis_schemes'][scheme]['postanalysis'] = False
                reanalyze_kinetics = False
                assign_hash = None
                for name in analysis_files:
                    arg_hash = None
                    if self.reanalyze == True:
                        reanalyze_kinetics = True
                        try:
                            os.remove(os.path.join(path, '{}.h5'.format(name)))
                        except:
                            pass
                    else:
                        try:
                            # Try to load the hash.  If we fail to load the hash or the file, we need to reload.
                            #if self.reanalyze == True:
                            #    raise ValueError('Reanalyze set to true.')
                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            arg_hash = self.__analysis_schemes__[scheme][name].attrs['arg_hash']
                            if name == 'assign':
                                assign_hash = arg_hash
                        except:
                            pass
                            # We shouldn't rely on this.
                            # self.reanalyze = True
                    if True:
                        if name == 'assign':
                            assign = w_assign.WAssign()

                            w_assign_config = { 'output': os.path.join(path, '{}.h5'.format(name))}
                            try:
                                w_assign_config.update(self.__settings['w_assign'])
                            except:
                                pass
                            try:
                                w_assign_config.update(self.__settings['analysis_schemes'][scheme]['w_assign'])
                            except:
                                pass
                            args = []
                            for key,value in w_assign_config.items():
                                if key != 'extra':
                                    args.append(str('--') + str(key).replace('_', '-'))
                                    args.append(str(value))
                            # This is for stuff like disabling correlation analysis, etc.
                            if 'extra' in list(w_assign_config.keys()):
                                # We're sorting to ensure that the order doesn't matter.
                                for value in sorted(w_assign_config['extra']):
                                    args.append(str('--') + str(value).replace('_', '-'))
                            # We're just calling the built in function.
                            # This is a lot cleaner than what we had in before, and far more workable.
                            args.append('--config-from-file')
                            args.append('--scheme-name')
                            args.append('{}'.format(scheme))
                            # Why are we calling this if we're not sure we're remaking the file?
                            # We need to load up the bin mapper and states and see if they're the same.
                            assign.make_parser_and_process(args=args)
                            import pickle
                            #new_hash = self.hash_args(args=args, path=path, extra=[self.niters, pickle.dumps(assign.binning.mapper), assign.states])
                            # We need to encode it properly to ensure that some OS specific thing doesn't kill us.  Same goes for the args, ultimately.
                            # Mostly, we just need to ensure that we're consistent.
                            new_hash = self.hash_args(args=args, path=path, extra=[int(self.niters), pickle.dumps(assign.binning.mapper).encode('base64'), str(assign.states).encode('base64')])
                            # Let's check the hash.  If the hash is the same, we don't need to reload.
                            if self.debug_mode == True:
                                print('{:<10}: old hash, new hash -- {}, {}'.format(name, arg_hash, new_hash))
                            if self.ignore_hash == False and (arg_hash != new_hash or self.reanalyze == True):
                                # If the hashes are different, or we need to reanalyze, delete the file.
                                try:
                                    os.remove(os.path.join(path, '{}.h5'.format(name)))
                                except:
                                    pass
                                print('Reanalyzing file {}.h5 for scheme {}.'.format(name, scheme))
                                #reanalyze_kinetics = True
                                # We want to use the work manager we have here.  Otherwise, just let the tool sort out what it needs, honestly.
                                assign.work_manager = self.work_manager

                                assign.go()
                                assign.data_reader.close()

                                # Stamp w/ hash, then reload as read only.
                                self.__analysis_schemes__[scheme][name] = self.stamp_hash(os.path.join(path, '{}.h5'.format(name)), new_hash)
                            del(assign)
                            # Update the assignment hash.
                            assign_hash = new_hash

                        # Since these are all contained within one tool, now, we want it to just... load everything.
                        if name == 'direct' or name == 'reweight':
                            assignment_file = self.__analysis_schemes__[scheme]['assign']
                            if name == 'direct':
                                analysis = w_direct.WDirect()
                            if name == 'reweight':
                                analysis = w_reweight.WReweight()
                            
                            analysis_config = { 'assignments': os.path.join(path, '{}.h5'.format('assign')), 'output': os.path.join(path, '{}.h5'.format(name)), 'kinetics': os.path.join(path, '{}.h5'.format(name))}

                            # Pull from general analysis options, then general SPECIFIC options for each analysis,
                            # then general options for that analysis scheme, then specific options for the analysis type in the scheme.

                            try:
                                analysis_config.update(self.__settings['kinetics'])
                            except:
                                pass
                            try:
                                analysis_config.update(self.__settings['w_{}'.format(name)])
                            except:
                                pass
                            try:
                                analysis_config.update(self.__settings['analysis_schemes'][scheme]['kinetics'])
                            except:
                                pass
                            try:
                                analysis_config.update(self.__settings['analysis_schemes'][scheme]['w_{}'.format(name)])
                            except:
                                pass

                            # We're pulling in a default set of arguments, then updating them with arguments from the west.cfg file, if appropriate, after setting the appropriate command
                            # Then, we call the magic function 'make_parser_and_process' with the arguments we've pulled in.
                            # The tool has no real idea it's being called outside of its actual function, and we're good to go.
                            args = ['all']
                            for key,value in analysis_config.items():
                                if key != 'extra':
                                    args.append(str('--') + str(key).replace('_', '-'))
                                    args.append(str(value))
                            # This is for stuff like disabling correlation analysis, etc.
                            if 'extra' in list(analysis_config.keys()):
                                for value in sorted(analysis_config['extra']):
                                    args.append(str('--') + str(value).replace('_', '-'))
                            # We want to not display the averages, so...
                            args.append('--disable-averages')
                            new_hash = self.hash_args(args=args, path=path, extra=[int(self.niters), assign_hash])
                            #if arg_hash != new_hash or self.reanalyze == True or reanalyze_kinetics == True:
                            if self.debug_mode == True:
                                print('{:<10}: old hash, new hash -- {}, {}'.format(name, arg_hash, new_hash))
                            if self.ignore_hash == False and (arg_hash != new_hash or reanalyze_kinetics == True):
                                try:
                                    os.remove(os.path.join(path, '{}.h5'.format(name)))
                                except:
                                    pass
                                print('Reanalyzing file {}.h5 for scheme {}.'.format(name, scheme))
                                analysis.make_parser_and_process(args=args)
                                # We want to hook into the existing work manager.
                                analysis.work_manager = self.work_manager

                                analysis.go()

                                # Open!
                                self.__analysis_schemes__[scheme][name] = self.stamp_hash(os.path.join(path, '{}.h5'.format(name)), new_hash)
                            del(analysis)

        # Make sure this doesn't get too far out, here.  We need to keep it alive as long as we're actually analyzing things.
        # self.work_manager.shutdown()
        print("")
        print("Complete!")

    @property
    def assign(self):
        return self.__analysis_schemes__[str(self.scheme)]['assign']

    @property
    def direct(self):
        """
        The output from w_kinavg.py from the current scheme.
        """
        return self.__analysis_schemes__[str(self.scheme)]['direct']


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
        if self.__settings['analysis_schemes'][str(self.scheme)]['postanalysis'] == True:
            return self.__analysis_schemes__[str(self.scheme)]['reweight']
        else:
            value = "This sort of analysis has not been enabled."
            current = { 'bin_prob_evolution': value, 'color_prob_evolution': value, 'conditional_flux_evolution': value, 'rate_evolution': value, 'state_labels': value, 'state_prob_evolution': value }
            current.update({ 'bin_populations': value, 'iterations': value })
            return current

    @property
    def scheme(self):
        '''
        Returns and sets what scheme is currently in use.
        To see what schemes are available, run:

            w.list_schemes

        '''
        # Let's do this a few different ways.
        # We want to return things about the DIFFERENT schemes, if possible.
        if self._scheme == None:
            self._scheme = WIPIScheme(scheme=self.__analysis_schemes__, name=self._schemename, parent=self, settings=self.__settings)

        # This just ensures that when we call it, it's clean.
        self._scheme.name = None
        return self._scheme

    @scheme.setter
    def scheme(self, scheme):
        self._future = None
        self._current = None
        self._past = None
        if scheme in self.__settings['analysis_schemes']:
            pass
        else:
            for ischeme, schemename in enumerate(self.__settings['analysis_schemes']):
                if ischeme == scheme:
                    scheme = schemename
        if self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None:
            self._schemename = scheme
        else:
            print("Scheme cannot be changed to scheme: {}; it is not enabled!".format(scheme))

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
        #print("The following schemes are available:")
        #print("")
        #for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
        #    print('{}. Scheme: {}'.format(ischeme, scheme))
        #print("")
        #print("Set via name, or via the index listed.")
        #print("")
        #print("Current scheme: {}".format(self.scheme))
        self._scheme.list_schemes

    @property
    def iteration(self):
        '''
        Returns/sets the current iteration.
        '''
        #print("The current iteration is {}".format(self._iter))
        return self._iter

    @iteration.setter
    def iteration(self, value):
        print("Setting iteration to iter {}.".format(value))
        if value <= 0:
            print("Iteration must begin at 1.")
            value = 1
        if value > self.niters:
            print("Cannot go beyond {} iterations!".format(self.niters))
            print("Setting to {}".format(self.niters))
            value = self.niters
        # We want to trigger a rebuild on our current/past/future bits.
        # The scheme should automatically reset to the proper iteration, but
        # future needs to be manually triggered.
        self._iter = value
        self._future = None
        return self._iter


    @property
    def current(self):
        '''
        The current iteration.  See help for __get_data_for_iteration__
        '''
        return self.scheme[self.scheme.scheme].current

    @property
    def past(self):
        '''
        The previous iteration.  See help for __get_data_for_iteration__
        '''
        return self.scheme[self.scheme.scheme].past


    def trace(self, seg_id):
        '''
        Runs a trace on a seg_id within the current iteration, all the way back to the beginning,
        returning a dictionary containing all interesting information:

            seg_id, pcoord, states, bins, weights, iteration, auxdata (optional)

        sorted in chronological order.


        Call with a seg_id.
        '''
        if seg_id >= self.current.walkers:
            print("Walker seg_id # {} is beyond the max count of {} walkers.".format(seg_id, self.current.walkers))
            return 1
        pi = self.progress.indicator
        with pi:
            pi.new_operation('Tracing scheme:iter:seg_id {}:{}:{}'.format(self.scheme, self.iteration, seg_id), self.iteration)
            current = { 'seg_id': [], 'pcoord': [], 'states': [], 'weights': [], 'iteration': [], 'bins': [] }
            keys = []
            try:
                current['auxdata'] = {}
                for key in list(self.current['auxdata'].keys()):
                    current['auxdata'][key] = []
                    key = []
            except:
                pass
            for iter in reversed(list(range(1, self.iteration+1))):
                iter_group = self.data_reader.get_iter_group(iter)
                particles = self.data_reader.data_manager.get_iter_summary(int(iter))['n_particles']
                current['pcoord'].append(iter_group['pcoord'][seg_id, :, :])
                current['states'].append(self.assign['trajlabels'][iter-1, seg_id,:])
                current['bins'].append(self.assign['assignments'][iter-1, seg_id,:])
                current['seg_id'].append(seg_id)
                current['weights'].append(iter_group['seg_index']['weight'][seg_id])
                current['iteration'].append(iter)
                try:
                    for key in keys:
                        current['auxdata'][key].append(iter_group['auxdata'][key][seg_id])
                except:
                    pass
                seg_id = iter_group['seg_index']['parent_id'][seg_id]
                if seg_id < 0:
                    # Necessary for steady state simulations.  This means they started in that iteration.
                    break
                pi.progress += 1
        current['seg_id'] = list(reversed(current['seg_id']))
        current['iteration'] = list(reversed(current['iteration']))
        current['states'] = np.concatenate(np.array(list(reversed(current['states']))))
        current['bins'] = np.concatenate(np.array(list(reversed(current['bins']))))
        current['weights'] = np.array(list(reversed(current['weights'])))
        current['pcoord'] = np.concatenate(np.array(list(reversed(current['pcoord']))))
        try:
            for key in keys():
                current['auxdata'][key] = np.concatenate(np.array(list(reversed(current['auxdata'][key]))))
        except:
            pass
        current['state_labels'] = self.assign['state_labels']
        for i in ['pcoord', 'states', 'bins', 'weights']:
            current[i] = WIPIDataset(raw=current[i], key=i)
            if i == 'weights':
                current[i].plotter = Plotter(np.log10(current[i].raw), str('log10 of ' + str(i)), iteration=current[i].raw.shape[0], interface=self.interface)
            else:
                current[i].plotter = Plotter(current[i].raw, i, iteration=current[i].raw.shape[0], interface=self.interface)
            current[i].plot = current[i].plotter.plot
        return WIPIDataset(raw=current, key=seg_id)

    @property
    def future(self, value=None):
        '''
        Similar to current/past, but keyed differently and returns different datasets.
        See help for Future.
        '''
        if self._future == None:
            self._future = self.Future(raw=self.__get_children__(), key=None)
            self._future.iteration = self.iteration+1
        return self._future

    class Future(WIPIDataset):

        # This isn't a real fancy one.
        def __getitem__(self, value):
            if isinstance(value, str):
                print(list(self.__dict__.keys()))
                try:
                    return self.__dict__['raw'][value]
                except:
                    print('{} is not a valid data structure.'.format(value))
            elif isinstance(value, int) or isinstance(value, np.int64):
                # Otherwise, we assume they're trying to index for a seg_id.
                #if value < self.parent.walkers:
                current = {}
                seg_items = ['weights', 'pcoord', 'auxdata', 'parents', 'seg_id', 'states']
                current['pcoord'] = self.__dict__['raw']['pcoord'][value]
                current['states'] = self.__dict__['raw']['states'][value]
                current['bins'] = self.__dict__['raw']['bins'][value]
                current['parents'] = self.__dict__['raw']['parents'][value]
                current['seg_id'] = self.__dict__['raw']['seg_id'][value]
                current['weights'] = self.__dict__['raw']['weights'][value]
                try:
                    current['auxdata'] = {}
                    for key in list(self.__dict__['raw']['auxdata'].keys()):
                        current['auxdata'][key] = self.__dict__['raw']['auxdata'][key][value]
                except:
                    pass
                current = WIPIDataset(current, 'Segment {} in Iter {}'.format(value, self.iteration))
                return current

    def __get_children__(self):
        '''
        Returns all information about the children of a given walker in the current iteration.
        Used to generate and create the future object, if necessary.
        '''
        
        if self.iteration == self.niters:
            print("Currently at iteration {}, which is the max.  There are no children!".format(self.iteration))
            return 0
        iter_data = __get_data_for_iteration__(value=self.iteration+1, parent=self)
        future = { 'weights': [], 'pcoord': [], 'parents': [], 'summary': iter_data['summary'], 'seg_id': [], 'walkers': iter_data['walkers'], 'states': [], 'bins': [] }
        for seg_id in range(0, self.current.walkers):
            children = np.where(iter_data['parents'] == seg_id)[0]
            if len(children) == 0:
                error = "No children for seg_id {}.".format(seg_id)
                future['weights'].append(error)
                future['pcoord'].append(error)
                future['parents'].append(error)
                future['seg_id'].append(error)
                future['states'].append(error)
                future['bins'].append(error)
            else:
                # Now, we're gonna put them in the thing.
                value = self.iteration+1 
                future['weights'].append(iter_data['weights'][children])
                future['pcoord'].append(iter_data['pcoord'][...][children, :, :])
                try:
                    aux_data = iter_data['auxdata'][...][children, :, :]
                    try:
                        future['aux_data'].append(aux_data)
                    except:
                        future['aux_data'] = aux_data
                except:
                    pass
                future['parents'].append(iter_data['parents'][children])
                future['seg_id'].append(iter_data['seg_id'][children])
                future['states'].append(self.assign['trajlabels'][value-1, children, :])
                future['bins'].append(self.assign['assignments'][value-1, children, :])
        return future

    def go(self):
        '''
        Function automatically called by main() when launched via the command line interface.
        Generally, call main, not this function.
        '''
        print("")
        print("Welcome to w_ipa (WESTPA Interactive Python Analysis) v. {}!".format(w.version))
        print("Run w.introduction for a more thorough introduction, or w.help to see a list of options.")
        print("Running analysis & loading files.")
        self.data_reader.open()
        self.analysis_structure()
        # Seems to be consistent with other tools, such as w_assign.  For setting the iterations.
        self.data_reader.open()
        self.niters = self.data_reader.current_iteration - 1
        self.iteration = self.niters
        try:
            print('Your current scheme, system and iteration are : {}, {}, {}'.format(w.scheme, os.getcwd(), w.iteration))
        except:
            pass

    @property
    def introduction(self):
        '''
        Just spits out an introduction, in case someone doesn't call help.
        '''
        help_string = '''
        Call as a dictionary item or a .attribute:

        w.past, w.current, w.future:
            
            {current}

        Raw schemes can be accessed as follows:

            w.scheme.{scheme_keys}

            and contain mostly the same datasets associated with w.

        The following give raw access to the h5 files associated with the current scheme

        w.west
        w.assign
        w.direct
        w.reweight

        OTHER:

        {w}

        '''.format(current=self.__format_keys__(self.current.__dir__(), split=' ', offset=12), scheme_keys=self.__format_keys__(list(self._scheme.raw.keys())),
                   w=self.__format_keys__(self.__dir__(), offset=8, max_length=0, split='', prepend='w.'))
        print(help_string)

    # Just a little function to be used with the introduction.
    def __format_keys__(self, keys, split='/', offset=0, max_length=80, prepend=''):
        rtn = ''
        run_length = 0
        for key in keys:
            rtn += prepend + str(key) + split
            run_length += len(str(key))
            if run_length >= max_length:
                run_length = offset
                rtn += '\n' + ' '*offset
        if rtn[-1] == split:
            return rtn[:-1]
        else:
            return rtn

    @property
    def help(self):
        ''' Just a minor function to call help on itself.  Only in here to really help someone get help.'''
        help(self)

    def _repr_pretty_(self, p, cycle):
        self.introduction
        return " "

    def __dir__(self):
        return_list = ['past', 'current', 'future']
        # For the moment, don't expose direct, reweight, or assign, as these are scheme dependent files.
        # They do exist, and always link to the current scheme, however.
        return_list += ['iteration', 'niters', 'scheme', 'list_schemes', 'bin_labels', 'state_labels', 'west', 'trace']
        return sorted(set(return_list))



west = WIPI()
w = west
if __name__ == '__main__':
    # We're gonna print some defaults.
    w.main()
    if w.analysis_mode == False:
        from IPython import embed, embed_kernel
        from IPython.lib.kernel import find_connection_file
        import IPython
        # We're using this to set magic commands.
        # Mostly, we're using it to allow tab completion of objects stored in dictionaries.
        try:
            # Worked on MacOS.  Probably just an older version.
            c = IPython.Config()
        except:
            # Seems to be necessary on Linux, and likely on newer installs.
            c = IPython.terminal.ipapp.load_default_config()
        c.IPCompleter.greedy = True
        embed(banner1='',
             exit_msg='Leaving w_ipa... goodbye.',
             config=c)
    print("")
