import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=FutureWarning)
import numpy as np
import h5py

# Must be run with the WEST wrapper.
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object
import westpa
import os, sys
import w_assign, w_direct, w_reweight
#warnings.filterwarnings('ignore')
import scipy.sparse as sp

from westtools import (WESTSubcommand, WESTParallelTool, WESTDataReader, WESTDSSynthesizer, BinMappingComponent, 
                       ProgressIndicatorComponent, IterRangeSelection, Plotter)


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
            w.current_scheme = OUTPUT FROM w.list_schemes

        To see the states and bins defined in the current analysis scheme:

            w.states
            w.bin_labels

        All information about the current iteration is available in a dictionary called 'current'.

            w.current.keys():
            walkers, summary, states, seg_id, weights, parents, kinavg, pcoord, bins, populations, and auxdata, if it exists.

        Populations prints the bin and state populations calculated by w_assign; it contains the following attributes:

            states, bins

            which can be called as w.current['populations'].states to return a numpy object.

        If postanalysis has been run, the following information is also available:

            instant_matrix, matrix (the aggregate matrix), kinrw

        Both the kinrw and kinavg key in 'current' have the following attributes:

            expected, error, flux, ferror, raw

            where raw is returned on a basic call.

        kinavg, states, and bins are pulled from the output from w_kinavg and w_assign; they always correspond to
        what is used in the current analysis scheme.  If you change the scheme, those, too, will change.

        You can look at the information for any walker by simply indexing according to that seg_id.

        Information about the previous iteration is available in the past dictionary, which contains the same information.
        It is keyed to use the current iteration's seg_id, such that if you're looking at walker 0 in the current iteration,
        w.past['pcoord'][0] will give you the progress coordinate for the parent of walker 0.  You can look at the actual
        walker seg_id in the previous iteration by
        
            w.past['parents'][0]

        The kinavg, assign, and kinetics file from the current state are available for raw access.  The postanalysis output
        is also available, should it exist:

            w.kinavg, w.assign, w.kinetics, w.matrix, and w.kinrw

        In addition, the function w.trace(seg_id) will run a trace over a seg_id in the current iteration and return a dictionary
        containing all pertinent information about that seg_id's history.  It's best to store this, as the trace can be expensive.

        Run help on any function or property for more information!

        Happy analyzing!
                
    '''

    def __init__(self):
        super(WIPI,self).__init__()
        self.data_reader = WESTDataReader()
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager

        self._iter = 1
        self.config_required = True
        self.version = ".99A"
        # Set to matplotlib if you want that.  But why would you?
        # Well, whatever, we'll just set it to that for now.
        self.interface = 'matplotlib'
        global iteration

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        rgroup = parser.add_argument_group('runtime options')
        rgroup.add_argument('--analysis-only', '-ao', dest='analysis_mode', action='store_true',
                             help='''Use this flag to run the analysis and return to the terminal.''')
        rgroup.add_argument('--reanalyze', '-ra', dest='reanalyze', action='store_true',
                             help='''Use this flag to delete the existing files and reanalyze.''')
        rgroup.add_argument('--terminal', '-t', dest='plotting', action='store_true',
                             help='''Plot output in terminal.''')
        # There is almost certainly a better way to handle this, but we'll sort that later.
        rgroup.add_argument('--f', '-f', dest='extra', default='blah',
                             help='''Temporary holding place for when this is called in a Jupyter notebook.''')
        
        parser.set_defaults(compression=True)

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.__config = westpa.rc.config
        self.__settings = self.__config['west']['analysis']
        for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
            if (self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None):
                self.scheme = scheme
        self.data_args = args
        self.analysis_mode = args.analysis_mode
        self.reanalyze = args.reanalyze
        if args.plotting:
            self.interface = 'text'

    def analysis_structure(self):
        # Make sure everything exists.
        try:
            os.mkdir(self.__settings['directory'])
        except:
            pass
        # Now, check to see whether they exist, and then load them.
        self.__analysis_schemes__ = {}
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
                for name in analysis_files:
                    if self.reanalyze == True:
                        try:
                            os.remove(os.path.join(path, '{}.h5'.format(name)))
                        except:
                            pass
                    try:
                        if self.reanalyze == True:
                            raise ValueError('Reanalyze set to true.')
                        self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                        # Try to actually load some data.
                        if name == 'assign':
                            test = self.__analysis_schemes__[scheme][name]['state_labels']
                        if name == 'direct':
                            # Comes from the kinetics analysis (old w_kinetics).  If we don't have this, we need to rerun it.
                            test = self.__analysis_schemes__[scheme][name]['durations']
                        if name == 'reweight':
                            # Comes from the flux matrix analysis.  Same as above.
                            test = self.__analysis_schemes__[scheme][name]['iterations']
                    except:
                        #self.data_reader.close()
                        print('Reanalyzing file {}.h5 for scheme {}.'.format(name, scheme))
                        try:
                            os.remove(os.path.join(path, '{}.h5'.format(name)))
                        except:
                            pass
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
                            for key,value in w_assign_config.iteritems():
                                args.append(str('--') + str(key))
                                args.append(str(value))
                            # We're just calling the built in function.
                            # This is a lot cleaner than what we had in before, and far more workable.
                            args.append('--config-from-file')
                            args.append('--scheme-name')
                            args.append('{}'.format(scheme))
                            assign.make_parser_and_process(args=args)
                            # We want to use the work manager we have here.  Otherwise, just let the tool sort out what it needs, honestly.
                            assign.work_manager = self.work_manager

                            assign.go()
                            assign.data_reader.close()
                            del(assign)

                            # It closes the h5 file.
                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')

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
                            for key,value in analysis_config.iteritems():
                                if key != 'extra':
                                    args.append(str('--') + str(key).replace('_', '-'))
                                    args.append(str(value))
                            # This is for stuff like disabling correlation analysis, etc.
                            if 'extra' in analysis_config.keys():
                                for value in analysis_config['extra']:
                                    args.append(str('--') + str(value).replace('_', '-'))
                            # We want to not display the averages, so...
                            args.append('--disable-averages')
                            analysis.make_parser_and_process(args=args)
                            # We want to hook into the existing work manager.
                            analysis.work_manager = self.work_manager

                            analysis.go()
                            del(analysis)

                            # Open!
                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')

        # Make sure this doesn't get too far out, here.  We need to keep it alive as long as we're actually analyzing things.
        self.work_manager.shutdown()
        print("")
        print("Complete!")

    @property
    def assign(self):
        return self.__analysis_schemes__[self.scheme]['assign']

    @property
    def direct(self):
        """
        The output from w_kinavg.py from the current scheme.
        """
        return self.__analysis_schemes__[self.scheme]['direct']

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
        if self.__settings['analysis_schemes'][self.scheme]['postanalysis'] == True:
            return self.__analysis_schemes__[self.scheme]['reweight']
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
        return self._scheme

    @scheme.setter
    def scheme(self, scheme):
        self._future = None
        if scheme in self.__settings['analysis_schemes']:
            pass
        else:
            for ischeme, schemename in enumerate(self.__settings['analysis_schemes']):
                if ischeme == scheme:
                    scheme = schemename
        if self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None:
            self._scheme = scheme
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
        self._iter = value
        self._future = None
        self._current = None
        self._past = None
        return self._iter

    @property
    def walkers(self):
        '''
        The number of walkers active in the current iteration.
        '''
        # Returns number of walkers for iteration X.  Assumes current iteration, but can go with different one.
        return self.current['summary']['n_particles']

    @property
    def aggregate_walkers(self):
        return self.west['summary']['n_particles'][:self.iteration].sum()

    # Returns the raw values, but can also calculate things based on them.
    class KineticsIteration(dict):
        def __init__(self, kin_h5file, index):
            # Check the start and stop, calculate the block size, and index appropriately.
            # While we could try and automatically generate this above, it's a little more consistent to try it here.
            # This should show the first block for which the current iteration has contributed data.
            self.step_iter = (kin_h5file['rate_evolution']['iter_stop'][0] - kin_h5file['rate_evolution']['iter_start'][0])[1,0]
            value = ((index-2) // self.step_iter)
            if value < 0:
                value = 0
            self.raw = kin_h5file['rate_evolution'][value, :, :]
            self.error = (self.raw['ci_ubound'] - self.raw['ci_lbound']) / (2*self.raw['expected'])
            self.expected = self.raw['expected']
            self.flux = kin_h5file['conditional_flux_evolution'][value, :, :]
            self.ferror = (self.flux['ci_ubound'] - self.flux['ci_lbound']) / (2*self.flux['expected'])
            self.flux = kin_h5file['conditional_flux_evolution'][value, :, :]['expected']
            self.__dict__ = { 'raw': self.raw, 'error': self.error, 'expected': self.expected, 'flux': self.flux, 'ferror': self.ferror }
        def __getattr__(self, attr):
            return self.__dict__[attr]
        def __repr__(self):
            return repr(self.raw)
        def __getitem__(self, value):
            return self.raw[value]
        def keys(self):
            a = []
            for i in self.raw.dtype.names:
                a.append(i)
            return a
        @property
        def names(self):
            return self.keys()

    # Returns the raw values, but can also calculate things based on them.
    class CalcPopIteration(dict):
        def __init__(self, sp_h5file, index):
            self.step_iter = (sp_h5file['state_pop_evolution']['iter_stop'][0] - sp_h5file['state_pop_evolution']['iter_start'][0])[1]
            value = ((index-1) // self.step_iter)
            if value < 0:
                value = 0
            self.raw = sp_h5file['state_pop_evolution'][value, :]
            self.error = (self.raw['ci_ubound'] - self.raw['ci_lbound']) / (2*self.raw['expected'])
            self.expected = self.raw['expected']
            self.__dict__ = { 'raw': self.raw, 'error': self.error, 'expected': self.expected }
        def __getattr__(self, attr):
            return self.__dict__[attr]
        def __repr__(self):
            return repr(self.raw)
        def __getitem__(self, value):
            return self.raw[value]
        def keys(self):
            a = []
            for i in self.raw.dtype.names:
                a.append(i)
            return a
        @property
        def names(self):
            return self.keys()

    class PopulationsIterations():
        def __init__(self, assign, current, scheme):
            nbins = assign['state_map'].shape[0]
            # We have to take the 'unknown' state into account
            nstates = assign['state_labels'].shape[0] + 1
            self.bins = np.histogram(current['bins'].flatten(), bins=range(0, nbins), weights=np.repeat(current['weights'], current['bins'].shape[1]))[0] / current['bins'].shape[1]
            self.states = np.histogram(current['states'].flatten(), bins=range(0, nstates + 1), weights=np.repeat(current['weights'], current['states'].shape[1]))[0] / current['states'].shape[1]
            self.scheme = scheme
        def __repr__(self):
            print("The following are populations from the assign.h5 file from scheme: ".format(self.scheme))
            print("")
            print("Bin Populations:")
            print(self.bins)
            print("")
            print("State Populations:")
            print(self.states)
            print("")
            print("Use the properties .states and .bins to access these.")
            #return repr((self.states, self.bins))
            return (" ")

    class __get_data_for_iteration__():
        def __init__(self, parent, value, seg_ids = None):
            '''
            This returns all important data for the current iteration.  It is optionally
            sorted by the seg_ids, which allows us to match parent-walker child pairs by
            only matching the current iteration.
            This is used internally.
            '''
            # We've classed this so that we can override some of the normal functions and allow indexing via seg_id
            iter_group = parent.data_reader.get_iter_group(value)
            current = {}
            if seg_ids == None:
                seg_ids = xrange(0, iter_group['seg_index']['weight'].shape[0])
            current['kinavg'] = parent.KineticsIteration(parent.direct, value)
            current['statepops'] = parent.CalcPopIteration(parent.direct, value)
            # Just make these easier to access.
            current['weights'] = iter_group['seg_index']['weight'][seg_ids]
            current['pcoord'] = iter_group['pcoord'][...][seg_ids, :, :]
            try:
                current['auxdata'] = {}
                for key in iter_group['auxdata'].keys():
                    current['auxdata'][key] = iter_group['auxdata'][key][...][seg_ids, :]
            except:
                pass
            current['parents'] = iter_group['seg_index']['parent_id'][seg_ids]
            current['summary'] = parent.data_reader.data_manager.get_iter_summary(int(value))
            current['seg_id'] = np.array(range(0, iter_group['seg_index'].shape[0]))[seg_ids]
            current['walkers'] = current['summary']['n_particles']
            current['states'] = parent.assign['trajlabels'][value-1, :current['walkers'], :][seg_ids]
            current['bins'] = parent.assign['assignments'][value-1, :current['walkers'], :][seg_ids]
            # Calculates the bin population for this iteration.
            nbins = parent.assign['state_map'].shape[0]
            # We have to take the 'unknown' state into account
            nstates = parent.assign['state_labels'].shape[0] + 1
            #current['pop_bins'] = np.histogram(current['bins'].flatten(), bins=range(0, nbins), weights=np.repeat(current['weights'], current['bins'].shape[1]))[0] / current['bins'].shape[1]
            #current['pop_states'] = np.histogram(current['states'].flatten(), bins=range(0, nstates + 1), weights=np.repeat(current['weights'], current['states'].shape[1]))[0] / current['states'].shape[1]
            current['populations'] = parent.PopulationsIterations(parent.assign, current, parent.scheme)
            current['plot'] = Plotter(parent.direct, parent.reweight, parent.iteration, parent.assign['bin_labels'], parent.assign['state_labels'], current['populations'].states, current['populations'].bins, parent.interface)
            try:
                # We'll make this not a sparse matrix...
                matrix = parent.reweight['iterations/iter_{:08d}'.format(value)]
                # Assume color.
                current['instant_matrix'] = sp.coo_matrix((matrix['flux'][...], (matrix['rows'][...], matrix['cols'][...])), shape=((nbins-1)*2, (nbins-1)*2)).todense()
                current['kinrw'] = parent.KineticsIteration(parent.reweight, value)
                # This feature is borked, as we've removed the code.  Ergo...
                #current['matrix'] = self.aggregate_matrix[value-1, :, :]
                current['rwstatepops'] = parent.CalcPopIteration(parent.reweight, value)
            except:
              # This analysis hasn't been enabled, so we'll simply return the default error message.
                current['instant_matrix'] = parent.reweight['bin_populations']
                current['kinrw'] = parent.reweight['rate_evolution']
                current['rwstatepops'] = parent.reweight['rate_evolution']
                current['matrix'] = parent.reweight['bin_populations']
            self.raw = current
        def __repr__(self):
            return repr(self.raw)
        def keys(self):
            return self.raw.keys()
        def __getitem__(self, value):
            active_items = ['kinavg', 'statepops', 'weights', 'pcoord', 'auxdata', 'parents', 'summary', 'seg_id', 'walkers', 'states', 'bins', 'populations', 'plot', 'instant_matrix', 'kinrw', 'matrix', 'rwstatepops']
            if value in active_items:
                return self.raw[value]
            else:
                current = {}
                seg_items = ['weights', 'pcoord', 'auxdata', 'parents', 'seg_id', 'states']
                #for i in seg_items:
                #    current[i] = self.raw[i]
                current['pcoord'] = self.raw['pcoord'][value, :, :]
                current['states'] = self.raw['states'][value, :]
                current['bins'] = self.raw['bins'][value, :]
                current['parents'] = self.raw['parents'][value]
                current['seg_id'] = self.raw['seg_id'][value]
                current['weights'] = self.raw['weights'][value]
                try:
                    current['auxdata'] = {}
                    for key in self.raw['auxdata'].keys():
                        current['auxdata'][key] = self.raw['auxdata'][key][value]
                except:
                    pass
                return current

    @property
    def current(self):
        '''
        All interesting data from the current iteration/scheme.  Whenever you change the scheme or iteration,
        this dictionary is automatically updated.
        Contains the following keys:

            kinavg, weights, pcoord, auxdata (optional), parents, summary, seg_id, walkers, states, bins

        kinavg, states, and bins refer to the output from w_kinavg and w_assign for this iteration
        and analysis scheme.  They are NOT dynamics bins, necessarily, but the bins defined in
        west.cfg.  If you change the analysis scheme, so, too, will the important values.
        '''
        if self._current == None:
            self._current = self.__get_data_for_iteration__(value=self.iteration, parent=self)
            return self._current
        else:
            return self._current

    @property
    def past(self):
        '''
        Returns the same information as current, but indexed according to the current iteration's seg_ids.
        That is, if w.iteration = 2, and you're interested in walker 0 in iteration 2, indexing any of the
        dictionary values here will return information about the PARENT of walker 0 in iteration 2.

        You can find the real seg_id by looking at the seg_id index.
        In addition,
            w.past['seg_id'][0] = w.current['parents'][0]

        by construction.
        '''
        if self.iteration > 1:
            if self._past == None:
                self._past = self.__get_data_for_iteration__(value=self.iteration - 1, seg_ids=self.current['parents'], parent=self)
                return self._past
            else:
                return self._past
        else:
            print("The current iteration is 1; there is no past.")


    def trace(self, seg_id):
        '''
        Runs a trace on a seg_id within the current iteration, all the way back to the beginning,
        returning a dictionary containing all interesting information:

            seg_id, pcoord, states, bins, weights, iteration, auxdata (optional)

        sorted in chronological order.
        '''
        # It should be noted that this is not a fast function, but was designed more as a 'proof of principle' of the generality of this approach.
        # It could, and most certainly should, have its speed increased.
        if seg_id >= self.walkers:
            print("Walker seg_id # {} is beyond the max count of {} walkers.".format(seg_id, self.walkers))
            return 1
        current = { 'seg_id': [seg_id], 'pcoord': [self.current['pcoord'][seg_id]], 'states': [self.current['states'][seg_id]], 'weights': [self.current['weights'][seg_id]], 'iteration': [self.iteration], 'bins': [self.current['bins'][seg_id]] }
        try:
            current['auxdata'] = {}
            for key in self.current['auxdata'].keys():
                current['auxdata'][key] = [self.current['auxdata'][key][seg_id]]
        except:
            pass
        parents = self.current['parents']
        for iter in reversed(range(1, self.iteration)):
            iter_data = self.__get_data_for_iteration__(value=iter, seg_ids=parents, parent=self)
            current['pcoord'].append(iter_data['pcoord'][seg_id, :, :])
            current['states'].append(iter_data['states'][seg_id, :])
            current['bins'].append(iter_data['bins'][seg_id, :])
            current['seg_id'].append(iter_data['seg_id'][seg_id])
            current['weights'].append(iter_data['weights'][seg_id])
            try:
                for key in self.current['auxdata'].keys():
                    current['auxdata'][key].append(iter_data['auxdata'][key][seg_id])
            except:
                pass
            current['iteration'].append(iter)
            seg_id = iter_data['seg_id'][seg_id]
            if seg_id < 0:
                # Necessary for steady state simulations.  This means they started in that iteration.
                break
            parents = self.__get_data_for_iteration__(value=iter, parent=self)['parents']
        current['seg_id'] = list(reversed(current['seg_id']))
        current['pcoord'] = np.concatenate(np.array(list(reversed(current['pcoord']))))
        current['states'] = np.concatenate(np.array(list(reversed(current['states']))))
        current['bins'] = np.concatenate(np.array(list(reversed(current['bins']))))
        current['weights'] = list(reversed(current['weights']))
        current['iteration'] = list(reversed(current['iteration']))
        try:
            for key in self.current['auxdata'].keys():
                current['auxdata'][key] = np.concatenate(np.array(list(reversed(current['auxdata'][key]))))
        except:
            pass
        #try:
        #    current['auxdata'] = np.array(current['auxdata'])
        #except:
        #    pass
        return current

    @property
    def future(self):
        if self._future == None:
            print("Running child analysis...")
            self.__get_children__()
        return self._future

    class Future():
        def __init__(self, rep={}):
            self.raw = rep
        def __repr__(self):
            return repr(self.raw)
        def keys(self):
            return self.raw.keys()
        def __getitem__(self, value):
            active_items = ['kinavg', 'statepops', 'weights', 'pcoord', 'auxdata', 'parents', 'summary', 'seg_id', 'walkers', 'states', 'bins']
            if value in active_items:
                return self.raw[value]
            else:
                current = {}
                seg_items = ['weights', 'pcoord', 'auxdata', 'parents', 'seg_id', 'states']
                current['pcoord'] = self.raw['pcoord'][value]
                current['states'] = self.raw['states'][value]
                current['bins'] = self.raw['bins'][value]
                current['seg_id'] = self.raw['seg_id'][value]
                current['weights'] = self.raw['weights'][value]
                current['parents'] = self.raw['parents'][value]
                try:
                    for key in self.raw['auxdata'].keys():
                        current['auxdata'][key] = self.raw['auxdata'][key][value]
                except:
                    pass
                return current

    def __get_children__(self):
        '''
        Returns all information about the children of a given walker in the current iteration.
        '''
        
        if self.iteration == self.niters:
            print("Currently at iteration {}, which is the max.  There are no children!".format(self.iteration))
            return 0
        iter_data = self.__get_data_for_iteration__(value=self.iteration+1, parent=self)
        self._future = self.Future(rep={ 'kinavg': iter_data['kinavg'], 'weights': [], 'pcoord': [], 'parents': [], 'summary': iter_data['summary'], 'seg_id': [], 'walkers': iter_data['walkers'], 'states': [], 'bins': [] })
        #self._future = { 'kinavg': iter_data['kinavg'], 'weights': [], 'pcoord': [], 'parents': [], 'summary': iter_data['summary'], 'seg_id': [], 'walkers': iter_data['walkers'], 'states': [], 'bins': [] }
        for seg_id in range(0, self.walkers):
            children = np.where(iter_data['parents'] == seg_id)[0]
            if len(children) == 0:
                error = "No children for seg_id {}.".format(seg_id)
                self._future['weights'].append(error)
                self._future['pcoord'].append(error)
                self._future['parents'].append(error)
                self._future['seg_id'].append(error)
                self._future['states'].append(error)
                self._future['bins'].append(error)
            else:
                # Now, we're gonna put them in the thing.
                value = self.iteration+1 
                self._future['weights'].append(iter_data['weights'][children])
                self._future['pcoord'].append(iter_data['pcoord'][...][children, :, :])
                try:
                    aux_data = iter_data['auxdata'][...][children, :, :]
                    try:
                        current['aux_data'].append(aux_data)
                    except:
                        current['aux_data'] = aux_data
                except:
                    pass
                self._future['parents'].append(iter_data['parents'][children])
                self._future['seg_id'].append(iter_data['seg_id'][children])
                self._future['states'].append(self.assign['trajlabels'][value-1, children, :])
                self._future['bins'].append(self.assign['assignments'][value-1, children, :])


    @property
    def aggregate_matrix(self):
        if self.__settings['analysis_schemes'][self.scheme]['postanalysis'] == True:
            try:
                if self.__analysis_schemes__[self.scheme]['aggregate_matrix'] == None:
                    print("Calculating aggregate matrix...")
                    self.__add_matrix__()
                return self.__analysis_schemes__[self.scheme]['aggregate_matrix']
            except:
                print("Calculating aggregate matrix...")
                self.__add_matrix__()
                return self.__analysis_schemes__[self.scheme]['aggregate_matrix']

    def __add_matrix__(self):
        # Hooks into the existing tools in w_postanalysis_reweight
        matrices = self.__analysis_schemes__[self.scheme]['reweight']
        nbins = matrices['bin_populations'].shape[1]
        self.__analysis_schemes__[self.scheme]['aggregate_matrix'] = np.zeros((self.niters, nbins, nbins))
        self.__analysis_schemes__[self.scheme]['total_pop'] = np.zeros((self.niters, nbins))
        print("Hey, this analysis is borked.")
        #matrix_accumulator = w_postanalysis_reweight.accumulate_statistics
        #normalize = w_postanalysis_reweight.normalize
        #total_fluxes = np.zeros((nbins, nbins))
        #total_obs = np.zeros((nbins, nbins))
        #for iter in range(2, self.niters+1):
        #    total_fluxes, total_obs, total_pop = matrix_accumulator(self.reweight, iter-1, iter, nbins, total_fluxes, total_obs)
        #    self.__analysis_schemes__[self.scheme]['aggregate_matrix'][iter-1][:] = normalize(total_fluxes)

    def go(self):
        self.data_reader.open()
        self.analysis_structure()
        # Seems to be consistent with other tools, such as w_assign.  For setting the iterations.
        self.niters = self.data_reader.current_iteration - 1
        self.iteration = 1
        if self.__settings['analysis_schemes'][self.scheme]['postanalysis'] == True:
            self.__analysis_schemes__[self.scheme]['aggregate_matrix'] = None

    def _help(self, item=None):
        if item == None:
            print(self.__doc__)
        else:
            print(item.__doc__)

    @property
    def introduction(self):
        self._help()

    @property
    def help(self):
        help_string = '''
        Call as a dictionary item, unless item is a .property; then simply call on the item itself

        w.past, w.current, w.future:
            
            weights, pcoord, seg_id, parents, auxdata, summary, walkers, states, bins, matrix, instant_matrix

                matrix          - aggregate transition matrix.
                instant_matrix  - instant transition matrix (uses current iteration only)
                bins            - bin assignments for walkers from current assignment file
                states          - state assignments for walkers from current assignment file

            kinavg, kinrw - call as is for native dataset, or:

                .expected, .error, .raw, .flux, .ferror
                expected, ci_ubound, ci_lbound, sterr, corrlen

            population.states, population.bin

        w.iteration     - Get/set current iteration
        w.niters        - Maximum number of iterations
        w.scheme        - Get/set current analysis scheme
        w.list_schemes  - Lists all analysis schemes, and current
        w.bin_labels    - pcoord values for bin assignments from current assignment file
        w.state_labels  - state labels for states from current assignment file

        The following give raw access to the h5 files associated with the current scheme

        w.kinavg
        w.kintrace
        w.assign
        w.west
        w.kinrw
        w.matrix

        w.trace()
        '''
        print(help_string)


west = WIPI()
w = west
if __name__ == '__main__':
    # We're gonna print some defaults.
    print("")
    print("Welcome to w_ipa (WESTPA Interactive Python Analysis) v. {}!".format(w.version))
    print("Run w.introduction for a more thorough introduction, or w.help to see a list of options.")
    print("Running analysis & loading files.")
    w.main()
    print('Your current scheme, system and iteration are : {}, {}, {}'.format(w.scheme, os.getcwd(), w.iteration))
    if w.analysis_mode == False:
        from IPython import embed, embed_kernel
        from IPython.lib.kernel import find_connection_file
        import IPython
        embed(banner1='',
             exit_msg='Leaving w_ipa... goodbye.')
        #cf=find_connection_file("*")
        #kern = w.work_manager.submit(embed_kernel())
        #print("We are running!")
        #print(dir(kern))
    print("")
