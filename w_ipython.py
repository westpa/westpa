import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
import numpy as np
import h5py

# Must be run with the WEST wrapper.
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object
import westpa
import os, sys
import w_assign, w_kinetics, w_kinavg, w_postanalysis_matrix, w_postanalysis_reweight
#warnings.filterwarnings('ignore')
import scipy.sparse as sp

from westtools import (WESTSubcommand, WESTParallelTool, WESTDataReader, WESTDSSynthesizer, BinMappingComponent, 
                       ProgressIndicatorComponent, IterRangeSelection)


class Kinetics(WESTParallelTool):
    '''
        Welcome to w_ipython!
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
        super(Kinetics,self).__init__()
        self.data_reader = WESTDataReader()
        self.wm_env.default_work_manager = self.wm_env.default_parallel_work_manager
        # From assign
        self.dssynth = WESTDSSynthesizer(default_dsname='pcoord')

        # From kinetics...
        self.iter_range = IterRangeSelection() 

        self._iter = 1
        self.config_required = True
        self.version = ".2A"
        global iteration


        #self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.dssynth.add_args(parser)
        self.iter_range.add_args(parser)
        rgroup = parser.add_argument_group('runtime options').add_mutually_exclusive_group()
        rgroup.add_argument('--analysis-only', '-ao', dest='analysis_mode', action='store_true',
                             help='''Use this flag to run the analysis and return to the terminal.''')
        rgroup.add_argument('--reanalyze', '-ra', dest='reanalyze', action='store_true',
                             help='''Use this flag to delete the existing files and reanalyze.''')
        
        parser.set_defaults(compression=True)

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.__config = westpa.rc.config
        self.__settings = self.__config['west']['w_ipython']
        for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
            if (self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None):
                self.scheme = scheme
        with self.data_reader:
            self.dssynth.h5filename = self.data_reader.we_h5filename
            self.dssynth.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)
        self.data_args = args
        self.analysis_mode = args.analysis_mode
        self.reanalyze = args.reanalyze

    def analysis_structure(self):
        #self.settings = self.config['west']['w_ipython']
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
                        analysis_files = ['assign', 'kintrace', 'kinavg', 'flux_matrices', 'kinrw']
                    else:
                        analysis_files = ['assign', 'kintrace', 'kinavg']
                except:
                    analysis_files = ['assign', 'kintrace', 'kinavg']
                    self.__settings['analysis_schemes'][scheme]['postanalysis'] = False
                for name in analysis_files:
                    if self.reanalyze == True:
                        os.remove(os.path.join(path, '{}.h5'.format(name)))
                    try:
                        #print('Loading {} from scheme: {}'.format(name, scheme))
                        self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                        # Try to actually load some data.
                        if name == 'assign':
                            test = self.__analysis_schemes__[scheme][name]['state_labels']
                        if name == 'kintrace':
                            test = self.__analysis_schemes__[scheme][name]['durations']
                        if name =='kinavg':
                            test = self.__analysis_schemes__[scheme][name]['rate_evolution']
                        if name == 'flux_matrices':
                            test = self.__analysis_schemes__[scheme][name]['bin_populations']
                        if name == 'kinrw':
                            test = self.__analysis_schemes__[scheme][name]['rate_evolution']
                    except:
                        self.data_reader.close()
                        print('Unable to load output from {}, or a re-run requested.'.format(name))
                        if name == 'assign':
                        # A lot of this is coded up to avoid the arg parser.  Probably not clean or happy, but it should work for now...
                            assign = w_assign.WAssign()

                            # Taken from w_assign
                            ystates = self.__settings['analysis_schemes'][scheme]['states']
                            states = []
                            for istate, ystate in enumerate(ystates):
                                state = {}
                                state['label'] = ystate.get('label', 'state{}'.format(istate))
                                # coords can be:
                                #  - a scalar, in which case it is one bin, 1-D
                                #  - a single list, which is rejected as ambiguous
                                #  - a list of lists, which is a list of coordinate tuples
                                coords = np.array(ystate['coords'])
                                if coords.ndim == 0:
                                    coords.shape = (1,1)
                                elif coords.ndim == 1:
                                    raise ValueError('list {!r} is ambiguous (list of 1-d coordinates, or single multi-d coordinate?)'
                                                     .format(ystate['coords']))
                                elif coords.ndim > 2:
                                    raise ValueError('coordinates must be 2-D')
                                state['coords'] = coords
                                states.append(state)
                            assign.states = states

                            assign.data_reader = WESTDataReader()
                            assign.data_reader.process_args(self.data_args)
                            assign.output_filename = os.path.join(path, '{}.h5'.format(name))

                            # Taken from bin mapper (core)
                            mapper = getattr(sys.modules['westpa.binning'], self.__settings['analysis_schemes'][scheme]['bins'][0]['type'])
                            if self.__settings['analysis_schemes'][scheme]['bins'][0]['type'] == 'RectilinearBinMapper':
                                boundary_lists = self.__settings['analysis_schemes'][scheme]['bins'][0]['boundaries']
                                for ilist, boundaries in enumerate(boundary_lists):
                                    boundary_lists[ilist] = map((lambda x: 
                                                                   float('inf') 
                                                                   if (x if isinstance(x, basestring) else '').lower() == 'inf' 
                                                                   else x), boundaries)
                            assign.binning.mapper = mapper(boundary_lists)

                            w_assign_config = { 'subsample': False }
                            try:
                                w_assign_config.update(self.__settings['w_assign'])
                            except:
                                pass
                            try:
                                w_assign_config.update(self.__settings['analysis_schemes'][scheme]['w_assign'])
                            except:
                                pass
                            assign.progress.process_args(self.args)
                            assign.subsample = w_assign_config['subsample']
                            assign.work_manager = self.work_manager
                            assign.dssynth = WESTDSSynthesizer(default_dsname='pcoord')
                            assign.dssynth.h5filename = self.data_reader.we_h5filename
                            assign.dssynth.process_args(self.data_args)
                            assign.go()
                            assign.data_reader.close()
                            del(assign)

                            # It closes the h5 file.
                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            #self.data_reader.open()

                        if name == 'kintrace':
                            assignment_file = self.__analysis_schemes__[scheme]['assign']
                            kintrace = w_kinetics.WKinetics()
                            trace = w_kinetics.KinTraceSubcommand(kintrace)
                            w_kinetics_config = { 'correl': False }
                            try:
                                w_kinetics_config.update(self.__settings['w_kinetics'])
                            except:
                                pass
                            try:
                                w_kinetics_config.update(self.__settings['analysis_schemes'][scheme]['w_kinetics'])
                            except:
                                pass
                            trace.progress.process_args(self.args)
                            # Reimplement process_args...
                            # ? This shouldn't be here.
                            trace.correl = w_kinetics_config['correl']
                            trace.assignments_file = assignment_file
                            trace.data_reader = WESTDataReader()
                            trace.data_reader.process_args(self.data_args)
                            trace.iter_range = self.iter_range
                            trace.output_file = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'w', creating_program=True)
                            h5io.stamp_creator_data(trace.output_file)
                            if not trace.iter_range.check_data_iter_range_least(trace.assignments_file):
                                raise ValueError('assignments do not span the requested iterations')
                            self.do_compression = True
                            trace.go()
                            trace.data_reader.close()

                            del(trace)

                            # Open!
                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            # It closes the h5 file.
                            #self.data_reader.open()

                        if name == 'kinavg':
                            ktrace = w_kinavg.WKinAvg()
                            ktrace.work_manager = self.work_manager
                            w_kinavg_config = { 'mcbs_alpha': 0.05, 'mcbs_nsets': 1000, 'evolution': 'cumulative', 'evol_window_frac': 1, 'step_iter': 1, 'bootstrap': True , 'correl': False, 'display_averages': False}
                            try:
                                w_kinavg_config.update(self.__settings['w_kinavg'])
                            except:
                                pass
                            try:
                                w_kinavg_config.update(self.__settings['analysis_schemes'][scheme]['w_kinavg'])
                            except:
                                pass
                            kinavg = w_kinavg.AvgTraceSubcommand(ktrace)
                            kinavg.kinetics_filename = os.path.join(path, '{}.h5'.format('kintrace'))
                            kinavg.assignments_filename = os.path.join(path, '{}.h5'.format('assign'))
                            #kinavg.data_reader = self.data_reader
                            kinavg.data_reader = WESTDataReader()
                            kinavg.data_reader.process_args(self.data_args)
                            kinavg.iter_range = self.iter_range
                            kinavg.mcbs_alpha = w_kinavg_config['mcbs_alpha']
                            kinavg.mcbs_acalpha = kinavg.mcbs_alpha
                            kinavg.mcbs_nsets = w_kinavg_config['mcbs_nsets']
                            kinavg.evolution_mode = w_kinavg_config['evolution']
                            kinavg.evol_window_frac = w_kinavg_config['evol_window_frac']
                            kinavg.iter_range.iter_step = w_kinavg_config['step_iter']
                            kinavg.mcbs_enable = w_kinavg_config['bootstrap']
                            kinavg.correl = w_kinavg_config['correl']
                            kinavg.display_averages = w_kinavg_config['display_averages']
                            with kinavg.data_reader:
                                kinavg.iter_range.process_args(self.args, default_iter_step=None)
                            if kinavg.iter_range.iter_step is None:
                                #use about 10 blocks by default
                                kinavg.iter_range.iter_step = max(1, (kinavg.iter_range.iter_stop - kinavg.iter_range.iter_start) // 10)
                            kinavg.output_filename = os.path.join(path, '{}.h5'.format(name))
                            kinavg.progress.process_args(self.args)
                            if kinavg.evol_window_frac <= 0 or kinavg.evol_window_frac > 1:
                                raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')
                            kinavg.dssynth = WESTDSSynthesizer(default_dsname='pcoord')
                            kinavg.dssynth.h5filename = self.data_reader.we_h5filename
                            kinavg.dssynth.process_args(self.data_args)

                            kinavg.go()
                            kinavg.data_reader.close()
                            del(kinavg)


                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')

                        if name == 'flux_matrices':
                            fmatrix = w_postanalysis_matrix.MatrixRw()
                            fmatrix.work_manager = self.work_manager

                            fmatrix.assignments_file = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format('assign')), 'r')
                            fmatrix.data_reader = WESTDataReader()
                            fmatrix.data_reader.process_args(self.data_args)
                            fmatrix.iter_range = self.iter_range
                            fmatrix.output_file = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'w', creating_program=True)
                            h5io.stamp_creator_data(fmatrix.output_file)
                            fmatrix.dssynth = WESTDSSynthesizer(default_dsname='pcoord')
                            fmatrix.dssynth.h5filename = self.data_reader.we_h5filename
                            fmatrix.dssynth.process_args(self.data_args)

                            matrix_config = { 'sampling_frequency': 'timepoint' }
                            try:
                                matrix_config.update(self.__settings['w_postanalysis_matrix'])
                            except:
                                pass
                            try:
                                matrix_config.update(self.__settings['analysis_schemes'][scheme]['w_postanalysis_matrix'])
                            except:
                                pass
                            fmatrix.progress.process_args(self.args)
                            fmatrix.sampling_frequency = matrix_config['sampling_frequency']


                            fmatrix.go()
                            fmatrix.data_reader.close()
                            del(fmatrix)

                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')

                        if name == 'kinrw':
                            reweight = w_postanalysis_reweight.WPostAnalysisReweightTool()
                            reweight.work_manager = self.work_manager

                            reweight.assignments_filename = os.path.join(path, '{}.h5'.format('assign'))
                            reweight.kinetics_filename = os.path.join(path, '{}.h5'.format('flux_matrices'))
                            reweight.data_reader = WESTDataReader()
                            reweight.data_reader.process_args(self.data_args)
                            reweight.iter_range = self.iter_range
                            reweight.output_filename = os.path.join(path, '{}.h5'.format(name))

                            reweight_config = { 'mcbs_alpha': 0.05, 'mcbs_nsets': 1000, 'evolution': 'cumulative', 'evol_window_frac': 1, 'step_iter': 1, 'bootstrap': True , 'correl': False, 'obs_threshold': 1}
                            try:
                                reweight_config.update(self.__settings['analysis_schemes'][scheme]['w_kinavg'])
                            except:
                                pass
                            try:
                                reweight_config.update(self.__settings['w_postanalysis_reweight'])
                            except:
                                pass
                            try:
                                reweight_config.update(self.__settings['analysis_schemes'][scheme]['w_postanalysis_reweight'])
                            except:
                                pass

                            reweight.progress.process_args(self.args)
                            reweight.mcbs_alpha = reweight_config['mcbs_alpha']
                            reweight.mcbs_acalpha = reweight.mcbs_alpha
                            reweight.mcbs_nsets = reweight_config['mcbs_nsets']
                            reweight.evolution_mode = reweight_config['evolution']
                            reweight.evol_window_frac = reweight_config['evol_window_frac']
                            reweight.iter_range.iter_step = reweight_config['step_iter']
                            reweight.mcbs_enable = reweight_config['bootstrap']
                            reweight.correl = reweight_config['correl']
                            reweight.obs_threshold = reweight_config['obs_threshold']
                            reweight.dssynth = WESTDSSynthesizer(default_dsname='pcoord')
                            reweight.dssynth.h5filename = self.data_reader.we_h5filename
                            reweight.dssynth.process_args(self.data_args)


                            reweight.go()
                            reweight.data_reader.close()
                            del(reweight)

                            self.__analysis_schemes__[scheme][name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                self.work_manager.shutdown()
        print("")
        print("Complete!")






    @property
    def assign(self):
        return self.__analysis_schemes__[self.scheme]['assign']

    @property
    def kinavg(self):
        """
        The output from w_kinavg.py from the current scheme.
        """
        return self.__analysis_schemes__[self.scheme]['kinavg']

    @property
    def kinetics(self):
        return self.__analysis_schemes__[self.scheme]['kintrace']

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
    def kinrw(self):
        if self.__settings['analysis_schemes'][self.scheme]['postanalysis'] == True:
            return self.__analysis_schemes__[self.scheme]['kinrw']
        else:
            value = "This sort of analysis has not been enabled."
            current = { 'bin_prob_evolution': value, 'color_prob_evolution': value, 'conditional_flux_evolution': value, 'rate_evolution': value, 'state_labels': value, 'state_prob_evolution': value }
            return current

    @property
    def matrix(self):
        if self.__settings['analysis_schemes'][self.scheme]['postanalysis'] == True:
            return self.__analysis_schemes__[self.scheme]['flux_matrices']
        else:
            value = "This sort of analysis has not been enabled."
            current = { 'bin_populations': value, 'iterations': value }
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
            w_ipython:
              directory: w_ipython.analysis
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

    @iteration.setter
    def iteration(self, value):
        print("Setting iteration to iter {}.".format(value))
        if value < 0:
            print("Iteration must begin at 1.")
            value = 1
        if value > self.niters:
            print("Cannot go beyond {} iterations!".format(self.niters))
            print("Setting to {}".format(self.niters))
            value = self.niters
        self._iter = value
        self._future = None
        return self._iter

    # Returns the raw values, but can also calculate things based on them.
    class KineticsIteration(dict):
        def __init__(self, kin_h5file, value):
            self.raw = kin_h5file['rate_evolution'][value - 1, :, :]
            self.error = (self.raw['ci_ubound'] - self.raw['ci_lbound']) / (2*self.raw['expected'])
            self.expected = self.raw['expected']
            self.flux = kin_h5file['conditional_flux_evolution'][value - 1, :, :]
            self.ferror = (self.flux['ci_ubound'] - self.flux['ci_lbound']) / (2*self.flux['expected'])
            self.flux = kin_h5file['conditional_flux_evolution'][value - 1, :, :]['expected']
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

    def __get_data_for_iteration__(self, value, seg_ids = None):
        '''
        This returns all important data for the current iteration.  It is optionally
        sorted by the seg_ids, which allows us to match parent-walker child pairs by
        only matching the current iteration.
        This is used internally.
        '''
        iter_group = self.data_reader.get_iter_group(value)
        current = {}
        if seg_ids == None:
            seg_ids = xrange(0, iter_group['seg_index']['weight'].shape[0])
        current['kinavg'] = self.KineticsIteration(self.kinavg, value)
        # Just make these easier to access.
        current['weights'] = iter_group['seg_index']['weight'][seg_ids]
        current['pcoord'] = iter_group['pcoord'][...][seg_ids, :, :]
        try:
            current['auxdata'] = iter_group['auxdata'][...][seg_ids, :, :]
        except:
            pass
        current['parents'] = iter_group['seg_index']['parent_id'][seg_ids]
        current['summary'] = self.data_reader.data_manager.get_iter_summary(value)
        current['seg_id'] = np.array(range(0, iter_group['seg_index'].shape[0]))[seg_ids]
        current['walkers'] = current['summary']['n_particles']
        current['states'] = self.assign['trajlabels'][value-1, :current['walkers'], :][seg_ids]
        current['bins'] = self.assign['assignments'][value-1, :current['walkers'], :][seg_ids]
        # Calculates the bin population for this iteration.
        nbins = self.assign['state_map'].shape[0]
        # We have to take the 'unknown' state into account
        nstates = self.assign['state_labels'].shape[0] + 1
        #current['pop_bins'] = np.histogram(current['bins'].flatten(), bins=range(0, nbins), weights=np.repeat(current['weights'], current['bins'].shape[1]))[0] / current['bins'].shape[1]
        #current['pop_states'] = np.histogram(current['states'].flatten(), bins=range(0, nstates + 1), weights=np.repeat(current['weights'], current['states'].shape[1]))[0] / current['states'].shape[1]
        current['populations'] = self.PopulationsIterations(self.assign, current, self.scheme)
        try:
            # We'll make this not a sparse matrix...
            matrix = self.matrix['iterations/iter_{:08d}'.format(value)]
            # Assume color.
            current['instant_matrix'] = sp.coo_matrix((matrix['flux'][...], (matrix['rows'][...], matrix['cols'][...])), shape=((nbins-1)*2, (nbins-1)*2)).todense()
            current['kinrw'] = self.KineticsIteration(self.kinrw, value)
            current['matrix'] = self.aggregate_matrix[value-1, :, :]
        except:
          # This analysis hasn't been enabled, so we'll simply return the default error message.
            current['instant_matrix'] = self.matrix['bin_populations']
            current['kinrw'] = self.kinrw['rate_evolution']
            #current['rw_error'] = self.kinrw['rate_evolution']
            #current['rw_expected'] = self.kinrw['rate_evolution']
            current['matrix'] = self.matrix['bin_populations']
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
        return self.__get_data_for_iteration__(self.iteration)

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
            return self.__get_data_for_iteration__(self.iteration - 1, self.current['parents'])
        else:
            print("The current iteration is 1; there is no past.")

    def trace(self, seg_id):
        '''
        Runs a trace on a seg_id within the current iteration, all the way back to the beginning,
        returning a dictionary containing all interesting information:

            seg_id, pcoord, states, bins, weights, iteration, auxdata (optional)

        sorted in chronological order.
        '''
        if seg_id >= self.walkers:
            print("Walker seg_id # {} is beyond the max count of {} walkers.".format(seg_id, self.walkers))
            return 1
        current = { 'seg_id': [seg_id], 'pcoord': [self.current['pcoord'][seg_id]], 'states': [self.current['states'][seg_id]], 'weights': [self.current['weights'][seg_id]], 'iteration': [self.iteration], 'bins': [self.current['bins'][seg_id]] }
        try:
            current['auxdata'] = [self.current['auxdata'][seg_id]]
        except:
            pass
        parents = self.current['parents']
        for iter in reversed(range(1, self.iteration)):
            iter_data = self.__get_data_for_iteration__(iter, parents)
            current['pcoord'].append(iter_data['pcoord'][seg_id, :, :])
            current['states'].append(iter_data['states'][seg_id, :])
            current['bins'].append(iter_data['bins'][seg_id, :])
            current['seg_id'].append(iter_data['seg_id'][seg_id])
            current['weights'].append(iter_data['weights'][seg_id])
            try:
                current['auxdata'].append(iter_data['auxdata'][seg_id])
            except:
                pass
            current['iteration'].append(iter)
            seg_id = iter_data['seg_id'][seg_id]
            if seg_id < 0:
                # Necessary for steady state simulations.  This means they started in that iteration.
                break
            parents = self.__get_data_for_iteration__(iter)['parents']
        current['seg_id'] = list(reversed(current['seg_id']))
        current['pcoord'] = np.concatenate(np.array(list(reversed(current['pcoord']))))
        current['states'] = np.concatenate(np.array(list(reversed(current['states']))))
        current['bins'] = np.array(list(reversed(current['states'])))
        current['weights'] = list(reversed(current['weights']))
        current['iteration'] = list(reversed(current['iteration']))
        try:
            current['auxdata'] = list(reversed(current['auxdata']))
        except:
            pass
        try:
            current['auxdata'] = np.array(current['auxdata'])
        except:
            pass
        return current

    @property
    def future(self):
        if self._future == None:
            print("Running child analysis...")
            self.__get_children__()
        return self._future

    def __get_children__(self):
        '''
        Returns all information about the children of a given walker in the current iteration.
        '''
        
        if self.iteration == self.niters:
            print("Currently at iteration {}, which is the max.  There are no children!".format(self.iteration))
            return 0
        iter_data = self.__get_data_for_iteration__(self.iteration+1)
        self._future = { 'kinavg': iter_data['kinavg'], 'weights': [], 'pcoord': [], 'parents': [], 'summary': iter_data['summary'], 'seg_id': [], 'walkers': iter_data['walkers'], 'states': [], 'bins': [] }
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
        matrices = self.__analysis_schemes__[self.scheme]['flux_matrices']
        nbins = matrices['bin_populations'].shape[1]
        self.__analysis_schemes__[self.scheme]['aggregate_matrix'] = np.zeros((self.niters, nbins, nbins))
        self.__analysis_schemes__[self.scheme]['total_pop'] = np.zeros((self.niters, nbins))
        matrix_accumulator = w_postanalysis_reweight.accumulate_statistics
        normalize = w_postanalysis_reweight.normalize
        total_fluxes = np.zeros((nbins, nbins))
        total_obs = np.zeros((nbins, nbins))
        for iter in range(2, self.niters+1):
            total_fluxes, total_obs, total_pop = matrix_accumulator(self.matrix, iter-1, iter, nbins, total_fluxes, total_obs)
            self.__analysis_schemes__[self.scheme]['aggregate_matrix'][iter-1][:] = normalize(total_fluxes)

    def go(self):
        self.data_reader.open()
        self.analysis_structure()
        self.data_reader.open()
        self.niters = self.kinavg['rate_evolution']['expected'].shape[0]
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


west = Kinetics()
w = west
if __name__ == '__main__':
    # We're gonna print some defaults.
    print("")
    print("Welcome to w_ipython v. {}!".format(w.version))
    print("Run w.introduction for a more thorough introduction, or w.help to see a list of options.")
    print("Running analysis & loading files.")
    w.main()
    print('Your current scheme, system and iteration are : {}, {}, {}'.format(w.scheme, os.getcwd(), w.iteration))
    if w.analysis_mode == False:
        from IPython import embed
        embed(banner1='',
             exit_msg='Leaving w_ipython... goodbye.')
    print("")
