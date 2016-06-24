import numpy as np
import h5py

# Must be run with the WEST wrapper.
from westpa import h5io
from westpa.h5io import WESTPAH5File
from westpa.extloader import get_object
import westpa
import os, sys
import w_assign, w_kinetics, w_kinavg
import warnings
warnings.filterwarnings('ignore')

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

        All information about the current iteration is available in a dictionary called 'current'.

            w.current.keys():
            walkers, summary, states, seg_id, weights, parents, kinavg, pcoord, bins, and auxdata, if it exists.

        kinavg, states, and bins are pulled from the output from w_kinavg and w_assign; they always correspond to
        what is used in the current analysis scheme.  If you change the scheme, those, too, will change.

        You can look at the information for any walker by simply indexing according to that seg_id.

        Information about the previous iteration is available in the past dictionary, which contains the same information.
        It is keyed to use the current iteration's seg_id, such that if you're looking at walker 0 in the current iteration,
        w.past['pcoord'][0] will give you the progress coordinate for the parent of walker 0.  You can look at the actual
        walker seg_id in the previous iteration by
        
            w.past['parents'][0]

        The kinavg, assign, and kinetics file from the current state are available for raw access from:

            w.kinavg, w.assign, and w.kinetics

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
        global iteration


        #self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.dssynth.add_args(parser)
        self.iter_range.add_args(parser)
        
        parser.set_defaults(compression=True)

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.__config = westpa.rc.config
        self.__settings = self.__config['west']['w_ipython']
        for ischeme, scheme in enumerate(self.__settings['analysis_schemes']):
            if (self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None):
                self._scheme = scheme
                break
        with self.data_reader:
            self.dssynth.h5filename = self.data_reader.we_h5filename
            self.dssynth.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)
        print("Arguments Processed!")

    def analysis_structure(self, rerun=False):
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
                path = os.path.join(os.getcwd(), self.__settings['directory'], scheme)
                try:
                    os.mkdir(path)
                except:
                    pass
                tmp = {}
                for name in ['assign', 'kintrace', 'kinavg']:
                    try:
                        if rerun == False:
                            print(os.path.join(path, '{}.h5'.format(name)), 'r')
                            tmp[name] = h5py.File(os.path.join(path, '{}.h5'.format(name)), 'r')
                        else:
                            raise Exception("Rerun is set to true, or the output cannot be loaded.")
                    except:
                        print('Unable to load output from {}'.format(name))
                        print('Rerunning...')
                        if name == 'assign':
                        #    # For the moment, let's hardcode options.
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

                            assign.data_reader = self.data_reader
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

                            # Fucking thing.  Just to get it to stop complaining...
                            assign.progress.process_args(self.args)
                            assign.work_manager = self.work_manager
                            assign.dssynth = self.dssynth
                            assign.go()
                            # It closes the h5 file.
                            tmp[name] = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            self.data_reader.open()

                        if name == 'kintrace':
                            assign = tmp['assign']
                            kintrace = w_kinetics.WKinetics()
                            trace = w_kinetics.KinTraceSubcommand(kintrace)
                            trace.progress.process_args(self.args)
                            # Reimplement process_args...
                            trace.assignments_file = assign
                            trace.data_reader = self.data_reader
                            trace.iter_range = self.iter_range
                            trace.output_file = h5io.WESTPAH5File(os.path.join(path, '{}.h5'.format(name)), 'w', creating_program=True)
                            h5io.stamp_creator_data(trace.output_file)
                            if not trace.iter_range.check_data_iter_range_least(trace.assignments_file):
                                raise ValueError('assignments do not span the requested iterations')
                            self.do_compression = True
                            trace.go()



                            # Open!
                            tmp[name] = h5py.File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            # It closes the h5 file.
                            self.data_reader.open()

                        if name == 'kinavg':
                            ktrace = w_kinavg.WKinAvg()
                            ktrace.work_manager = self.work_manager
                            w_kinavg_config = { 'mcbs_alpha': 0.05, 'mcbs_nsets': 1000, 'evolution': 'cumulative', 'evol_window_frac': 1, 'step_iter': 1, 'bootstrap': True }
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
                            kinavg.data_reader = self.data_reader
                            kinavg.iter_range = self.iter_range
                            kinavg.mcbs_alpha = w_kinavg_config['mcbs_alpha']
                            kinavg.mcbs_acalpha = kinavg.mcbs_alpha
                            kinavg.mcbs_nsets = w_kinavg_config['mcbs_nsets']
                            kinavg.evolution_mode = w_kinavg_config['evolution']
                            kinavg.evol_window_frac = w_kinavg_config['evol_window_frac']
                            kinavg.iter_range.iter_step = w_kinavg_config['step_iter']
                            kinavg.bootstrap = w_kinavg_config['bootstrap']
                            with kinavg.data_reader:
                                kinavg.iter_range.process_args(self.args, default_iter_step=None)
                            if kinavg.iter_range.iter_step is None:
                                #use about 10 blocks by default
                                kinavg.iter_range.iter_step = max(1, (kinavg.iter_range.iter_stop - kinavg.iter_range.iter_start) // 10)
                            kinavg.output_filename = os.path.join(path, '{}.h5'.format(name))
                            kinavg.progress.process_args(self.args)
                            if kinavg.evol_window_frac <= 0 or kinavg.evol_window_frac > 1:
                                raise ValueError('Parameter error -- fractional window defined by --window-frac must be in (0,1]')
                            kinavg.dssynth = self.dssynth

                            kinavg.go()


                            tmp[name] = h5py.File(os.path.join(path, '{}.h5'.format(name)), 'r')
                            # It closes the h5 file.
                            self.data_reader.open()
                self.__analysis_schemes__[scheme] = tmp





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

    #@property
    #def west(self):
    #    return self.

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
        if self.__settings['analysis_schemes'][scheme]['enabled'] == True or self.__settings['analysis_schemes'][scheme]['enabled'] == None:
            self._scheme = scheme
        else:
            print("Scheme cannot be changed to scheme: {}; it is not enabled!".format(scheme))

    #def enable_scheme(self, scheme):
    #    self.__settings['analysis_schemes'][scheme]['enabled'] = True
    #    self.analysis_structure()

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
        for scheme in self.__settings['analysis_schemes']:
            print(scheme)

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

    @iteration.setter
    def iteration(self, value):
        print("Setting iteration to iter {}.".format(value))
        if value < 0:
            raise ValueError("Iteration must begin at 1.")
        if value > self.niters:
            print("Cannot go beyond {} iterations!".format(self.niters))
            print("Setting to {}".format(self.niters))
            value = self.niters
        self._iter = value
        self._future = None
        return self._iter

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
        current['kinavg'] = self.kinavg['rate_evolution'][value - 1, :, :]
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
        current['states'] = self.assign['trajlabels'][value-1, :current['walkers'], :]
        current['bins'] = self.assign['assignments'][value-1, :current['walkers'], :]
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
        
        iter_data = self.__get_data_for_iteration__(self.iteration+1)
        self._future = { 'kinavg': iter_data['kinavg'], 'weights': [], 'pcoord': [], 'parents': [], 'summary': iter_data['summary'], 'seg_id': [], 'walkers': iter_data['walkers'], 'states': [], 'bins': [] }
        for seg_id in range(0, self.walkers):
            if self.iteration == self.niters:
                return 0
            else:
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




    def go(self):
        self.data_reader.open()
        self.analysis_structure()
        self.niters = self.kinavg['rate_evolution']['expected'].shape[0]
        self.iteration = 1
        print(self.__doc__)

west = Kinetics()
w = west
if __name__ == '__main__':
    w.main()


