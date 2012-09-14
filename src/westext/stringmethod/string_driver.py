from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np
import types

import west
from west.util import extloader
from west.util.miscfn import vgetattr
from westext.stringmethod import WESTStringMethod, DefaultStringMethod
from west.pcoords import VoronoiRegionSet


class StringDriver(object):
    def __init__(self, sim_manager):
        super(StringDriver, self).__init__()

        if sim_manager.work_manager.mode != sim_manager.work_manager.MODE_MASTER:
            return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system

        # Parameters from config file
        self.windowsize = west.rc.config.get_int('stringmethod.windowsize', 10)
        self.update_interval = west.rc.config.get_int('stringmethod.update_interval', 10)
        self.initial_update = west.rc.config.get_int('stringmethod.initial_update', 20)
        self.priority = west.rc.config.get_int('stringmethod.priority', 0)

        self.write_avg_pos = west.rc.config.get_bool('stringmethod.write_avgpos', True)
        self.do_update = west.rc.config.get_bool('stringmethod.do_update', True)
        self.init_from_data = west.rc.config.get_bool('stringmethod.init_from_data', True)

        self.dfunc = self.get_dfunc_method()

        # Load method to calculate average position in a bin
        # If the method is defined in an external module, correctly bind it
        ap = self.get_avgpos_method()
        if hasattr(ap,'im_class'):
            self.get_avgpos = ap
        else:
            self.get_avgpos = types.MethodType(ap,self)

        # Get initial set of string centers
        centers = self.get_initial_centers()

        try:
            sm_params = self.system.sm_params
        except AttributeError as e:
            log.error('String Driver Error: system does not define sm_params. \
                        This is required and should be added to the system definition; {}'.format(e))
            raise

        # Initialize the string
        str_method = self.get_string_method()

        try:
            self.strings = str_method(centers, **sm_params)
        except (TypeError, AssertionError) as e:
            log.error('String Driver Error: Failed during initialization of string method: {}'.format(e))
            raise

        # Update the Region Set
        self.update_regionset()

        # Register callback
        sim_manager.register_callback(sim_manager.prepare_new_segments, self.prepare_new_segments, self.priority)

        west.rc.pstatus('-westext.stringmethod -----------------\n')
        west.rc.pstatus('windowsize: {}\n'.format(self.windowsize))
        west.rc.pstatus('update interval: {}\n'.format(self.update_interval))
        west.rc.pstatus('initial update: {}\n'.format(self.initial_update))
        west.rc.pstatus('priority: {}\n'.format(self.priority))
        west.rc.pstatus('write average positions: {}\n'.format(self.write_avg_pos))
        west.rc.pstatus('do update: {}\n'.format(self.do_update))
        west.rc.pstatus('initialize from WE data: {}\n'.format(self.init_from_data))
        west.rc.pstatus('----------------------------------------\n') 
        west.rc.pflush()

    def dfunc(self):
        raise NotImplementedError

    def get_avgpos(self,n_iter):
        raise NotImplementedError

    def get_dfunc_method(self):
        methodname = west.rc.config.require('stringmethod.dfunc_method')
        pathinfo = west.rc.config.get_pathlist('stringmethod.module_path', default=None)
        dfunc_method = extloader.get_object(methodname, pathinfo)

        log.debug('loaded stringmethod dfunc method {!r}'.format(dfunc_method))

        return dfunc_method

    def get_avgpos_method(self):
        methodname = west.rc.config.require('stringmethod.avgpos_method')
        if methodname.lower() == 'cartesian':
            avgpos_method = self.avgpos_cartesian
        else:
            pathinfo = west.rc.config.get_pathlist('stringmethod.module_path', default=None)
            avgpos_method = extloader.get_object(methodname, pathinfo)

        log.debug('loaded stringmethod avgpos method {!r}'.format(avgpos_method))

        return avgpos_method

    def get_string_method(self):
        methodname = west.rc.config.require('stringmethod.string_method')
        if methodname.lower() == 'default':
            str_method = DefaultStringMethod
        else:
            pathinfo = west.rc.config.get_pathlist('stringmethod.module_path', default=None)
            str_method = extloader.get_object(methodname, pathinfo)

        assert issubclass(str_method, WESTStringMethod)
        log.debug('loaded stringmethod string method {!r}'.format(str_method))

        return str_method

    def get_initial_centers(self):
        self.data_manager.open_backing()
        
        with self.data_manager.lock:
            n_iter = max(self.data_manager.current_iteration - 1, 1)
            iter_group = self.data_manager.get_iter_group(n_iter)

            # First attempt to initialize string from data rather than system
            centers = None
            if self.init_from_data and 'stringmethod' in iter_group:
                log.info('Attempting to initialize stringmethod from data')

                try:
                    centers = iter_group['stringmethod']['centers'][...]
                except:
                    log.warning('Initializing string centers from data failed; Using definition in system instead.')
                    centers = self.system.curr_region_set.centers
            else:
                log.info('Initializing string centers from system definition')
                centers = self.system.curr_region_set.centers

        self.data_manager.close_backing()

        return centers

    def update_regionset(self):
        '''Update the Region Set using the current string'''

        west.rc.pstatus('westext.stringmethod: Updating bin definitions\n')
        west.rc.pflush()

        bins = self.system.curr_region_set.get_all_bins()
        target_counts = vgetattr('target_count', bins,np.int)

        try:
            self.system.curr_region_set = VoronoiRegionSet(self.dfunc,self.strings.centers)
            bins = self.system.curr_region_set.get_all_bins()
            for bi,bin in enumerate(bins):
                bin.target_count = target_counts[bi]
        except (ValueError,TypeError) as e:
            log.error('StringDriver Error: Failed updating region set: {}'.format(e))
            raise

    def avgpos_cartesian(self, n_iter):
        '''Get average position of replicas in each bin as of n_iter for the
        the user selected update interval'''

        ncenters = self.system.curr_region_set.ncenters
        ndim = self.system.pcoord_ndim

        avg_pos = np.zeros((ncenters, ndim), dtype=self.system.pcoord_dtype)
        sum_bin_weight = np.zeros((ncenters,), dtype=self.system.pcoord_dtype)

        start_iter = max(n_iter - min(self.windowsize, n_iter), 2)
        stop_iter = n_iter + 1

        for n in xrange(start_iter, stop_iter):
            with self.data_manager.lock:
                iter_group = self.data_manager.get_iter_group(n)
                seg_index = iter_group['seg_index'][...]

                pcoords = iter_group['pcoord'][:,-1,:]  # Only read final point
                bin_indices = self.system.curr_region_set.map_to_all_indices(pcoords)
                weights = seg_index['weight']

                pcoord_w = pcoords * weights[:,np.newaxis]
                uniq_indices = np.unique(bin_indices)

                for indx in uniq_indices:
                    avg_pos[indx,:] += pcoord_w[bin_indices == indx].sum(axis=0)

                sum_bin_weight += np.bincount(bin_indices.astype(np.int),weights=weights,minlength=ncenters)

        # Some bins might have zero samples so exclude to avoid divide by zero
        occ_ind = np.nonzero(sum_bin_weight)
        avg_pos[occ_ind] /= sum_bin_weight[occ_ind][:,np.newaxis]

        return avg_pos,sum_bin_weight

    def prepare_new_segments(self):

        n_iter = self.sim_manager.n_iter

        with self.data_manager.lock:
            iter_group = self.data_manager.get_iter_group(n_iter)

            ncenters = self.system.curr_region_set.ncenters
            ndim = self.system.pcoord_ndim

            try:
                del iter_group['stringmethod']
            except KeyError:
                pass

            sm_iter_group = iter_group.create_group('stringmethod')
            centers_ds = sm_iter_group.create_dataset('centers',shape=(ncenters,ndim),dtype=self.system.pcoord_dtype)
            sm_global_group = self.data_manager.h5file.require_group('stringmethod')
            last_update = long(sm_global_group.attrs.get('last_update', 0))

        if n_iter - last_update < self.update_interval or n_iter < self.initial_update or not self.do_update:
            log.debug('Not updating string this iteration')
            with self.data_manager.lock:
                centers_ds[...] = self.strings.centers
            return
        else:
            log.debug('Updating string - n_iter: {}'.format(n_iter))

        west.rc.pstatus('-westext.stringmethod -----------------\n')
        west.rc.pstatus('westext.stringmethod: Calculating average position in string images\n')
        west.rc.pflush()

        with self.data_manager.lock:
            avg_pos,sum_bin_weight = self.get_avgpos(n_iter)

            if self.write_avg_pos:
                avg_pos_ds = sm_iter_group.create_dataset('avgpos',shape=(ncenters,ndim),dtype=np.float64)
                avg_pos_ds[...] = avg_pos

        west.rc.pstatus('westext.stringmethod: Updating string\n')
        west.rc.pflush()

        self.strings.update_string_centers(avg_pos,sum_bin_weight)

        west.rc.pstatus('westext.stringmethod: String lengths: {}\n'.format(self.strings.length))
        west.rc.pflush()

        with self.data_manager.lock:
            centers_ds[...] = self.strings.centers

        # Update the bin definitions
        self.update_regionset()

        sm_global_group.attrs['last_update'] = n_iter
