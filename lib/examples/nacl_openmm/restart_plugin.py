import logging
log = logging.getLogger(__name__)


class RestartDriver(object):
    def __init__(self, sim_manager, plugin_config):
        super(RestartDriver, self).__init__()

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.priority = plugin_config.get('priority', 0)

        # Register callback
        sim_manager.register_callback(sim_manager.pre_propagation, self.pre_propagation, self.priority)

    def pre_propagation(self):

        segments = list(self.sim_manager.incomplete_segments.values())
        n_iter = self.sim_manager.n_iter

        if n_iter == 1:
            return

        parent_iter_group = self.data_manager.get_iter_group(n_iter - 1)

        # Get parent ids for segments
        parent_ids = [seg.parent_id for seg in segments]

        # Get a list of unique parent ids and collect restart data for each
        unique_parent_ids = set(parent_ids)
        restart_data = {segid: {} for segid in unique_parent_ids}

        for dsname in ['coord', 'veloc']:
            try:
                dsinfo = self.data_manager.dataset_options[dsname]
            except KeyError:
                raise KeyError('Data set {} not found'.format(dsname))

            ds = parent_iter_group[dsinfo['h5path']]

            for seg_id in unique_parent_ids:
                restart_data[seg_id][dsname] = ds[seg_id][-1,...]

        for segment in segments:
            segment.data['restart_coord'] = restart_data[segment.parent_id]['coord']
            segment.data['restart_veloc'] = restart_data[segment.parent_id]['veloc']
