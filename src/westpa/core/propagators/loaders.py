import logging
import numpy as np
import pickle
import tarfile
from io import BytesIO

import westpa
from westpa.core.trajectory import load_mdtraj, load_netcdf, load_mdanalysis
from westpa.core.h5io import safe_extract

log = logging.getLogger(__name__)


def pcoord_loader(fieldname, pcoord_return_filename, destobj, single_point):
    """Read progress coordinate data into the ``pcoord`` field on ``destobj``.
    An exception will be raised if the data is malformed.  If ``single_point`` is true,
    then only one (N-dimensional) point will be read, otherwise system.pcoord_len points
    will be read.
    """

    system = westpa.rc.get_system_driver()

    assert fieldname == 'pcoord'

    pcoord = np.loadtxt(pcoord_return_filename, dtype=system.pcoord_dtype)

    if single_point:
        expected_shape = (system.pcoord_ndim,)
        if pcoord.ndim == 0:
            pcoord.shape = (1,)
    else:
        expected_shape = (system.pcoord_len, system.pcoord_ndim)
        if pcoord.ndim < 2:
            pcoord.shape = expected_shape
    if pcoord.shape != expected_shape:
        raise ValueError(
            'progress coordinate data has incorrect shape {!r} [expected {!r}] Check pcoord.err or seg_logs for more '
            'information.'.format(pcoord.shape, expected_shape)
        )
    destobj.pcoord = pcoord


def aux_data_loader(fieldname, data_filename, segment, single_point):
    data = np.loadtxt(data_filename)
    segment.data[fieldname] = data
    if data.nbytes == 0:
        raise ValueError('could not read any data for {}'.format(fieldname))


def numpy_data_loader(fieldname, coord_file, segment, single_point):
    log.debug('using numpy_data_loader')
    data = np.load(coord_file, allow_pickle=True)
    segment.data[fieldname] = data
    if data.nbytes == 0:
        raise ValueError('could not read any data for {}'.format(fieldname))


def pickle_data_loader(fieldname, coord_file, segment, single_point):
    log.debug('using pickle_data_loader')
    with open(coord_file, 'rb') as fo:
        data = pickle.load(fo)
    segment.data[fieldname] = data
    if data.nbytes == 0:
        raise ValueError('could not read any data for {}'.format(fieldname))


def mdtraj_trajectory_loader(fieldname, coord_folder, segment, single_point):
    '''Load data from the trajectory return using MDTraj. ``coord_folder`` should be the path to a folder
    containing trajectory files. ``segment`` is the ``Segment`` object that the data is associated with.
    Please see ``load_mdtraj`` for more details. ``single_point`` is not used by this loader.'''
    try:
        data = load_mdtraj(coord_folder)
        segment.data['iterh5/trajectory'] = data
    except Exception as e:
        log.warning('could not read any {} data for HDF5 Framework: {}'.format(fieldname, str(e)))


def netcdf_trajectory_loader(fieldname, coord_folder, segment, single_point):
    '''Load Amber .nc data from the trajectory return. ``coord_folder`` should be the path to a folder
    containing trajectory files. ``segment`` is the ``Segment`` object that the data is associated with.
    Please see ``load_netcdf`` for more details. ``single_point`` is not used by this loader.'''
    try:
        data = load_netcdf(coord_folder)
        segment.data['iterh5/trajectory'] = data
    except Exception as e:
        log.warning('Falling back to default loader for {}: {}'.format(fieldname, str(e)))
        mdtraj_trajectory_loader(fieldname, coord_folder, segment, single_point)


def mdanalysis_trajectory_loader(fieldname, coord_folder, segment, single_point):
    '''Load data from the trajectory return. ``coord_folder`` should be the path to a folder
    containing trajectory files. ``segment`` is the ``Segment`` object that the data is associated with.
    Please see ``load_mdanalysis`` for more details. ``single_point`` is not used by this loader.'''
    try:
        data = load_mdanalysis(coord_folder)
        segment.data['iterh5/trajectory'] = data
    except Exception as e:
        log.warning('Falling back to default loader for {}: {}'.format(fieldname, str(e)))
        mdtraj_trajectory_loader(fieldname, coord_folder, segment, single_point)


def restart_loader(fieldname, restart_folder, segment, single_point):
    '''Load data from the restart return. The loader will tar all files in ``restart_folder``
    and store it in the per-iteration HDF5 file. ``segment`` is the ``Segment`` object that
    the data is associated with. ``single_point`` is not used by this loader.'''
    try:
        with BytesIO() as d:
            with tarfile.open(mode='w:gz', fileobj=d) as t:
                t.add(restart_folder, arcname='.')

            segment.data['iterh5/restart'] = d.getvalue() + b'\x01'  # add tail protection
    except Exception as e:
        log.warning('could not read any {} data for HDF5 Framework: {}'.format(fieldname, str(e)))
    finally:
        d.close()


def restart_writer(path, segment):
    '''Prepare the necessary files from the per-iteration HDF5 file to run ``segment``.'''
    try:
        restart = segment.data.pop('iterh5/restart', None)
        # Making an exception for start states in iteration 1
        if restart is None:
            raise ValueError('restart data is not present')

        with BytesIO(restart[:-1]) as d:  # remove tail protection
            with tarfile.open(fileobj=d, mode='r:gz') as t:
                safe_extract(t, path=path)
    except ValueError as e:
        log.warning('could not write HDF5 Framework restart data for {}: {}'.format(str(segment), str(e)))
        if segment.n_iter == 1:
            log.warning(
                'In iteration 1. Assuming this is a start state and proceeding to skip reading restart from per-iteration HDF5 file for {}'.format(
                    str(segment)
                )
            )
    except Exception as e:
        log.warning('could not write HDF5 Framework restart data for {}: {}'.format(str(segment), str(e)))


def seglog_loader(fieldname, log_file, segment, single_point):
    '''Load data from the log return. The loader will tar all files in ``log_file``
    and store it in the per-iteration HDF5 file. ``segment`` is the ``Segment`` object that
    the data is associated with. ``single_point`` is not used by this loader.'''
    try:
        with BytesIO() as d:
            with tarfile.open(mode='w:gz', fileobj=d) as t:
                t.add(log_file, arcname='seg.log')

            segment.data['iterh5/log'] = d.getvalue() + b'\x01'  # add tail protection
    except Exception as e:
        log.warning('could not read any data for {}: {}'.format(fieldname, str(e)))


def seglog_writer(path, segment):
    '''Untar the log file from segment.'''
    try:
        seglog = segment.data.pop('iterh5/log', None)

        # Raising an exception for start states in iteration 1
        if seglog is None:
            raise ValueError('Log file is not present for segment: {}'.format(str(segment)))

        with BytesIO(seglog[:-1]) as d:  # remove tail protection
            with tarfile.open(fileobj=d, mode='r:gz') as t:
                safe_extract(t, path=path)
    except Exception as e:
        log.warning('could not extract HDF5 Framework log file for {}: {}'.format(str(segment), str(e)))


# Dictionary with all the possible aux dataset loaders
data_loaders = {
    'default': aux_data_loader,
    'auxdata_loader': aux_data_loader,
    'aux_data_loader': aux_data_loader,
    'numpy_loader': numpy_data_loader,
    'npy_loader': numpy_data_loader,
    'numpy_data_loader': numpy_data_loader,
    'npy_data_loader': numpy_data_loader,
    'pickle_loader': pickle_data_loader,
    'pickle_data_loader': pickle_data_loader,
}


# Dictionary with all the possible trajectory loaders
trajectory_loaders = {
    'default': mdtraj_trajectory_loader,
    'trajectory_loader': mdtraj_trajectory_loader,
    'MDTraj_trajectory_loader': mdtraj_trajectory_loader,
    'mdtraj_trajectory_loader': mdtraj_trajectory_loader,
    'amber_trajectory_loader': netcdf_trajectory_loader,
    'netcdf_trajectory_loader': netcdf_trajectory_loader,
    'mda_trajectory_loader': mdanalysis_trajectory_loader,
    'MDAnalysis_trajectory_loader': mdanalysis_trajectory_loader,
    'mdanalysis_trajectory_loader': mdanalysis_trajectory_loader,
}
