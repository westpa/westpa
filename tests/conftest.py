import pytest
import os
import glob
from shutil import copyfile, copy

import numpy as np
from scipy.io import netcdf_file

import westpa
from westpa.core.h5io import WESTIterationFile

REFERENCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'refs')

H5_FILENAME = 'west.h5'
CFG_FILENAME = 'west.cfg'
refs_cfg_file = 'west_ref.cfg'
refs_h5_file = 'west_ref.h5'

# TODO: Is putting this here, outside of any function bad practice?
#   Need it here so that clear_state doesn't take an argument...
STARTING_PATH = os.getcwd()


def copy_ref(dest_dir):
    for filename in glob.glob(os.path.join(REFERENCE_PATH, '*.*')):
        copy(filename, dest_dir)


def clear_state():
    os.chdir(STARTING_PATH)
    if 'WEST_SIM_ROOT' in os.environ:
        del os.environ['WEST_SIM_ROOT']
    westpa.rc = westpa.core._rc.WESTRC()


@pytest.fixture
def ref_3iter(request, tmp_path):
    """Fixture that prepares a simulation directory with a completed 3-iteration WESTPA,
    west.h5, plus the config file west.cfg"""

    test_dir = str(tmp_path)
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_3iter.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_cfg(request, tmp_path):
    """Fixture that prepares a simulation directory with a populated west.cfg file."""

    test_dir = str(tmp_path)
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)

    if np.__version__ <= '2.0.0':
        copyfile(os.path.join(REFERENCE_PATH, 'west_init_numpy1.h5'), "west_init_ref.h5")
    else:
        copyfile(os.path.join(REFERENCE_PATH, 'west_init_numpy2.h5'), "west_init_ref.h5")

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME
    request.cls.ref_h5_filepath = 'west_init_ref.h5'

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_initialized(request, tmp_path):
    """
    Fixture that prepares a simulation directory with an initialized WESTPA system,
    west.h5, plus the config file west.cfg
    """

    test_dir = str(tmp_path)

    os.chdir(test_dir)
    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_init_numpy2.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_50iter(request, tmp_path):
    """
    Fixture that prepares a simulation directory with a completed 50-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    test_dir = str(tmp_path)

    os.chdir(test_dir)
    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.cfg'), CFG_FILENAME)

    analysis_path = f'{test_dir}/ANALYSIS/TEST'
    os.makedirs(analysis_path, exist_ok=True)
    copyfile(os.path.join(REFERENCE_PATH, 'assign_ref.h5'), f'{analysis_path}/assign.h5')
    copyfile(os.path.join(REFERENCE_PATH, 'direct_ref.h5'), f'{analysis_path}/direct.h5')

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_multi(request, tmp_path):
    """
    Fixture that prepares a simulation directory for w_multi_west, including a master
    folder with sub folders 01, 02, 03 containing west_aux_ref.h5 renamed as west.h5.
    """

    test_dir = str(tmp_path)

    os.chdir(test_dir)
    copy_ref(test_dir)
    os.mkdir('01')
    os.mkdir('02')
    os.mkdir('03')

    copyfile(os.path.join(REFERENCE_PATH, 'west_aux_ref.h5'), "01/west.h5")
    copyfile(os.path.join(REFERENCE_PATH, 'west_aux_ref.h5'), "02/west.h5")
    copyfile(os.path.join(REFERENCE_PATH, 'west_aux_ref.h5'), "03/west.h5")

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_multi_noaux(request, tmp_path):
    """Fixture that prepares a simulation directory for w_multi_west, including a master
    folder with sub folders 01, 02, 03 containing west_aux_ref.h5 renamed as west.h5."""

    test_dir = str(tmp_path)

    os.chdir(test_dir)
    copy_ref(test_dir)
    os.mkdir('01')
    os.mkdir('02')
    os.mkdir('03')

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), "01/west.h5")
    copyfile(os.path.join(REFERENCE_PATH, 'west_aux_ref.h5'), "02/west.h5")
    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), "03/west.h5")

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_idtype(request, tmp_path):
    """Fixture that prepares the west.h5 file and also links in the "correct" istate dtype array."""

    test_dir = str(tmp_path)
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'ref_dtype.pickle'), 'ref_dtype.pickle')

    request.cls.h5_filepath = H5_FILENAME
    request.cls.correct_pkl = 'ref_dtype.pickle'

    os.environ['WEST_SIM_ROOT'] = test_dir

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_executable(request, tmp_path):
    """Fixture that prepares a simulation directory with a populated west_executable.cfg file."""

    test_dir = str(tmp_path)
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_executable.cfg'), CFG_FILENAME)

    # Create a "temp" file
    with open(CFG_FILENAME, 'r') as f, open('west_implicit.cfg', 'w') as g:
        for line in f.readlines()[:25]:
            g.write(line)

    request.cls.cfg_filepath = CFG_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture(scope='function')
def west_iteration_file(request, tmp_path):
    os.chdir(tmp_path)
    request.cls.h5_iter_file_path = tmp_path / 'WESTITERFILE.h5'

    request.cls.rng = rng = np.random.default_rng()
    request.cls.dummy_data = {'iterh5/trajectory': rng.uniform(low=-3, high=3, size=(4, 5, 3))}

    # Initialize and close the file
    WESTIterationFile(request.cls.h5_iter_file_path, mode='w').close()

    request.addfinalizer(clear_state)


@pytest.fixture
def traj_setup(request, tmp_path):
    """Fixture for testing the trajectory reading capabilities of the HDF5 Framework"""

    test_dir = str(tmp_path)

    os.chdir(test_dir)

    traj_file_path = os.path.join(REFERENCE_PATH, 'ntl9.nc')
    top_file_path = os.path.join(REFERENCE_PATH, 'ntl9_reference.pdb')

    copyfile((traj_file_path), 'ntl9.nc')
    copyfile(top_file_path, 'ntl9_reference.pdb')

    request.cls.current_path = tmp_path

    with netcdf_file(traj_file_path) as rootgrp:
        request.cls.ref_coords = rootgrp.variables['coordinates'][()].copy() / 10
        # Not all simulations have periodic boundaries
        # request.cls.ref_lengths = rootgrp.variables['cell_lengths'][()]
        # request.cls.ref_angles = rootgrp.variables['cell_angles'][()]
        request.cls.ref_time = rootgrp.variables['time'][()].copy()


@pytest.fixture
def ref_mab(request, tmp_path):
    """Fixture that prepares an rc/sim_manager/WESTSystem from west_mab.cfg"""

    test_dir = str(tmp_path)

    os.chdir(test_dir)
    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_mab.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()
    westpa.rc.read_config(filename='west.cfg')

    request.cls.tmpdir = test_dir


@pytest.fixture
def nacl_restart_files(request, tmp_path):
    request.cls.test_dir = tmp_path
    request.cls.return_dir = tmp_path / 'restart_return'
    request.cls.write_dir = tmp_path / 'restart_write'

    request.cls.nacl_restart_files = ['nacl.prmtop', 'nacl.ncrst']

    os.chdir(tmp_path)
    os.mkdir(request.cls.return_dir)
    os.mkdir(request.cls.write_dir)

    for file in request.cls.nacl_restart_files:
        copyfile(os.path.join(REFERENCE_PATH, file), request.cls.return_dir / file)
