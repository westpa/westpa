import pytest
import os

from shutil import copyfile

from westpa import rc

REFERENCE_PATH = os.path.join(os.path.dirname(__file__), '../refs')


H5_FILENAME = 'west.h5'
CFG_FILENAME = 'west.cfg'
refs_cfg_file = 'west_ref.cfg'
refs_h5_file = 'west_ref.h5'

# TODO: Is putting this here, outside of any function bad practice?
#   Need it here so that clear_state doesn't take an argument...
STARTING_PATH = os.getcwd()

@pytest.fixture
def ref_3iter(request):
    """
    Fixture that prepares a simulation directory with a completed 3-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    os.chdir(REFERENCE_PATH)

    copyfile('west_3iter.h5', H5_FILENAME)
    copyfile(REF_CFG_FILENAME, CFG_FILENAME)

    request.cls.cfg_filepath = os.path.join(REFERENCE_PATH, CFG_FILENAME)
    request.cls.h5_filepath = os.path.join(REFERENCE_PATH, H5_FILENAME)

    os.environ['WEST_SIM_ROOT'] = REFERENCE_PATH

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_cfg(request):
    """
    Fixture that prepares a simulation directory with a populated west.cfg file.
    """

    os.chdir(REFERENCE_PATH)

    copyfile(REF_CFG_FILENAME, CFG_FILENAME)

    request.cls.cfg_filepath = os.path.join(REFERENCE_PATH, CFG_FILENAME)
    request.cls.h5_filepath = os.path.join(REFERENCE_PATH, H5_FILENAME)
    request.cls.ref_h5_filepath = os.path.join(REFERENCE_PATH, 'west_init_ref.h5')

    os.environ['WEST_SIM_ROOT'] = REFERENCE_PATH

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_initialized(request):
    """
    Fixture that prepares a simulation directory with an initialized WESTPA system,
    west.h5, plus the config file west.cfg
    """

    os.chdir(REFERENCE_PATH)

    copyfile('west_init_ref.h5', H5_FILENAME)
    copyfile(REF_CFG_FILENAME, CFG_FILENAME)

    request.cls.cfg_filepath = os.path.join(REFERENCE_PATH, CFG_FILENAME)
    request.cls.h5_filepath = os.path.join(REFERENCE_PATH, H5_FILENAME)
    request.cls.REFERENCE_PATH = REFERENCE_PATH

    os.environ['WEST_SIM_ROOT'] = REFERENCE_PATH

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_50iter(request):
    """
    Fixture that prepares a simulation directory with a completed 50-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    os.chdir(REFERENCE_PATH)

    copyfile(refs_h5_file, H5_FILENAME)
    copyfile(refs_cfg_file, CFG_FILENAME)

    request.cls.cfg_filepath = os.path.join(REFERENCE_PATH, CFG_FILENAME)
    request.cls.h5_filepath = os.path.join(REFERENCE_PATH, H5_FILENAME)

    os.environ['WEST_SIM_ROOT'] = REFERENCE_PATH

    request.addfinalizer(clear_state)


def clear_state():

    os.remove(CFG_FILENAME)
    os.remove(H5_FILENAME)

    os.chdir(STARTING_PATH)

    os.environ['WEST_SIM_ROOT'] = ''

    rc._sim_manager = None
    rc._system = None
    rc._data_manager = None
    rc._we_driver = None
    rc._propagator = None
