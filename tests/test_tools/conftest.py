import pytest
import os

from shutil import copyfile

from westpa import rc

refs_dir = os.path.join(os.path.dirname(__file__), 'refs')


H5_FILENAME = 'west.h5'
CFG_FILENAME = 'west.cfg'
refs_cfg_file = 'west_ref.cfg'
refs_h5_file = 'west_ref.h5'

# TODO: Is putting this here, outside of any function bad practice?
#   Need it here so that clear_state doesn't take an argument...
STARTING_PATH = os.getcwd()


@pytest.fixture
def ref_50iter(request):
    """
    Fixture that prepares a simulation directory with a completed 50-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    os.chdir(refs_dir)

    copyfile(refs_h5_file, H5_FILENAME)
    copyfile(refs_cfg_file, CFG_FILENAME)

    request.cls.cfg_filepath = os.path.join(refs_dir, CFG_FILENAME)
    request.cls.h5_filepath = os.path.join(refs_dir, H5_FILENAME)

    os.environ['WEST_SIM_ROOT'] = refs_dir

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
