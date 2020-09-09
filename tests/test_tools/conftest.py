import pytest
import os

from shutil import copyfile

from westpa import rc


@pytest.fixture
def ref_3iter(request):
    """
    Fixture that prepares a simulation directory with a completed 3-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    starting_path = os.getcwd()
    odld_path = os.path.dirname(__file__) + '/ref'

    os.chdir(odld_path)
    copyfile('west_3iter.h5', 'west.h5')
    copyfile('west_init_ref.cfg', 'west.cfg')

    request.cls.h5_filepath = os.path.join(odld_path, 'west.h5')
    request.cls.cfg_filepath = os.path.join(odld_path, 'west.cfg')
    os.environ['WEST_SIM_ROOT'] = odld_path

    def fin():

        os.remove('west.cfg')
        os.remove('west.h5')
        os.chdir(starting_path)

        rc._sim_manager = None
        rc._system = None
        rc._data_manager = None
        rc._we_driver = None
        rc._propagator = None

    request.addfinalizer(fin)


@pytest.fixture
def ref_cfg(request):
    """
    Fixture that prepares a simulation directory with a populated west.cfg file.
    """

    starting_path = os.getcwd()

    odld_path = os.path.dirname(__file__) + '/ref'

    os.chdir(odld_path)

    copyfile('west_init_ref.cfg', 'west.cfg')

    request.cls.cfg_filepath = os.path.join(odld_path, 'west.cfg')
    request.cls.h5_filepath = os.path.join(odld_path, 'west.h5')
    request.cls.ref_h5_filepath = os.path.join(odld_path, 'west_ref.h5')

    os.environ['WEST_SIM_ROOT'] = odld_path

    def fin():

        os.remove('west.cfg')
        os.remove('west.h5')
        os.chdir(starting_path)
        os.environ['WEST_SIM_ROOT'] = ''

        rc._sim_manager = None
        rc._system = None
        rc._data_manager = None
        rc._we_driver = None
        rc._propagator = None

    request.addfinalizer(fin)


@pytest.fixture
def ref_initialized(request):
    """
    Fixture that prepares a simulation directory with an initialized WESTPA system,
    west.h5, plus the config file west.cfg
    """

    starting_path = os.getcwd()

    odld_path = os.path.dirname(__file__) + '/ref'

    os.chdir(odld_path)
    copyfile('west_ref.h5', 'west.h5')
    copyfile('west_init_ref.cfg', 'west.cfg')

    os.environ['WEST_SIM_ROOT'] = odld_path
    request.cls.cfg_filepath = os.path.join(odld_path, 'west.cfg')
    request.cls.h5_filepath = os.path.join(odld_path, 'west.h5')
    request.cls.odld_path = odld_path

    def fin():

        print("Cleaning up w_run")
        os.remove('west.cfg')
        os.remove('west.h5')
        os.chdir(starting_path)
        os.environ['WEST_SIM_ROOT'] = ''

        rc._sim_manager = None
        rc._system = None
        rc._data_manager = None
        rc._we_driver = None
        rc._propagator = None

    request.addfinalizer(fin)
