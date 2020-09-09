import pytest
import os

from shutil import copyfile

from westpa.core.yamlcfg import YAMLConfig
from westpa import rc


@pytest.fixture
def initialized_west_ref(request):

    """
    Setup code to copy the reference .cfg file
    """

    starting_path = os.getcwd()

    odld_path = os.path.dirname(__file__) + '/ref'

    os.chdir(odld_path)

    copyfile('west_init_ref.cfg', 'west.cfg')

    request.cls.cfg_path = os.path.join(odld_path, 'west.cfg')
    request.cls.odld_path = odld_path
    os.environ['WEST_SIM_ROOT'] = odld_path

    def fin():

        os.remove('west.cfg')
        os.chdir(starting_path)
        os.environ['WEST_SIM_ROOT'] = ''

    request.addfinalizer(fin)


@pytest.fixture
def w_run_fixture(request):

    starting_path = os.getcwd()

    odld_path = os.path.dirname(__file__) + '/ref'

    rc.config = YAMLConfig()

    os.chdir(odld_path)
    copyfile('west_ref.h5', 'west.h5')
    copyfile('west_init_ref.cfg', 'west.cfg')

    os.environ['WEST_SIM_ROOT'] = odld_path
    request.cls.cfg_filepath = os.path.join(odld_path, 'west.cfg')
    request.cls.h5_filepath = os.path.join(odld_path, 'west.h5')
    request.cls.odld_path = odld_path

    def fin():

        os.remove('west.cfg')
        os.remove('west.h5')
        os.chdir(starting_path)
        os.environ['WEST_SIM_ROOT'] = ''

    request.addfinalizer(fin)
