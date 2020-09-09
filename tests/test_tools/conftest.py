import pytest
import os

from shutil import copyfile


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
