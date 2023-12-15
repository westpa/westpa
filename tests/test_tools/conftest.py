import pytest
import os
import glob

from shutil import copyfile, copy
import tempfile

import westpa

REFERENCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'refs')

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


@pytest.fixture
def ref_3iter(request):
    """
    Fixture that prepares a simulation directory with a completed 3-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    test_dir = tempfile.mkdtemp()
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_3iter.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_cfg(request, tmpdir):
    """
    Fixture that prepares a simulation directory with a populated west.cfg file.
    """

    test_dir = str(tmpdir)
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.h5'), "west_init_ref.h5")

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME
    request.cls.ref_h5_filepath = 'west_init_ref.h5'

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_initialized(request, tmpdir):
    """
    Fixture that prepares a simulation directory with an initialized WESTPA system,
    west.h5, plus the config file west.cfg
    """

    test_dir = str(tmpdir)

    os.chdir(test_dir)
    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_init_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_50iter(request, tmpdir):
    """
    Fixture that prepares a simulation directory with a completed 50-iteration WESTPA,
    west.h5, plus the config file west.cfg
    """

    test_dir = str(tmpdir)

    os.chdir(test_dir)
    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.cfg'), CFG_FILENAME)

    request.cls.cfg_filepath = CFG_FILENAME
    request.cls.h5_filepath = H5_FILENAME

    os.environ['WEST_SIM_ROOT'] = test_dir
    westpa.rc = westpa.core._rc.WESTRC()

    request.addfinalizer(clear_state)


@pytest.fixture
def ref_multi(request, tmpdir):
    """
    Fixture that prepares a simulation directory for w_multi_west, including a master
    folder with sub folders 01, 02, 03 containing west_aux_ref.h5 renamed as west.h5.
    """

    test_dir = str(tmpdir)

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
def ref_multi_noaux(request, tmpdir):
    """
    Fixture that prepares a simulation directory for w_multi_west, including a master
    folder with sub folders 01, 02, 03 containing west_aux_ref.h5 renamed as west.h5.
    """

    test_dir = str(tmpdir)

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
def ref_idtype(request):
    """
    Fixture that prepares the west.h5 file and also links in the "correct" istate dtype array.
    """
    test_dir = tempfile.mkdtemp()
    os.chdir(test_dir)

    copy_ref(test_dir)

    copyfile(os.path.join(REFERENCE_PATH, 'west_ref.h5'), H5_FILENAME)
    copyfile(os.path.join(REFERENCE_PATH, 'ref_dtype.pickle'), 'ref_dtype.pickle')

    request.cls.h5_filepath = H5_FILENAME
    request.cls.correct_pkl = 'ref_dtype.pickle'

    os.environ['WEST_SIM_ROOT'] = test_dir

    request.addfinalizer(clear_state)


def clear_state():
    os.chdir(STARTING_PATH)
    del os.environ['WEST_SIM_ROOT']
    westpa.rc = westpa.core._rc.WESTRC()
