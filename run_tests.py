import os
import sys
import subprocess

from tempfile import gettempdir

from argparse import ArgumentParser


TMP_DIR = "TMPDIR"
PYTHON_PATH = "PYTHONPATH"
# MPI tests barf if setting a different temp directory on OSX
PREFERRED_TEMP_PATH = "/tmp"

# Paths that need to be added to PYTHONPATH
SRC_PATH = "src"
WORK_MANAGER_PATH = "lib/wwmgr"
WEST_TOOLS_PATH = "lib/west_tools/"

TEST_DIRECTORIES = [
    "src/west/",
    "lib/west_tools/",
    "lib/wwmgr/", 
]


def test_directory(directory: str, timeout: int = 60) -> bool:
    """
    Run a the tests for a specific directory using nosetest

    :param string directory: Directory to run nosetest in
    :param integer timeout: How long before cancelling the tests for a directory
    :rtype: bool
    """
    cur_dir = os.getcwd()
    working_dir = os.path.join(cur_dir, directory)
    if not os.path.isdir(working_dir):
        print("Unable to find directory '{}'".format(working_dir))
        return
    print("Running tests for directory '{}'".format(working_dir))
    env = os.environ 
    work_managers_path = os.path.join(cur_dir, WORK_MANAGER_PATH)
    src_path = os.path.join(cur_dir, SRC_PATH)
    west_tools_path = os.path.join(cur_dir, WEST_TOOLS_PATH)
    env[PYTHON_PATH] = os.pathsep.join([work_managers_path, src_path, west_tools_path])
    if os.path.isdir(PREFERRED_TEMP_PATH):
        env[TMP_DIR] = PREFERRED_TEMP_PATH
    else: 
        env[TMP_DIR] = gettempdir()
    try:
        proc = subprocess.run(["nosetests"], cwd=working_dir, check=False, env=env, timeout=timeout)
        return proc.returncode == 0
    except Exception:
        return False

def main():
    parser = ArgumentParser(description="Run Westpa tests") 
    parser.add_argument("--directory", default=None)
    args = parser.parse_args()
    success = True
    if args.directory:
        success &= test_directory(args.directory)
    else:
        for directory in TEST_DIRECTORIES:
            success &= test_directory(directory)
    if not success:
        print("Not all tests passed")
        sys.exit(1)



if __name__ == "__main__":
    main()
