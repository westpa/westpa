import os

template = '''\

# source this file to set the environment variables in your shell.

export WEST_ROOT={WESTROOT}
export WEST_BIN={WESTBIN}
export WEST_PYTHON={WESTPYTHON}

export PATH=$WEST_BIN:$PATH
\n'''

with open('westpa.sh', 'w') as f:
    f.write(template.format(WESTROOT=os.environ.get('PWD'),
                            WESTBIN=os.path.join(os.environ.get('PWD'), 'bin'),
                            WESTPYTHON=os.environ.get('WEST_PYTHON','')))
