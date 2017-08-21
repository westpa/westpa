import os

template_sh = '''\

# source this file to set the environment variables in your shell.

export WEST_ROOT={WESTROOT}
export WEST_BIN={WESTBIN}
export WEST_PYTHON={WESTPYTHON}

export PATH=$WEST_BIN:$PATH
\n'''

template_csh = '''\

# source this file to set the environment variables in your shell.

setenv WEST_ROOT "{WESTROOT}"
setenv WEST_BIN "{WESTBIN}"
setenv WEST_PYTHON "{WESTPYTHON}"

\n'''

with open('westpa.sh', 'w') as f:
    f.write(template_sh.format(WESTROOT=os.environ.get('PWD'),
                            WESTBIN=os.path.join(os.environ.get('PWD'), 'bin'),
                            WESTPYTHON=os.environ.get('WEST_PYTHON','')))

with open('westpa.csh', 'w') as f:
    f.write(template_csh.format(WESTROOT=os.environ.get('PWD'),
                            WESTBIN=os.path.join(os.environ.get('PWD'), 'bin'),
                            WESTPYTHON=os.environ.get('WEST_PYTHON','')))
    string = 'setenv PATH "' + os.path.join(os.environ.get('PWD'), 'bin') +\
            ':$PATH"\n'

    f.write(string)
