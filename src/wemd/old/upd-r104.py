import sys, os, shutil, re
import cPickle as pickle
import wemd

sim_config_file = 'sim_config.pkl'
sim_config = {}

runtime_config = wemd.rc.read_config()
if 'data.sim_config' not in runtime_config:
    sys.stdout.write('adding sim_config entry to %s\n' % runtime_config['__file__'])
    backup_rcfile = '%s.orig' % runtime_config['__file__']
    shutil.copy2(runtime_config['__file__'], backup_rcfile)
    rcfile_in = open(backup_rcfile, 'rt')
    rcfile_out = open(runtime_config['__file__'], 'wt')
    
    line = rcfile_in.readline()
    while line and line.strip() != '[data]':
        rcfile_out.write(line)
        line = rcfile_in.readline()
    else:
        rcfile_out.write(line)
        rcfile_out.write('sim_config = %s\n' % sim_config_file)
    line = rcfile_in.readline()
    while line:
        rcfile_out.write(line)
        line = rcfile_in.readline()
    rcfile_out.close()
    rcfile_in.close()
    
    runtime_config = wemd.rc.read_config()
    
sim_config['backend.driver'] = runtime_config['backend.driver']
sim_config['bins.type'] = 'fixed' # if this isn't right, "rebin" will fix it

sys.stdout.write('creating %s\n' % sim_config_file)
pickle.dump(sim_config, open(sim_config_file, 'wb'), pickle.HIGHEST_PROTOCOL)
sys.stdout.write('run wemdctl.py rebin now\n')

