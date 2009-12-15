import os, sys, subprocess, re
import numpy
from optparse import OptionParser

from wemd.util.binfiles import UniformTimeData

parser = OptionParser(usage='mkdist.py TRAJ_DIR GROUPS',
                      description = 'stitch together a bunch of GROMACS trajectories')
parser.add_option('-n', '--index-file', dest='index_file',
                  help='GROMACS index (.ndx) file')
parser.add_option('-r', '--rc-file', dest='rc_file',
                  help='GROMACS run control (.tpr) file')
parser.add_option('-t', '--timestep',  dest='timestep', type='float',
                  help='MD timestep (default: read from first distance calc)')
parser.add_option('--text-output', dest='text_output', default='dist.txt',
                  help='filename for text output (default: dist.txt)')
parser.add_option('--g_dist-verbose', dest='gdist_verbose',
                  action = 'store_true',
                  help = 'do not suppress g_dist output')
(opts, args) = parser.parse_args()

if len(args) != 2:
    parser.print_usage(sys.stderr)
    sys.exit(2)
else:
    traj_dir = args[0]
    dist_groups = args[1]

if not opts.rc_file:
    sys.stderr.write('a GROMACS run control (.tpr) file must be provided\n')
    parser.print_usage(sys.stderr)
    sys.exit(2)
    
if not opts.index_file:
    sys.stderr.write('a GROMACS index file (.ndx) must be provided\n')
    parser.print_usage(sys.stderr)
    sys.exit(2)
    
reTrajName = re.compile(r'(.*?)(\d+)\.xtc$')

trajfiles = []
for fn in os.listdir(traj_dir):
    m = reTrajName.match(fn)
    if not m: continue

    groups = m.groups()
    
    fn_prefix = groups[0]
    ntraj = int(groups[1])
    trajfiles.append((ntraj, fn))
trajfiles = [i[1] for i in sorted(trajfiles)]

txt_out_fmt = '%-20.15g    %-12.8g    %-20.15g\n'
txt_out = open(opts.text_output, 'wt')

dt = opts.timestep
last_dist_txt = None
ti = -1

base_args = ['g_dist', '-noxvgr', '-s', opts.rc_file, '-n', opts.index_file]

for (i, trajfile) in enumerate(trajfiles):
    basename = trajfile[:-4]
    distfilename = 'dist-' + basename + '.xvg'
    trajfile_fullpath = os.path.join(traj_dir, trajfile)
    
    sys.stdout.write('g_traj %s -> %s\n' % (trajfile_fullpath, distfilename))
    sys.stdout.flush()
    if opts.gdist_verbose:
        child_output = None
    else:
        child_output = open('/dev/null', 'wb')
    proc = subprocess.Popen(base_args + ['-f', trajfile_fullpath, 
                                         '-o', distfilename],
                            stdin = subprocess.PIPE,
                            stdout = child_output,
                            stderr = subprocess.STDOUT)                  
    proc.stdin.write(dist_groups + '\n')
    rc = proc.wait()
    
    if rc != 0:
        sys.stderr.write('g_traj exited with code %r; aborting\n' % rc)
        sys.exit(1)

    dist_in = open(distfilename, 'rt')

    if dt is None:
        line1 = dist_in.readline()
        line2 = dist_in.readline()
        
        fields1 = line1.split()
        fields2 = line2.split()
        
        dt = float(fields2[0]) - float(fields1[0])
        sys.stdout.write('using dt = %g\n' % dt)
        dist_in.seek(0)
        
    sys.stdout.write('processing %s\n' % distfilename)

    if i > 0:
        line = dist_in.readline()
        fields = line.split()
        if fields[1] != last_dist_txt:
            sys.stderr.write('coordinates do not overlap correctly\n')
            sys.exit(1)

    for line in dist_in:
        ti += 1
        t = dt * ti
        fields = line.split()
        last_dist_txt = fields[1]
        dist = float(fields[1])

        txt_out.write(txt_out_fmt % (t, dist, weight))

    txt_out.flush()
    dist_in.close()

txt_out.close()
