import os, sys, subprocess, re
from optparse import OptionParser

parser = OptionParser(usage='mkdist.py -n INDEX_FILE -r RC_FILE GROUP_1 GROUP_2 TRAJ_FILES...',
                      description = 'stitch together a bunch of GROMACS trajectories')
parser.add_option('-d', '--distance-program', dest='distance_program',
                  help='use DISTANCE_PROGRAM to calculate distance (default: '
                      +'g_dist)',
                  default='g_dist')
parser.add_option('-n', '--index-file', dest='index_file',
                  help='GROMACS index (.ndx) file')
parser.add_option('-r', '--rc-file', dest='rc_file',
                  help='GROMACS run control (.tpr) file')
parser.add_option('-t', '--timestep',  dest='timestep', type='float',
                  help='MD timestep (default: read from first distance calc)')
parser.add_option('-o', '--output', '--text-output', 
                  dest='text_output', default='dist.txt',
                  help='Filename for text output (default: dist.txt)')
parser.add_option('-s', '--sort-files', dest='sort_files', action='store_true',
                  help='Sort trajectory files')
parser.add_option('-N', '--sort-numerically', dest='sort_numerically',
                  action='store_true',
                  help='Sort trajectory files numerically rather than '
                      +'lexicographically')
parser.add_option('-p', '--filename-pattern', dest='filename_pattern',
                  default=r'.*?(\d+)\.xtc$',
                  help='Regular expression used to extract sort key from '
                      +'filenames when sorting numerically. Must contain '
                      +'exactly one group matching an integer value. '
                      +"(Default: '.*?(\d+)\.xtc$')")
parser.add_option('-v', '--verbose', dest='verbose',
                  action = 'store_true',
                  help = 'Do not suppress child output')
(opts, args) = parser.parse_args()

if len(args) < 4:
    parser.print_help(sys.stderr)
    sys.exit(2)
else:
    try:
        group1 = int(args.pop(0))
        group2 = int(args.pop(0))
    except ValueError, e:
        sys.stderr.write('groups must be integers\n')
        parser.print_help(sys.stderr)
        sys.exit(2)
    xtc_files = args
        
if not opts.rc_file:
    sys.stderr.write('a GROMACS run control (.tpr) file must be provided\n')
    parser.print_help(sys.stderr)
    sys.exit(2)
    
if not opts.index_file:
    sys.stderr.write('a GROMACS index file (.ndx) must be provided\n')
    parser.print_help(sys.stderr)
    sys.exit(2)
    
if opts.sort_files or opts.sort_numerically:
    if opts.sort_numerically:
        re_get_key = re.compile(opts.filename_pattern)
        try:
            xtc_files = [item[1] for item in 
                         sorted((int(re_get_key.search(xtc_file).group(1)), 
                                 xtc_file)
                                for xtc_file in xtc_files)]
        except AttributeError:
            sys.stderr.write('at least one file did not match the key pattern\n')
            sys.exit(1)
    else:
        xtc_files = list(sorted(xtc_files))
        
sys.stdout.write(('Generating distances for groups %d and %d in the following '
                  +'files:\n') % (group1, group2))
for xtc_file in xtc_files:
    sys.stdout.write('  %s\n' % xtc_file)

txt_out_fmt = '%-20.15g    %-12.8g\n'
txt_out = open(opts.text_output, 'wt')

dt = opts.timestep
last_dist_txt = None
ti = -1

base_args = [opts.distance_program, '-noxvgr', '-s', opts.rc_file, '-n', opts.index_file]

if opts.distance_program == 'g_mindist':
    output_arg = '-od'
else:
    output_arg = '-o'

for (i, trajfile) in enumerate(xtc_files):
    basename = os.path.basename(trajfile[:-4])
    distfilename = 'dist-' + basename + '.xvg'    
    sys.stdout.write('%s %s -> %s\n' % (opts.distance_program, trajfile, 
                                        distfilename))
    sys.stdout.flush()
    if opts.verbose:
        child_output = None
    else:
        child_output = open('/dev/null', 'wb')
    proc = subprocess.Popen(base_args + ['-f', trajfile, 
                                         output_arg, distfilename],
                            stdin = subprocess.PIPE,
                            stdout = child_output,
                            stderr = subprocess.STDOUT)                  
    proc.stdin.write('%s %s\n' % (group1, group2))
    rc = proc.wait()
    
    if rc != 0:
        sys.stderr.write('%s exited with code %r; aborting\n' 
                         % (opts.distance_program, rc))
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

        txt_out.write(txt_out_fmt % (t, dist))

    txt_out.flush()
    dist_in.close()

txt_out.close()
