#!/usr/bin/env python2.7

import sys, argparse

parser = argparse.ArgumentParser(description='''\
separate two residues by moving all atoms of the first by a specified amount
in the negative x direction''')
parser.add_argument('-o', '--output', default='sep.gro', help='output file (default %(default)s).')
parser.add_argument('distance', type=float, help='new separation (in A)')
parser.add_argument('input', help='input data')
args = parser.parse_args()

if args.output == '-':
    outfile = sys.stdout
else:
    outfile = file(args.output, 'wt')
    
if args.input == '-':
    infile = sys.stdin
else:
    infile = file(args.input, 'rt')
    
# Read and write title card
outfile.write(infile.readline())

# Read and write number of atoms
line = infile.readline()
outfile.write(line)
n_atoms = int(line.strip())
resnum = 0
for iatom in xrange(n_atoms):
    line = infile.readline()
    fields = [int(line[:5]),line[5:10],line[10:15],int(line[15:20]),
              float(line[20:28]),float(line[28:36]),float(line[36:44])]
    
    if iatom == 0:
        resnum = fields[0]
    
    if fields[0] == resnum:
        fields[4] -= args.distance/10.0
        
    outfile.write('{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(*fields))

# Write box info
outfile.write(infile.readline())

outfile.close()
infile.close()

