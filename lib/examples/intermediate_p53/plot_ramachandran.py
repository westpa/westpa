#!/usr/bin/env python
import numpy
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors

phi = []
psi = []
for i in range(1,41):
  string = "iterations/iter_" + str(i).zfill(8) + "/auxdata/dihedral_2"
  f = h5py.File("west.h5")
  angles = f[string]
  phi = phi + list(angles[:,-1,0])
  psi = psi + list(angles[:,-1,1])
plt.hist2d(phi, psi, bins=100, norm=mcolors.PowerNorm(0.30), cmap='plasma')
plt.xlabel(r'$\phi$ (Degrees)')
plt.ylabel(r'$\psi$ (Degrees)')
#plt.title("Ramachandran Plot of p53 Peptide")
plt.colorbar()
plt.savefig("p53_ramachandran.pdf", dpi=600)
