#
# This is an example script to calculate dihedral angles from a PDB.
# This script runs on POKY BUILD 02/06/21e or newer
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

in_file = s.open_filedialog('Select the PDB ensemble file', 
                            'mmCIF (*.cif);; PDB (*.pdb);; Any (*)', '')
if in_file == '':
  raise SystemExit

import os.path
if not os.path.exists(in_file):
  raise SystemExit

from pymol import cmd

cmd.delete('for_dihe')
cmd.load(in_file, 'for_dihe')
nmodel = cmd.count_states('for_dihe')
cmd.split_states('for_dihe')

# get min, max residue number
nmin, nmax = 10**5, -10**5
atom_list = cmd.get_model("for_dihe & name CA").atom
for atm in atom_list:
  nmin = min(int(atm.resi), nmin)
  nmax = max(int(atm.resi), nmax)

import numpy as np
phi_angles = np.zeros((nmodel, nmax-nmin))
psi_angles = np.zeros((nmodel, nmax-nmin))
phi_means = np.zeros(nmax-nmin)
psi_means = np.zeros(nmax-nmin)
phi_devs = np.zeros(nmax-nmin)
psi_devs = np.zeros(nmax-nmin)

# phi/psi angle.
for i in range(1, nmodel+1):
  for j in range(nmin+1, nmax): # ignore N', C'
    prefix = 'for_dihe///' + str(j-1)
    prefix2 = 'for_dihe///' + str(j)
    prefix3 = 'for_dihe///' + str(j+1)
    phi_angles[i-1, j-nmin] = cmd.get_dihedral(prefix + '/C', 
      prefix2 + '/N', prefix2 + '/CA', prefix2 + '/C', state=i)
    psi_angles[i-1, j-nmin] = cmd.get_dihedral(prefix2 + '/N', 
      prefix2 + '/CA', prefix2 + '/C', prefix3 + '/N', state=i)
cmd.delete('for_dihe*')

# calculate deviation
from scipy.stats import circstd, circmean
for i in range(nmin+1, nmax):
  phi_means[i-nmin] = circmean(phi_angles[:, i-nmin], high=180, low=-180)
  psi_means[i-nmin] = circmean(psi_angles[:, i-nmin], high=180, low=-180)
  phi_devs[i-nmin] = circstd(phi_angles[:, i-nmin], high=180, low=-180)
  psi_devs[i-nmin] = circstd(psi_angles[:, i-nmin], high=180, low=-180)

print('\n\n--------------------------------------------')
print('* POKY DIHEDRAL ANGLES')
print('* PDB FILE: ' + in_file)
print('--------------------------------------------')
for i in range(1, nmodel+1):
  print(f'* MODEL {i}')  
  print('ResID  PHI     PSI')
  for j in range(nmin+1, nmax):
    print('%4d - %8.3f %8.3f' % (i, phi_angles[i-1, j-nmin], 
                                 psi_angles[i-1, j-nmin]))
  print('--------------------------------------------')

if nmodel > 1:
  print(f'* MODEL ENSEMBLE')
  print('ResID  PHI MEAN PHI DEV. PSI MEAN PSI DEV. - ORDERNESS')
  olist = []
  for i in range(nmin+1, nmax):
    if phi_devs[i-nmin] > 24 or  psi_devs[i-nmin] > 24:
      postfix = 'not ordered'
    else:
      postfix = 'ordered'
      olist.append(i)
    print('%4d - %8.3f %8.3f %8.3f %8.3f - %s' % (i, 
                      phi_means[i-nmin], phi_devs[i-nmin], 
                      psi_means[i-nmin], psi_devs[i-nmin], postfix))

  print('--------------------------------------------')
  # format ordered regions
  from operator import itemgetter
  from itertools import groupby
  ranges = ''
  for k,g in groupby(enumerate(olist),lambda x:x[0]-x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int,group))
    if len(group) != 1:
      ranges += '%d-%d+' % (group[0], group[-1])
    else:
      ranges += '%d+' % (group[0])
  print('Ordered regions: ' + ranges[:-1])
