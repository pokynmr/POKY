#
# This is an example script to plot a ramachandran plot from a PDB.
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

from sputil import process_pymol_residue_range
# model number
if nmodel > 1:
  temp = s.show_inputdialog('Model number',
        f'Select the model numbers in PDB between 1 and {nmodel} (e.g. 1-{nmodel})',
        f'1-{nmodel}')
  temp2 = process_pymol_residue_range(temp)
  if temp2[0] == True:
    model_numbers = temp2[1]
else:
  model_numbers = [1,]

# range
residue_numbers = s.show_inputdialog('Residue number',
        f'Select the residue numbers between {nmin+1} and {nmax-1} (e.g. {nmin+1}-{nmax-1})',
        f'{nmin+1}-{nmax-1}')
temp = process_pymol_residue_range(residue_numbers)
if temp[0] == True:
  res_list = temp[1]
else:
  s.show_message('Error', 'Wrong residue numbers.')
  raise SystemError

annotate = s.show_message_yes_no('Annotate',
                                 'Do you want to annotate?')

# plot
import matplotlib.pyplot as plt
plt.xlabel("Phi ($\phi$)")
plt.ylabel("Psi ($\psi$)")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
c_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
for nm in model_numbers:
  for nres in res_list:
    phi = phi_angles[nm-1, nres-nmin]
    psi = psi_angles[nm-1, nres-nmin]
    plt.scatter(phi, psi, c=c_list[nm % len(c_list)])
    if annotate:
      plt.annotate(f'#{nm}:{nres}', [phi, psi])

plt.show(block=False)