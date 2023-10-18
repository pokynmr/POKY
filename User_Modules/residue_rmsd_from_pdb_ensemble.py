#
# This is an example script to calculate BB & HA RMSD per residue from PDB ensemble file.
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
                            'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')
if in_file == '':
  raise SystemExit

import os.path
if not os.path.exists(in_file):
  raise SystemExit

from pymol import cmd

cmd.delete('for_rmsd')
cmd.load(in_file, 'for_rmsd')
model = cmd.get_model('for_rmsd')
residues = model.get_residues()
first_resi = int(model.atom[0].resi)
last_resi = first_resi + len(residues)

# Backbone
bb_rmsds = []
for i in range(first_resi, last_resi + 1):
  fit_result = cmd.intra_rms(f'for_rmsd & i. {i} & n. N+CA+C+O', 1)
  bb_rmsds.append(fit_result[1:])

# All heavy atoms
ha_rmsds = []
for i in range(first_resi, last_resi + 1):
  fit_result = cmd.intra_rms(f'for_rmsd & i. {i} & n. N*+C*+O*+S*', 1)
  ha_rmsds.append(fit_result[1:])

import numpy as np

print('\n\n* POKY ENSEMBLE RESIDUE RMSD')
print('* PDB FILE: ' + in_file)
print('   RESIDUE    BB_RMSD_MEAN    BB_RMSD_STD    HA_RMSD_MEAN    HA_RMSD_STD')

for i in range(len(residues)):
  bb_rmsd_mean = np.mean(bb_rmsds[i])
  bb_rmsd_std = np.std(bb_rmsds[i])
  ha_rmsd_mean = np.mean(ha_rmsds[i])
  ha_rmsd_std = np.std(ha_rmsds[i])

  m = '%10d    ' % (first_resi + i)
  m += '%12.3f   ' % (bb_rmsd_mean)
  m += '%12.3f    ' % (bb_rmsd_std)
  m += '%12.3f   ' % (ha_rmsd_mean)
  m += '%12.3f   ' % (ha_rmsd_std)
  print(m)

cmd.delete('for_rmsd')