#
# This is an example script to calculate BB & HA RMSD from PDB ensemble file.
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

ordered = s.show_inputdialog('Ordered residues', 
                  'Type ordered residues (e.g. 2-5+10-30)', '2-5+10-30')

oreg = ordered.replace(',', '+')    # to alleviate a small bug.
ordered = oreg.replace('+', ', ')   # more common in papers

from pymol import cmd

cmd.delete('for_rmsd')
cmd.load(in_file, 'for_rmsd')
nmodel = cmd.count_states('for_rmsd')
cmd.split_states('for_rmsd')

# Backbone
bb_rmsds = []
for i in range(2, nmodel+1):
  fit_result = cmd.align('for_rmsd_0001 & i. %s & n. N+CA+C+O' % (oreg), 
                        'for_rmsd_%04d' % (i), 
                        cutoff=0)
  bb_rmsds += [fit_result[0],]

# All heavy atoms
ha_rmsds = []
for i in range(2, nmodel+1):
  fit_result = cmd.align('for_rmsd_0001 & i. %s & n. N*+C*+O*+S*' % (oreg), 
                        'for_rmsd_%04d' % (i), 
                        cutoff=0)
  ha_rmsds += [fit_result[0],]

import numpy as np
bb_rmsd_mean = np.mean(bb_rmsds)
bb_rmsd_std = np.std(bb_rmsds)
ha_rmsd_mean = np.mean(ha_rmsds)
ha_rmsd_std = np.std(ha_rmsds)

print('\n\n* POKY ENSEMBLE RMSD')
print('* PDB FILE: ' + in_file)
print('* Ordered residues by CYRANGE: ' + ordered)
print('MODEL        BB RMSD  HA RMSD')

for i in range(2, nmodel+1):
  m = '%2d <-> 1    ' % (i)
  m += '%6.3f   ' % (bb_rmsds[i-2])
  m += '%6.3f ' % (ha_rmsds[i-2])
  print(m)
print('--------------------------------------------')
print('MEAN RMSD : %6.3f   %6.3f' % (bb_rmsd_mean, ha_rmsd_mean))
print('STDEV RMSD: %6.3f   %6.3f' % (bb_rmsd_std, ha_rmsd_std))

if s.show_message_yes_no('Save', 'Do you want to save the superimposed structure?'):
  out_path = s.save_filedialog('Save as', 'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')
  if out_path == '':
    raise SystemExit
  cmd.intra_fit('for_rmsd & i. %s & n. N*+C*+O*+S*' % (oreg))
  cmd.save(out_path, 'for_rmsd', state=0)
  s.show_message('Finished', 'Finished.')

cmd.delete('for_rmsd*')