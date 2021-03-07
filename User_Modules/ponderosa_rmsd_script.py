#
# This is an example script to calculate BB & HA RMSD from PONDEROSA results
# This script runs on POKY BUILD 02/06/21e or newer
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

pdir = s.open_directorydialog('Select the PONDEROSA result directory', '')
if pdir == '':
  raise SystemExit

import os.path
pbdir = os.path.join(pdir, 'BestEvaluated')
pdbpath = os.path.join(pbdir, 'final_water_refined.pdb')
cypath = os.path.join(pbdir, 'cyrange_region.txt_p')
if not os.path.exists(pdbpath):
  pdbpath = os.path.join(pbdir, 'final.pdb')

if not os.path.exists(pdbpath) or not os.path.exists(cypath):
  raise SystemExit # does not seem like a ponderosa result directory

f = open(cypath,'r')
ordered = f.read()
f.close()
oreg = ordered.replace(',', '+')    # to alleviate a small bug.
ordered = oreg.replace('+', ', ')   # more common in papers

from pymol import cmd

cmd.load(pdbpath, 'for_rmsd')
cmd.split_states('for_rmsd')

# Backbone
bb_rmsds = []
for i in range(2, 21):
  fit_result = cmd.align('for_rmsd_0001 & i. %s & n. N+CA+C+O' % (oreg), 
                        'for_rmsd_%04d' % (i), 
                        cutoff=0)
  bb_rmsds += [fit_result[0],]

# All heavy atoms
ha_rmsds = []
for i in range(2, 21):
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
print('* PDB FILE: ' + pdbpath)
print('* Ordered residues by CYRANGE: ' + ordered)
print('MODEL        BB RMSD  HA RMSD')

for i in range(2, 21):
  m = '%2d <-> 1    ' % (i)
  m += '%6.3f   ' % (bb_rmsds[i-2])
  m += '%6.3f ' % (ha_rmsds[i-2])
  print(m)
print('--------------------------------------------')
print('MEAN RMSD : %6.3f   %6.3f' % (bb_rmsd_mean, ha_rmsd_mean))
print('STDEV RMSD: %6.3f   %6.3f' % (bb_rmsd_std, ha_rmsd_std))
