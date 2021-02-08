#
# This is an example script to create a distance constraint file.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Distance Constraint Generator')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# Parameters
spec_names = s.show_spectrumselectiondialog('Select spectra to use', 1)
spec_list = spec_names.split('\t')

file_type = 'upl'           # format: upl for PONDEROSA/CYANA, tbl for XPLOR/CNS
only_selected = False       # use only selected peaks. If not, consider all peaks.
approximation_mode = 1/6    # approximation mode to be used. r^??? (e.g. 1/3, 1/4, 1/6)
reference_dist = 4.5        # NOE center calibration. It can be None. ssNMR carbons ~7
max_dist = 5.5              # Maximum distance. ssNMR carbons ~8
min_dist = 2.4              # Minimum distance. ssNMR carbons ~ 4.4?
vdw = 1.8                   # vdw distance. ssNMR carbons ~3.0
tolerance = 1.4             # Extra distance added to the calculated distance (1.4x default)
no_diagonal = True          # Exclude diagonal peaks?
res = {'1H': 0.02,
       '13C': 0.3,
       '15N': 0.3}          # Resolution used for excluding diagonal peaks

import __main__
s = __main__.main_session
proj = s.project

# where to save your file
outName = s.save_filedialog('Save your new restraint file as...',
            'DIANA UPL file (*.upl);; XPLOR TBL file (*.tbl);; Any (*)', '')

# Processing start
from sputil import name_to_spectrum, sort_peaks_by_assignment, parse_assignment_entirely
import numpy as np

A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

import nomenclature
upl_list = []
for sp in spec_list:
  spec = name_to_spectrum(sp, s)
  if spec == None:
    print('%s spectrum does not exist.' % (sp))
    raise SystemError

  dMin, dMax = 10**13, -10**13
  peaks = spec.spectrum().peak_list()
  if len(peaks) == 0: continue

  peak_list = []

  for peak in peaks:
    nuc = spec.nuclei
    freq = peak.frequency
    if peak.is_assigned == 0 or peak.data_height == None: continue
    if no_diagonal and abs(freq[0] - freq[1]) < res[nuc[0]]: continue
    if spec.spectrum().dimension == 2:
      dMax = max(dMax, abs(peak.data_height))
      dMin = min(dMin, abs(peak.data_height))
    elif spec.spectrum().dimension == 3:
      if no_diagonal and abs(freq[0] - freq[2]) < res[nuc[0]]: continue
      if no_diagonal and abs(freq[1] - freq[2]) < res[nuc[1]]: continue
      dMax = max(dMax, abs(peak.data_height))
      dMin = min(dMin, abs(peak.data_height))

    peak_list.append( [peak.resonances(), max(1, abs(peak.data_height)),] )

  if len(peak_list) == 0: continue

  # start calculating calibration constant
  if dMin == 0.0 or abs(-1.0+(dMax / dMin)**approximation_mode) < 0.0001):
    constant = 1
  else:
    constant = (max_dist - min_dist) / (-1.0 + (dMax / dMin)**approximation_mode)

  # calculate upper limits
  pairs, upl_limits = [], []
  for resn, height in peak_list:
    upl_limit = min_dist + (-1.0 + (dMax / height )**approximation_mode) * constant
    # process nomenclature
    atm = resn[0][2].replace('*','').replace('#','').replace('CQ', 'C')
    if atm[0] in ['Q', 'M']: atm = 'H' + atm[1:]
    atm2 = resn[1][2].replace('*','').replace('#','').replace('CQ', 'C')
    if atm2[0] in ['Q', 'M']: atm2 = 'H' + atm2[1:]
    try:
      atm3 = resn[1][2].replace('*','').replace('#','').replace('CQ', 'C')
      if atm3[0] in ['Q', 'M']: atm3 = 'H' + atm3[1:]
    except:
      pass
    # nuclei type should match. H/Q/M-H/Q/M vs C/CQ-C/CQ
    try:
      if nuc[0] == nuc[1] or '1H' not in nuc[0:2]:
        pairs.append( (A_dict[resn[0][0]], resn[0][0], resn[0][1], atm,
                      A_dict[resn[1][0]], resn[1][0], resn[1][1], atm2) )
        heights.append(height)
      elif nuc[0] == nuc[2] or '1H' not in [nuc[0], nuc[2]]:
        pairs.append( (A_dict[resn[0][0]], resn[0][0], resn[0][1], atm,
                      A_dict[resn[2][0]], resn[2][0], resn[2][1], atm2) )
        heights.append(height)
      elif nuc[1] == nuc[2] or '1H' not in [nuc[1], nuc[2]]:
        pairs.append( (A_dict[resn[1][0]], resn[1][0], resn[1][1], atm,
                      A_dict[resn[2][0]], resn[2][0], resn[2][1], atm2) )
        upl_limits.append(upl_limit)
    except:
      continue

  # reference to the distance set by a user
  if reference_dist != None:
    avg_upl = sum(upl_limits)/len(upl_limits)
    diff = reference_dist - avg_upl
    upl_limits = list(map(lambda x: x + diff, upl_limits)

  #add to upl_list
  for i in range(len(pairs)):
    pair, upl = pairs[i], upl_limit[i]

    if file_type == 'upl':
      line = '%3d %3s %4s %5d %3s %4s %8.2f\n' % (pair[2], pair[0], pair[3],
                                pair[6],pair[4], pair[7], upl * tolerance)
    elif file_type == 'tbl':
      pair[3] = nomenclature.translate(pair[1], pair[3], 'DIANA', 'XPLOR')
      pair[6] = nomenclature.translate(pair[5], pair[6], 'DIANA', 'XPLOR')
      if len(pair[3]) > 1: pair[3] += '*'
      if len(pair[6]) > 1: pair[6] += '*'

      line = 'assign ( resid %3d and name %4s ) ' % (pair[2], pair[3] )
      line += '( resid %3d and name %4s ) ' % (pair[5], pair[6] )
      line += '%5.2f %5.2f %5.2f' % (upl, upl - vdw, 0.)

    upl_list.append(line)

# write
f = open(outName, 'w')
for line in upl_list:
  f.write(line+'\n')
  print(line)
f.close()
