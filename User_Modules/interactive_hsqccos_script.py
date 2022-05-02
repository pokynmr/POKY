#
# This is an example script for performing interactive HSQCcos.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: May 2, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

# POKY libraries
import __main__
s = __main__.main_session

import numpy as np
from sputil import name_to_spectrum

print('\n\n\n------------------------------------------------------')
print('POKY Interactive HSQCcos Analysis')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

if len(s.selected_peaks()) != 1:
  print('One peak from the reference spectrum needs to be selected.')
  raise SystemError

ref_peak = s.selected_peaks()[0]
ref_spec = ref_peak.spectrum

specnames = s.show_spectrumselectiondialog('Select a spectrum', 1)
specname_list = specnames.split('\t')
if len(specname_list) < 3:
  print('Select multiple spectra. Abort.')
  raise SystemError
spec_list = list(map(lambda x: name_to_spectrum(x, s), 
                              specname_list))

# Set R cutoff
try:
  Rcutoff = float(s.show_inputdialog('R cutoff','R cutoff: ','0.95'))
except:
  print('R cutoff must be a floating value.')
  raise SystemError

# Set first peak as the reference peak
peaks = [ref_peak, ]
for peak in ref_spec.peak_list():
  if peak not in peaks:
    peaks.append(peak)

data_list = []

# Read the data in, and mean center and pareto scale
for peak in peaks:
  data_list.append([peak, []])
  for spec in spec_list:
    hts = spec.data_height(peak.frequency)
    data_list[-1][1].append(hts)
  data_list[-1][1] = np.array(data_list[-1][1])
  # mean centered
  avg = np.average(data_list[-1][1])
  mc_data = np.subtract(data_list[-1][1], avg) 
  # Pareto scale
  std = np.std(mc_data)
  data_list[-1][1] = np.divide(mc_data, np.sqrt(std))

# R calculation against the reference peak
# Select the peak if over R cutoff
print('Reference peak %.3f, %.3f' % (ref_peak.frequency))
print('------------------------------------------------------')
ref_peak, ref_data = data_list[0]
for peak_and_data in data_list[1:]:
  peak, data = peak_and_data
  R = np.corrcoef(ref_data, data)
  if R[0][1] >= Rcutoff:
    peak.selected = 1
    state = 'selected'
  else:
    state = 'not selected'
  print('Peak %.3f, %.3f: R corrcoef %.3f- %s' % \
    (peak.frequency[0], peak.frequency[1], R[0][1], state))
    
print('------------------------------------------------------')  
cplbl = s.show_inputdialog('Label','Label selected.','Compound 1')
if cplbl != '':
  for peak in s.selected_peaks():
    peak.show_label(cplbl)
