#
# This is an example script for performing interactive HSQCcos.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: May 22, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# HSQCcos method citation:
# Rudd TR, Macchi E, Muzi L, Ferro M, Gaudesi D, Torri G, Casu B, 
# Guerrini M, Yates EA. 
# Unravelling Structural Information from Complex Mixtures Utilizing 
# Correlation Spectroscopy Applied to HSQC Spectra. 
# Anal Chem. 2013 Aug 6;85(15):7487–7493. 
# https://doi.org/10.1021/ac4014379
#

# POKY libraries
import __main__
s = __main__.main_session

import numpy as np
import nmrglue as ng
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

specnames = s.show_spectrumselectiondialog('Select multiple spectra', 1)
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

# Get min-max of spectra for normalization
range_list = []
for spec in spec_list:        
  dic, data = ng.sparky.read(spec.data_path)
  #dic, data = ng.sparky.read_lowmem(spec.data_path)
  dmax, dmin = np.max(data), np.min(data)
  range_list.append([dmin, dmax, dmax - dmin])
  
# Set first peak as the reference peak
peaks = ref_spec.peak_list()
peaks.remove(ref_peak)
peaks = [ref_peak, ] + peaks

# Also make sure ref_spec is first.
spec_list.remove(ref_spec)
spec_list = [ref_spec, ] + spec_list

# Read the data in, and mean center and pareto scale
for i in range(len(peaks)):
  peak = peaks[i]  
  hts_list = []
  for j in range(len(spec_list)):
    spec = spec_list[j]
    # Normalized height
    hts = (spec.data_height(peak.frequency) - range_list[j][0]) \
            / range_list[j][2] 
    hts_list.append(hts)
  if i == 0:
    data_list = np.array(hts_list)
  else:
    data_list = np.vstack( [data_list, hts_list] )
  
# mean centered and pareto scale
for i in range(len(data_list[0])):
  # mean center
  avg = np.average(data_list[:, i])
  mc_data = np.subtract(data_list[:, i], avg) 
  
  # pareto scale
  std = np.std(mc_data)
  if std != 0:
    data_list[:, i] = np.divide(mc_data, np.sqrt(std))

# R calculation against the reference peak
# Select the peak if over R cutoff
print('Reference peak %.3f, %.3f' % (ref_peak.frequency))
print('------------------------------------------------------')
ref_data = data_list[0]
for i in range(1, len(data_list)):
  peak = peaks[i]
  data = data_list[i]
  #peak, data = peak_and_data
  R = np.corrcoef(ref_data, data)
  if abs(R[0][1]) >= Rcutoff:
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