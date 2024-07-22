#
# This is an example script to predict Kd from NMR titration data.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Kd prediction by CSP or intensity/volume ratios')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project

# PARAMETERS
# place spectrum names, protein conc, and ligand conc.
# Molar ratio can be used and base conc. can be set below.
spec_list = [
            ['protein_ligand_1_0',  1, 0],
            ['protein_ligand_1_1',  1, 1],
            ['protein_ligand_1_2',  1, 2],
            ['protein_ligand_1_4',  1, 4],
            ['protein_ligand_1_8',  1, 8],
            ]

# Real protein concentration in case ratio used above
protein_conc = 500
unit = 'uM'

# the number of equivalent binding sites
N = 1

# The number of fitting trials
ntrial = 1000

# Initial guess for Kd
Kd_guess = 200

# Mode: 'CS', 'I', or 'V'
diff_mode = 'CS'

print('Protein base concentration: %f %s' % (protein_conc, unit))

# If you want to make a plot, set True below
plot = True

# If you are running many peaks, high number can be risky
max_peaks_for_ind_plot = 10 # maximum peak counts for individual plots

##################################################
if len(s.selected_peaks()) < 1:
  s.show_message('Error', 'At least one peak should be selected.')
  quit()
print('Experiments used:')
for sp in spec_list: print(sp)

import sys
from sputil import name_to_spectrum, find_peak_by_assignment, sort_peaks_by_assignment
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as plt

sp_list = list(map(lambda sp: name_to_spectrum(sp[0], s), spec_list))

# Calculate euclidean distance between two peaks
def euc_dist(peak, peak2, ratio):
  d1 = (peak.frequency[0] - peak2.frequency[0])**2
  d2 = (peak.frequency[1] - peak2.frequency[1])**2
  if peak.frequency[0] > peak.frequency[1]:
    return (d1 + d2 / ratio / ratio) ** .5
  return (d2 + d1 / ratio / ratio) ** .5

# Calculate intensity ratio between two peaks
# mode 0: intensity, mode 1: volume
def Iratio(peak, peak2, mode=0):
  try:
    if mode == 0:
      return peak2.data_height / peak.data_height
    elif mode == 1:
      return peak2.volume / peak.volume
  except:
    print('This peak failed.')
    return 0.
try:
  none_idx = sp_list.index(None)
  s.show_message('Error', '%s spectrum could not be found.' %
          (spec_list[none_idx][0]))
  sys.exit()
except:
  pass

peaks = sort_peaks_by_assignment(s.selected_peaks(), False)


print('-------------------------')
print('Seq. %9s %9s' % ('Kd', 'Err.'))
print('-------------------------')
Kd_list, E_list, N_list = [], [], []
for ref_peak in peaks:
  peak_list = list(map(lambda sp:
    find_peak_by_assignment(sp, ref_peak.assignment), sp_list))

  try:
    none_idx = peak_list.index(None)
    s.show_message('Error', '%s assignment not found in %s' %
            (ref_peak.assignment, sp_list[none_idx].name))
    sys.exit()
  except:
    pass

  # Reinitialize a holo-peak
  ref_peak = peak_list[0]
  nseq = ref_peak.resonances().group.number
  print(f'{ref_peak.frequency}, {nseq}')

  Ps = np.array(list(map(lambda sp: sp[1]*protein_conc, spec_list)))
  Ls = np.array(list(map(lambda sp: sp[2]*protein_conc, spec_list)))
  P2Ls = np.array(list(map(lambda sp: sp[2]/sp[1], spec_list)))

  # 15N to 1H ratio. Gly 0.2, Others 0.14
  B = 0.14
  if ref_peak.resonances()[0].group.symbol == 'G':
    B = 0.2

  if diff_mode == 'CS': 
    # Calculate distances over all experiments
    dCSs = np.array(list(map(lambda p: euc_dist(ref_peak, p, B), peak_list)))
  elif diff_mode == 'I':
    # Calculate intensity ratios over all experiments
    dCSs = np.array(list(map(lambda p: Iratio(ref_peak, p, 0), peak_list)))
  elif diff_mode == 'V':
    # Calculate intensity ratios over all experiments
    dCSs = np.array(list(map(lambda p: Iratio(ref_peak, p, 1), peak_list)))

  # Maximum difference
  dCSmax = max(dCSs)

  # Model Becker et al. Kd model
  # dCS = dCSmax [ (N Pt + Lt + Kd) - [ (N Pt + Lt + Kd)^2 - 4 N Pt Lt ]^.5 ] / 2 N Pt
  def func(Kd, Pt, Lt):
    dCS = dCSmax * ((N*Pt+Lt+Kd) - ((N*Pt+Lt+Kd)**2 - 4*N*Pt*Lt)**.5) / (2*N*Pt)
    return dCS

  # Subtract from model and return residual
  def residual(Kd, Pt, Lt, dCS):
      return dCS - func(Kd, Pt, Lt)


  Kd, cov, infodict, mesg, ier = optimize.leastsq(
    residual, Kd_guess, args=(Ps, Ls, dCSs), full_output=True)

  E = np.sqrt(np.diag(cov))
  Kd_list.append(Kd)
  E_list.append(E) 
  N_list.append(nseq)

  print('%4d %9.4f %9.4f' % (nseq, Kd, E))
  
  s.show_message('Info', 'Dissociation constant (Kd): %.4f %s predicted from %s' %
                  (Kd[0], unit, ref_peak.assignment))
  print('Dissociation constant (Kd): %.4f %s predicted from %s' %
          (Kd[0], unit, ref_peak.assignment))

  # individual plotting. only less than ten.
  if plot and len(peaks) < max_peaks_for_ind_plot:
    xlabel = 'Molar Ratio'
    if diff_mode == 'CS':
      ylabel = 'CSP (ppm)'
    elif diff_mode == 'I':
      ylabel = 'Iratio'
    elif diff_mode == 'V':
      ylabel = 'Vratio'
    title = ref_peak.assignment
    plt.figure()
    plt.plot(P2Ls, dCSs, 'bo', label='data', markersize=5)
    plt.plot(P2Ls, func(Kd, Ps, Ls), 'b-',
         label='fit: Kd=%.4f %s' % (Kd, unit))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.pause(0.1)
    plt.show(block=False)

minKd, maxKd = min(Kd_list), max(Kd_list)
avgKd, stdKd = sum(Kd_list) / len(Kd_list), np.std(Kd_list)

if len(Kd_list) != 0:
  print('-------------------------')
  print('Max: %.4f %s' % (maxKd, unit))
  print('Min: %.4f %s' % (minKd, unit))
  print('Avg: %.4f %s (+/- %.4f)' % (avgKd, unit, stdKd))
  print('-------------------------')

plt.figure()
plt.errorbar(N_list, Kd_list, yerr = E_list, fmt='O', capsize=5)
plt.xlabel('Sequence')
plt.ylabel(f'Kd ({unit})')
plt.axhline(y=avgKd, color='g', linestyle='--')
plt.title('POKY Kd Plot')
plt.pause(0.1)
plt.show(block=False)
