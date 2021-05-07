#
# This is an example script to calculate T1/T2/TauC from relaxation data.
# Also approximate M.W. from TauC for a rigid protein.
# This code will calculate TauC per residue unlike T1_T2_TauC_script.py.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY T1/T2/TauC calculation')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

#################
# User parameters

# Two ways to set experiments.
# 1. Set manually
# set experiments for T1 (ms)
T1_spec_list = [ # place spectrum names and (time) parameters. e.g.
            #['protein_0_1_T1',  0],
            #['protein_50_2_T1',  50],
            #['protein_100_3_T1',  100],
            #['protein_200_4_T1',  200],
            #['protein_300_5_T1',  300],
            ]
# set experiments for T2 (ms)
T2_spec_list = [ # place spectrum names and (time) parameters. e.g.
            #['protein_10_1_T2',  10],
            #['protein_30_2_T2',  30],
            #['protein_50_3_T2',  50],
            #['protein_70_4_T2',  70],
            #['protein_100_5_T2',  100],
            ]

# 2. Automatic setup- spectrum name must include the time parameter
T1_names = s.show_spectrumselectiondialog('Select T1 spectra', 1)
T2_names = s.show_spectrumselectiondialog('Select T2 spectra', 1)
T1_list = T1_names.split('\t')
T2_list = T2_names.split('\t')
import re
for specname in T1_list:
  t = re.search('[0-9]+', specname)
  T1_spec_list.append([specname, int(t.group(0))])
T1_spec_list.sort(key = lambda x: x[1])
for specname in T2_list:
  t = re.search('[0-9]+', specname)
  T2_spec_list.append([specname, int(t.group(0))])
T2_spec_list.sort(key = lambda x: x[1])
###############

# calculate T1/T2
def exp_func(x, a, b):
    return a * np.exp(-1 * x / b)

from sputil import name_to_spectrum, sort_peaks_by_assignment
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

# (time) parameters
xdataT1 = np.array(list(map(lambda sp: sp[1], T1_spec_list)))
xdataT2 = np.array(list(map(lambda sp: sp[1], T2_spec_list)))

# get spectrum instances first
spT1_list = list(map(lambda sp: name_to_spectrum(sp[0], s), T1_spec_list))
spT2_list = list(map(lambda sp: name_to_spectrum(sp[0], s), T2_spec_list))

if None in spT1_list:
  none_name = T1_spec_list[spT1_list.index(None)][0]
  print('T1 experiment (' + none_name + ') does not exist in your spectrum list.')
  raise SystemExit

if None in spT2_list:
  none_name = T2_spec_list[spT2_list.index(None)][0]
  print('T2 experiment (' + none_name + ') does not exist in your spectrum list.')
  raise SystemExit

# Calculation function definition
def calcT(sp_list, xdata, t1t2):
  ref_peaks = spT1_list[0].peak_list()
  sorted_peaks = sort_peaks_by_assignment(ref_peaks, False)
  fit_result_list, tc_used_list = [], []
  for ref_peak in sorted_peaks:
    pos = ref_peak.position

    # residue number
    nres = ref_peak.resonances()[0].group.number
    # peak intensities
    ydata = np.array(list(map(lambda sp: sp.data_height(pos), sp_list)))

    # fitting
    popt, pcov = curve_fit(exp_func, xdata, ydata,
                          bounds=([max(ydata), 0], [10*max(ydata), 10000]))
    point_sd = list(np.sqrt(np.diag(pcov)))
    fit_result_list.append([nres, popt[1], point_sd[1], ref_peak.assignment])

  # print out results
  print(t1t2 + ' Relaxation Results')
  print('%-5s %-12s %-12s %-16s' % \
    ('SeqID', t1t2 + ' decay', 'Deviation', 'Assignment') )
  for fr in fit_result_list:
    line = '%-5d %-12.3f %-12.3f %-16s' % (fr[0], fr[1], fr[2], fr[3])
    print(line)

  xdata2 = np.array(list(map(lambda y: y[0], fit_result_list)))
  ydata = np.array(list(map(lambda y: y[1], fit_result_list)))
  ysddata = np.array(list(map(lambda y: y[2], fit_result_list)))
  return xdata2, ydata, ysddata, np.average(ydata), np.std(ydata), \
                    fit_result_list

# Plotting function definition
def plotT(x, y, y2, t1t2):
  xlabel = 'Residue Number'
  ylabel = '%s decay rate' % (t1t2)
  title = '%s relaxation' % (t1t2)

  plt.figure()
  plt.errorbar(x, y, yerr=y2, fmt='bo', markersize=5)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.ylim([0, max(y+y2) * 1.1])
  plt.pause(0.1)
  plt.show(block=False)

####
# calculate T1 and plot if requested
x, y, y2, T1avg, T1std, T1all = calcT(spT1_list, xdataT1, 'T1')
plotT(x, y, y2, 'T1')

# calculate T2 and plot if requested
x, y, y2, T2avg, T2std, T2all = calcT(spT2_list, xdataT2, 'T2')
plotT(x, y, y2, 'T2')


# http://www.nmr2.buffalo.edu/nesg.wiki/NMR_determined_Rotational_correlation_time
# NESG's empirical correlation
# between TauC (ns) and M.W. (kDa) for rigid proteins.
#
# TauC = 0.58 * MW + 0.472

# Tc (ns) = 1 / (4*pi*Nresonance_freq) * (6*T1/T2 - 7)**.5 * 1000
# This means, when TauC is determined, MW can be predicted.

T1all_dict = {x[0]:x[1:] for x in T1all} 
T2all_dict = {x[0]:x[1:] for x in T2all} 

nucidx = spT1_list[0].nuclei.index('15N')
Nfreq = spT1_list[0].hz_per_ppm[nucidx]

nres1 = min([x[0] for x in T1all])
nres2 = max([x[0] for x in T1all])

print('\n* TauC/MW Results')
print('%-5s %-12s %-12s' % ('SeqID', 'TauC (ns)', 'M.W. (kDa)') )
res_list, TauC_list, MW_list = [], [], [] 
for i in range(nres1, nres2+1):
  try:
    t1data = T1all_dict[i]
    t2data = T2all_dict[i]
  except:
    continue

  TauC = 1 / (4 * 3.141592653589793238 * Nfreq) * \
          (6 * t1data[0]/t2data[0] - 7) ** 0.5 * 1000
  MW = (TauC - 0.472) / 0.58
  res_list.append(i)
  TauC_list.append(TauC)
  MW_list.append(MW)
  line = '%-5d %-12.3f %-12.3f' % (i, TauC, MW)
  print(line)

TauCavg = np.average(TauC_list)
TauCstd = np.std(TauC_list)
MWavg = np.average(MW_list)
MWstd = np.std(MW_list)

print('\nNitrogen resonance: %f.3f' % (Nfreq))
print('T1 average: %.3f (+/-%.3f)' % (T1avg, T1std))
print('T2 average: %.3f (+/-%.3f)' % (T2avg, T2std))
print('TauC (ns): %.3f (+/-%.3f)' % (TauCavg, TauCstd))
print('Approx. M.W. (kDa): %.3f (+/-%.3f)' % (MWavg, MWstd))

plt.figure()
plt.plot(res_list, TauC_list, 'bo')
plt.xlabel('Residue Number')
plt.ylabel('TauC (ns)')
plt.title('TauC prediction')
plt.ylim([0, max(TauC_list) * 1.1])
plt.pause(0.1)
plt.show(block=False)

plt.figure()
plt.plot(res_list, MW_list, 'ro')
plt.xlabel('Residue Number')
plt.ylabel('Molecular Weight (kDa)')
plt.title('M.W. prediction')
plt.ylim([0, max(MW_list) * 1.1])
plt.pause(0.1)
plt.show(block=False)
