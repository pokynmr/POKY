#
# This is an example script to calculate T1/T2/TauC from relaxation data.
# Also approximate M.W. from TauC for a rigid protein.
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
plotT1 = True # if T1 plotting preferred
plotT2 = True # if T2 plotting preferred
plotMW = True # if approx. M.W. plotting preferred
minH = 8.5    # H ppm higher than this value will be considered for TauC      
maxSTD = 150  # T deviation less than this will be considered for TauC
excl_residues = [] # Exclude residues for TauC calculation.
                  # This can be like this, too.
                  # list(range(1, 10)) + list(range(60, 70))

# Assignment needed for the first spectrum of each list

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
    if min(pos) > minH and point_sd[1] < maxSTD and nres not in excl_residues:
      tc_used_list.append(popt[1])
    
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
                    np.average(tc_used_list), np.std(tc_used_list)

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
x, y, y2, T1avg, T1std, T1avg2, T1std2 = calcT(spT1_list, xdataT1, 'T1')
if plotT1:
  plotT(x, y, y2, 'T1')

# calculate T2 and plot if requested
x, y, y2, T2avg, T2std, T2avg2, T2std2 = calcT(spT2_list, xdataT2, 'T2')
if plotT2:
  plotT(x, y, y2, 'T2')

# http://www.nmr2.buffalo.edu/nesg.wiki/NMR_determined_Rotational_correlation_time
# NESG's empirical correlation
# between TauC (ns) and M.W. (kDa) for rigid proteins.
#
# TauC = 0.58 * MW + 0.472

# Tc (ns) = 1 / (4*pi*Nresonance_freq) * (6*T1/T2 - 7)**.5 * 1000
# This means, when TauC is determined, MW can be predicted.

nucidx = spT1_list[0].nuclei.index('15N')
Nfreq = spT1_list[0].hz_per_ppm[nucidx]
TauC = 1 / (4 * 3.141592653589793238 * Nfreq) * \
        (6 * T1avg2/T2avg2 - 7) ** 0.5 * 1000
MW = (TauC - 0.472) / 0.58
print('T1 average (all): %.3f (+/-%.3f)' % (T1avg, T1std))
print('T2 average (all): %.3f (+/-%.3f)' % (T2avg, T2std))
print('T1 average (>%.1fppm): %.3f (+/-%.3f)' % (minH, T1avg2, T1std2))
print('T2 average (>%.1fppm): %.3f (+/-%.3f)' % (minH, T2avg2, T2std2))

print('Nitrogen resonance: %f.3f' % (Nfreq))
print('TauC (ns): %f' % (TauC))
print('Approx. M.W. (kDa): %f' % (MW))

# Plot MW on NESG data
if plotMW:
  NESG_data = '''PsR76A  (NC5)	7.2	478.0	128.0	5.10
  VfR117 (NC)	11.2	605.0	119.0	6.30
  SyR11 (NC5)	12.4	630.0	104.0	7.10
  ER541-37-162 (NC5)	15.8	729.0	66.5	10.0
  ER540 (NC5)	18.8	909.0	66.5	11.3
  SoR190 (NC)	13.8	697.5	100.9	7.70
  TR80 (NC5)	10.5	612.8	102.9	7.00
  Ubiquitin (NC)	9.0	441.8	144.6	4.40
  HR2873B (NC)	10.7	492.0	115.0	5.70
  B-domain (NC)	7.2	423.5	153.3	4.05
  BcR97A (NC)	13.1	705.8	80.6	8.80
  PfR193A (NC)	13.6	733.9	80.9	9.00
  MvR76 (NC)	20.2	1015.0	64.5	12.2
  DvR115G (NC)	10.9	608.7	115.6	6.50
  MrR110B (NC5)	11.8	707.0	99.2	7.80
  VpR247 (NC5)	12.5	661.2	88.3	8.05
  BcR147A (NC)	11.9	645.0	104.0	7.20
  WR73 (NC5)	21.9	1261.0 41.3 13.0
  NsR431C (NC5)	16.8	855.5	71.2	10.6
  StR82 (NC)	9.2	537.3	100.4	6.6'''

  ln = NESG_data.split('\n')
  MW_list = list(map(lambda x: float(x.strip().split()[2]), ln))
  Tc_list = list(map(lambda x: float(x.strip().split()[5]), ln))
  xdata = np.array(sorted(MW_list))
  ydata = np.array(list(map(lambda i: Tc_list[MW_list.index(xdata[i])],
                      range(len(xdata)))))

  def func(x, a, b):
    return a * x + b
  popt, pcov = curve_fit(func, xdata, ydata)

  plt.figure()
  plt.plot(xdata, ydata, 'bo', label='NESG data')
  plt.plot(xdata, func(xdata, *popt), 'k-',
            label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
  plt.plot([MW,], [TauC], 'ro', label='%.3f kDa' % (MW), markersize = 7)
  plt.xlabel('MW (kDa)')
  plt.ylabel('TauC (ns)')
  plt.legend()
  plt.pause(0.1)
  plt.show(block=False)
