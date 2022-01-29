#
# This is an example script to calculate R1/R2 from relaxation data.
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
print('POKY R1/R2 calculation')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

#################
# User parameters
plotR1 = True # if R1 plotting preferred
plotR2 = True # if R2 plotting preferred

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

# calculate R1/R2
def exp_func(x, a, b):
    return a * np.exp(-0.001 * x * b)

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
def calcR(sp_list, xdata, r1r2):
  ref_peaks = spT1_list[0].peak_list()
  sorted_peaks = sort_peaks_by_assignment(ref_peaks, False)
  fit_result_list = []
  for ref_peak in sorted_peaks:
    pos = ref_peak.position

    # residue number
    nres = ref_peak.resonances()[0].group.number
    # peak intensities
    ydata = np.array(list(map(lambda sp: sp.data_height(pos), sp_list)))
    while np.mean(ydata) > 10**11:
      ydata = ydata ** .5
      
    # fitting
    popt, pcov = curve_fit(exp_func, xdata, ydata)
    point_sd = list(np.sqrt(np.diag(pcov)))
    fit_result_list.append([nres, popt[1], point_sd[1], ref_peak.assignment])
    
  # print out results
  print(r1r2 + ' Relaxation Results')
  print('%-5s %-12s %-12s %-16s' % \
    ('SeqID', r1r2, 'Deviation', 'Assignment') )
  for fr in fit_result_list:
    line = '%-5d %-12.3f %-12.3f %-16s' % (fr[0], fr[1], fr[2], fr[3])
    print(line)

  xdata2 = np.array(list(map(lambda y: y[0], fit_result_list)))
  ydata = np.array(list(map(lambda y: y[1], fit_result_list)))
  ysddata = np.array(list(map(lambda y: y[2], fit_result_list)))
  return xdata2, ydata, ysddata, np.average(ydata), np.std(ydata)

# Plotting function definition
def plotR(x, y, y2, r1r2):
  xlabel = 'Residue Number'
  ylabel = '%s' % (r1r2)
  title = '%s relaxation' % (r1r2)

  plt.figure()
  plt.errorbar(x, y, yerr=y2, fmt='bo', markersize=5)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.ylim([0, max(y+y2) * 1.1])
  plt.pause(0.1)
  plt.show(block=False)

####
# calculate R1 and plot if requested
x, y, y2, R1avg, R1std = calcR(spT1_list, xdataT1, 'R1')
if plotR1:
  plotR(x, y, y2, 'R1')

# calculate R2 and plot if requested
x, y, y2, R2avg, R2std = calcR(spT2_list, xdataT2, 'R2')
if plotR2:
  plotR(x, y, y2, 'R2')

print('R1 average (all): %.3f (+/-%.3f)' % (R1avg, R1std))
print('R2 average (all): %.3f (+/-%.3f)' % (R2avg, R2std))
