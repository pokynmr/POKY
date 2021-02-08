#
# This is an example script to show an T1/T2/SEA line chart from assigned peaks.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY T-decay bar chart by residue')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project

spec_list = [ # place spectrum names and (time) parameters. e.g.
             #['relaxT_10_1_N15_T2', 10],
             #['relaxT_30_2_N15_T2', 30],
             #['relaxT_50_3_N15_T2', 50],
             #['relaxT_70_4_N15_T2', 70],
             #['relaxT_90_5_N15_T2', 90],
             #['relaxT_110_6_N15_T2', 110],
            ]

def t_exp_func(x, a, b):
    return a * np.exp(-1 * x / b)

def sea_exp_func(x, a, b):
    return a * (1 - np.exp(-1 * x / b))

exp_func = t_exp_func
#exp_func = sea_exp_func # if sea-hsqc experiments

from sputil import name_to_spectrum
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

# (time) parameters
xdata = np.array(list(map(lambda sp: sp[1], spec_list)))

# get spectrum instances first
sp_list = list(map(lambda sp: name_to_spectrum(sp[0], s), spec_list))

if None in sp_list:
  none_name = spec_list[sp_list.index(None)][0]
  print(none_name + ' does not exist in your spectrum list.')
  raise SystemExit

# collect reference peaks and positions from the first spectrum
# assignment needed!
ref_peaks = sp_list[0].peak_list()
# ref_peaks = sp_list[0].selected_peaks() # if only selected peaks preferred

# sort peaks by assignment
from sputil import sort_peaks_by_assignment
sorted_peaks = sort_peaks_by_assignment(ref_peaks, False)

fit_result_list = []
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

xdata = np.array(list(map(lambda y: y[0], fit_result_list)))
ydata = np.array(list(map(lambda y: y[1], fit_result_list)))
ysddata = np.array(list(map(lambda y: y[2], fit_result_list)))

# plotting
xlabel = 'Residue Number'
ylabel = 'T-decay rate'
title = 'T-relaxation'

plt.figure()
plt.errorbar(xdata, ydata, yerr=ysddata, fmt='bo', markersize=5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.pause(0.1)
plt.show(block=False)

# print out results
print('%-5s %-12s %-12s %-16s' % ('SeqID', 'T-decay', 'Deviation', 'Assignment') )
for fr in fit_result_list:
  line = '%-5d %-12.3f %-12.3f %-16s' % (fr[0], fr[1], fr[2], fr[3])
  print(line)
