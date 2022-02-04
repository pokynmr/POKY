#
# This is an example script to extract Kex from SEA-TROSY data.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# SEA-TROSY (Solvent Exposed Amides with TROSY): A Method to Resolve 
# the Problem of Spectral Overlap in Very Large Proteins
# J. Am. Chem. Soc. 2001, 123, 19, 4633–4634
# Pellechia et al.
# https://doi.org/10.1021/ja005850t

import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY SEA-TROSY analysis')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

#################
# User parameters
plotK = True  # if R1 plotting preferred
minErr = 30   # minimum deviation requirement for plotting

# Assignment needed for the first spectrum of each list

# Two ways to set experiments.
# 1. Set manually
# set SEA-TROSY experiments (ms)
ST_spec_list = [ # place spectrum names and (time) parameters. e.g.
            #['protein_10_1_T1',  10],
            #['protein_20_2_T1',  20],
            #['protein_30_3_T1',  30],
            #['protein_40_4_T1',  40],
            #['protein_60_5_T1',  60],
            ]

# 2. Automatic setup- spectrum name must include the time parameter
ST_names = s.show_spectrumselectiondialog('Select SEA-TROSY spectra' + \
            ' The first number appears in each spectrum name' + \
            ' must be time parameter [ms].', 1)
ST_list = ST_names.split('\t')

import re
for specname in ST_list:
  t = re.search('[0-9]+', specname)
  ST_spec_list.append([specname, int(t.group(0))])
ST_spec_list.sort(key = lambda x: x[1])
###############

# calculate Kex
def exp_func(x, a, b, c):
    return a * (1 - b * np.exp(-1.0 * x / c))

from sputil import name_to_spectrum, sort_peaks_by_assignment
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

# (time) parameters
xdataST = np.array(list(map(lambda sp: sp[1], ST_spec_list)))

# get spectrum instances first
spST_list = list(map(lambda sp: name_to_spectrum(sp[0], s), ST_spec_list))

if None in spST_list:
  none_name = ST_spec_list[spST_list.index(None)][0]
  print('SEA-TROSY experiment (' + none_name + ') does not exist in your spectrum list.')
  raise SystemExit

# Calculation function definition
def calcKex(sp_list, xdata):
  ref_peaks = spST_list[0].peak_list()
  sorted_peaks = sort_peaks_by_assignment(ref_peaks, False)
  fit_result_list, fit_result_list2 = [], []
  for ref_peak in sorted_peaks:
    pos = ref_peak.position

    # residue number
    nres = ref_peak.resonances()[0].group.number
    # peak intensities
    ydata = np.array(list(map(lambda sp: sp.data_height(pos), sp_list)))
    while np.mean(ydata) > 10**11:
      ydata = ydata ** .5
    # fitting
    popt, pcov = curve_fit(exp_func, xdata, ydata, 
                       bounds=[[max(ydata), 0, 0], [10*max(ydata), 2, 1000] ],
                       maxfev=1000000000)
    point_sd = list(np.sqrt(np.diag(pcov)))
    fit_result_list.append([nres, popt, point_sd, ref_peak.assignment])
    if point_sd[2] < minErr:
      fit_result_list2.append([nres, popt, point_sd, ref_peak.assignment])   
    
  # print out results
  print('SEA-TROSY Exchange Rates')
  print('%-5s %-12s %-12s %-12s %-12s %-16s' % \
    ('SeqID', 'Kex', 'Deviation', 'I0', 'B', 'Assignment') )
  for fr in fit_result_list:
    line = '%-5d %-12.3f %-12.3f %-12.3f %-12.3f %-16s' % \
            (fr[0], fr[1][2], fr[2][2], fr[1][0], fr[1][1], fr[3])
    print(line)

  xdata2 = np.array(list(map(lambda y: y[0], fit_result_list)))
  ydata = np.array(list(map(lambda y: y[1][2], fit_result_list)))
  ysddata = np.array(list(map(lambda y: y[2][2], fit_result_list)))
  
  xpdata2 = np.array(list(map(lambda y: y[0], fit_result_list2)))
  ypdata = np.array(list(map(lambda y: y[1][2], fit_result_list2)))
  ypsddata = np.array(list(map(lambda y: y[2][2], fit_result_list2)))
  
  return xdata2, ydata, ysddata, xpdata2, ypdata, ypsddata

# Plotting function definition
def plotKex(x, y, y2):
  xlabel = 'Residue Number'
  ylabel = 'Kex (ms)'
  title = 'SEA-TROSY Exchange Rates'

  plt.figure()
  plt.errorbar(x, y, yerr=y2, fmt='bo', markersize=5)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.ylim([0, max(y+y2) * 1.1])
  plt.pause(0.1)
  plt.show(block=False)

####
# calculate Kex and plot if requested
x, y, y2, xp, yp, yp2 = calcKex(spST_list, xdataST)
if plotK:
  plotKex(xp, yp, yp2)
