#
# This is an example script to show an exponential fitting of a selected peak.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Exponential Fitting')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project


if len(s.selected_peaks()) != 1:
  print('One peak should be selected.')
  quit()

spec_list = [ # place spectrum names and (time) parameters. e.g.
            # ['relaxT_10_1_N15_T2', 10],
            # ['relaxT_30_2_N15_T2', 30],
            # ['relaxT_50_3_N15_T2', 50],
            # ['relaxT_70_4_N15_T2', 70],
            # ['relaxT_90_5_N15_T2', 90],
            # ['relaxT_110_6_N15_T2', 110],
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

# set reference peak and position
ref_peak = s.selected_peaks()[0]
pos = ref_peak.position

# (time) parameters
xdata = np.array(list(map(lambda sp: sp[1], spec_list)))

# peak intensities
sp_list = list(map(lambda sp: name_to_spectrum(sp[0], s), spec_list))

if None in sp_list:
  none_name = spec_list[sp_list.index(None)][0]
  print(none_name + ' does not exist in your spectrum list.')

ydata = np.array(list(map(lambda sp: sp.data_height(pos), sp_list)))

# fitting
popt, pcov = curve_fit(exp_func, xdata, ydata,
                      bounds=([max(ydata), 0], [10*max(ydata), 10000]))
point_sd = list(np.sqrt(np.diag(pcov)))

# plotting
xlabel = 'Time Parameter (ms)'
ylabel = 'Intensity (I)'
title = 'Curvefit for ' + ref_peak.assignment

plt.figure()
plt.plot(xdata, ydata, '.r:', label='data')
plt.plot(xdata, exp_func(xdata, *popt), 'b-',
     label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend()
plt.pause(0.1)
plt.show(block=False)
