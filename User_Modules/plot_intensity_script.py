#
# This is an example script to plot sorted peak intensities.
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
print('POKY Plot Peak Intensities')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

spec_name = s.show_spectrumselectiondialog('Select a spectrum', 0)

from sputil import name_to_spectrum, sort_peaks_by_height
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

spec = name_to_spectrum(spec_name, s)
if spec == None:
  raise SystemExit

peaks = sort_peaks_by_height(spec.peak_list())
xdata = np.array(list(range(1, len(peaks)+1)))
ydata = np.array(list(map(lambda p: p.data_height, peaks)))
lnydata = np.array(list(map(lambda p: np.log(p.data_height), peaks)))

# plot
xlabel = 'Peaks'
ylabel = 'Intensities (I)'
lnylabel = 'Logarithmic Intensities (ln I)'
title = 'Poky Intensity Plot'
lntitle = 'Poky Logarithmic Intensity Plot'

plt.figure()
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.plot(xdata, ydata)
plt.title(title)
plt.pause(0.1)

plt.figure()
plt.xlabel(xlabel)
plt.ylabel(lnylabel)
plt.plot(xdata, lnydata)
plt.title(title)
plt.pause(0.1)

plt.show(block=False)
