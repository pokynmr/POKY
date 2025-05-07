#
# This is an example script to plot/print 1D slice of a selected 2D peak.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session
peaks = s.selected_peaks()

if len(peaks) != 1:
  s.show_message('Error', 'Please select one 2D peak.')
  raise SystemError

peak = peaks[0]
spec = peak.spectrum
if spec.dimension != 2:
  s.show_message('Error', 'Please select one 2D peak.')
  raise SystemError

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')


dim, dim2 = 0, 1
region = spec.region

# dim 1
nucleus = spec.nuclei[dim]
npoint = spec.data_size[dim]
ppm_per_pt = spec.spectrum_width[dim] / npoint

# dim 2
nucleus2 = spec.nuclei[dim2]
npoint2 = spec.data_size[dim2]
ppm_per_pt2 = spec.spectrum_width[dim2] / npoint2

xdata = np.array(list(map(lambda x: region[1][dim] - ppm_per_pt * x,
                                  range(npoint))))
ydata = np.array(list(map(lambda x: region[1][dim2] - ppm_per_pt2 * x,
                                  range(npoint2))))

# intensity
idata = [0.,] * npoint
idata2 = [0.,] * npoint2
idata = list(map(lambda x: spec.data_height((xdata[x], peak.frequency[1])), 
                 range(npoint)))
idata2 = list(map(lambda x: spec.data_height((peak.frequency[0], ydata[x])), 
                 range(npoint2)))
    
# plotting
def plot(nucleus, xdata, ydata, bt):
  xlabel = nucleus + ' (ppm)'
  ylabel = 'Data Height'
  title = 'POKY 1D Slice ' + bt
    
  plt.figure()
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  try:
    from colormap import get_contour_hex_color
    c = get_contour_hex_color(spec)
    plt.plot(xdata, ydata, color=c, linestyle='-')
  except:
    plt.plot(xdata, ydata, 'b-')
  for i, x in enumerate(xdata):
    y = ydata[i]
    print('%8.3f %8.3f' % (x, y))
  
  ax = plt.gca()
  ax.invert_xaxis()

  plt.pause(0.1)
  plt.show(block=False)

print('\n\n------------------------------')
print('POKY 1D Slice At %.3f ppm' % (peak.frequency[1]))
print('------------------------------\n\n')
plot(nucleus, xdata, idata, 'At %.3f ppm' % (peak.frequency[1]))

print('\n\n------------------------------')
print('POKY 1D Slice At %.3f ppm' % (peak.frequency[0]))
print('------------------------------\n\n')
plot(nucleus2, ydata, idata2, 'At %.3f ppm' % (peak.frequency[0]))
