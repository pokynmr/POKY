#
# This is an example script to plot 1D projection of current 2D view window.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session
view = s.selected_view()

if view == None:
  print('Cannot find a view.')
  raise SystemError

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

spec = view.spectrum
if spec.dimension != 2:
  print('Only supports 2D data now.')
  raise SystemError

dim, dim2 = 0, 1
vregion, region = view.region, spec.region

# dim 1
nucleus = spec.nuclei[dim]
npoint = spec.data_size[dim]
vnpoint = int((vregion[1][dim] - vregion[0][dim]) / \
          (region[1][dim] - region[0][dim]) * npoint + 0.5)
ppm_per_pt = spec.spectrum_width[dim] / npoint

# dim 2
nucleus2 = spec.nuclei[dim2]
npoint2 = spec.data_size[dim2]
vnpoint2 = int((vregion[1][dim2] - vregion[0][dim2]) / \
          (region[1][dim2] - region[0][dim2]) * npoint2 + 0.5)
ppm_per_pt2 = spec.spectrum_width[dim2] / npoint2

xdata = np.array(list(map(lambda x: vregion[1][dim] - ppm_per_pt * x,
                                  range(vnpoint))))
ydata = np.array(list(map(lambda x: vregion[1][dim2] - ppm_per_pt2 * x,
                                  range(vnpoint2))))

# intensity
idata2d = np.zeros((vnpoint, vnpoint2))
for i in range(vnpoint):
  for j in range(vnpoint2):
    idata2d[i, j] = spec.data_height((xdata[i], ydata[j]))

# sum along axis
idata = np.sum(idata2d, axis=1)
idata2 = np.sum(idata2d, axis=0)

# plotting
def plot(nucleus, xdata, ydata, bt):
  xlabel = nucleus + ' (ppm)'
  ylabel = 'Data Height'
  title = 'Poky 1D Projection between ' + bt
    
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
  ax = plt.gca()
  ax.invert_xaxis()

  plt.pause(0.1)
  plt.show(block=False)

plot(nucleus, xdata, idata, '%.3f and %.3f' % (vregion[0][1], vregion[1][1]))
plot(nucleus2, ydata, idata2, '%.3f and %.3f' % (vregion[0][0], vregion[1][0]))
