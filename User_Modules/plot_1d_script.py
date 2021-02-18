#
# This is an example script to plot 1D slice of a selected peak.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project


if len(s.selected_peaks()) != 1:
  print('One peak should be selected.')
  raise SystemExit

peak = s.selected_peaks()[0]
pos = peak.position
spec = peak.spectrum

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

for dim in range(len(pos)):
  nucleus = spec.nuclei[dim]
  region, npoint = spec.region, spec.data_size[dim]
  ppm_per_pt = spec.spectrum_width[dim] / (npoint - 1)
  xdata = np.array(list(map(lambda x: region[1][dim] - ppm_per_pt * x,
                                    range(npoint))))
  if len(pos) == 2:
    if dim == 0:
      ydata = np.array(list(map(lambda x:
                  spec.data_height((xdata[x], pos[1])), range(npoint))))
    else:
      ydata = np.array(list(map(lambda x:
                  spec.data_height((pos[0], xdata[x])), range(npoint))))

  elif len(pos) == 3:
    if dim == 0:
      ydata = np.array(list(map(lambda x:
                  spec.data_height((xdata[x], pos[1], pos[2])), range(npoint))))
    elif dim == 1:
      ydata = np.array(list(map(lambda x:
                  spec.data_height((pos[0], xdata[x], pos[2])), range(npoint))))
    else:
      ydata = np.array(list(map(lambda x:
                  spec.data_height((pos[0], pos[1], xdata[x])), range(npoint))))

  # plotting
  proj_pos = list(pos)
  proj_pos.pop(dim)
  xlabel = nucleus + ' (ppm)'
  ylabel = 'Data Height'
  title = 'Poky 1D Slice at'
  for i in range(len(proj_pos)):
    title += ' %.2f ppm' % (proj_pos[i])
  if peak.is_assigned:
    title += ' for ' + peak.resonances()[0].group.name

  plt.figure()
  plt.plot(xdata, ydata, 'b-')
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.pause(0.1)
  plt.show(block=False)
