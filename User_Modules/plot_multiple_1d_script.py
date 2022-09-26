#
# This is an example script to plot multiple 1D slices of a selected peak.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
#
# Runs BUILD 08/25/2022c or newer
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
peak_spec = peak.spectrum

specnames = s.show_spectrumselectiondialog('Select spectra', 1)
specname_list = specnames.split('\t')
if len(specname_list) == 0:
  raise SystemExit

from sputil import name_to_spectrum
sp_list = list(map(lambda sp: name_to_spectrum(sp, s), specname_list))

if None in sp_list:
  none_name = specname_list[sp_list.index(None)]
  print('Spectrum (' + none_name + ') does not exist in your spectrum list.')
  raise SystemExit

from colormap import get_contour_hex_color
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

for dim in range(len(pos)):
  # plotting parameters
  proj_pos = list(pos)
  proj_pos.pop(dim)
  xlabel = peak_spec.nuclei[dim] + ' (ppm)'
  ylabel = 'Data Height'
  title = 'Poky 1D Slice at'
  for i in range(len(proj_pos)):
    title += ' %.2f ppm' % (proj_pos[i])
  if peak.is_assigned:
    title += ' for ' + peak.resonances()[0].group.name
  plt.figure()

  for spec in sp_list:
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
    c = get_contour_hex_color(spec)
    plt.plot(xdata, ydata, color=c, linestyle='-')

  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  ax = plt.gca()
  ax.invert_xaxis()

  plt.pause(0.1)
  plt.show(block=False)