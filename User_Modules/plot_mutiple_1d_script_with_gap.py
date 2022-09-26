#
# This is an example script to plot 1D slice of selected multiple peaks.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#
#   1. Select multiple peaks from different spectra. 
#      Hold SHIFT key while selecting a peak not to lose a peak selection.
#
#   2. In Poky Notepad,
#     File -> Run Python Module
#
#   To distinguish peaks by elevating intensity, modify the intensity_gap.
#

# intensity gap between peaks
intensity_gap = 0 #10**5

import __main__
s = __main__.main_session
proj = s.project

if len(s.selected_peaks()) < 1:
  print('Peak should be selected.')
  raise SystemExit

peaks = s.selected_peaks()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')
cmap = 'bgrcmyb'
mmap = ['-', '--', '-.', ':', '.',',', 'o', 'v', '^', '<', '>', '1', '2', 
  '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']

for i in range(len(peaks)):
  peak = peaks[i]
  pos = peak.position
  spec = peak.spectrum

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
    
    ydata += intensity_gap * i
    # plotting
    if i == 0:
      proj_pos = list(pos)
      proj_pos.pop(dim)
      xlabel = nucleus + ' (ppm)'
      ylabel = 'Data Height'
      title = 'Poky 1D Slice at'
      for j in range(len(proj_pos)):
        title += ' %.2f ppm' % (proj_pos[j])
      if peak.is_assigned:
        title += ' for ' + peak.resonances()[0].group.name
      
    plt.figure(str(dim+1)+ ' dimension slice')
    if i == 0:
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)
      plt.title(title)
    try:
      from colormap import get_contour_hex_color
      c = get_contour_hex_color(peak.spectrum)
      plt.plot(xdata, ydata, color=c, linestyle='-')
    except:
      plt.plot(xdata, ydata, cmap[i % len(cmap)] + mmap[i // len(cmap) % len(mmap)])
    ax = plt.gca()
    ax.invert_xaxis()

plt.pause(0.1)
plt.show(block=False)