#
# This is an example script to conduct ChemEx analysis.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

# This needs more development

import __main__
s = __main__.main_session
proj = s.project


if len(s.selected_peaks()) != 1:
  print('One peak should be selected.')
  raise SystemExit

# If more than two spectra are given, automatically detect the lowest dip
# However, the first spectrum is the saturated one giving I0
spec_list = [ # place spectrum names and (time) parameters. e.g.
            ['relaxT_10_1_N15_T2', 10],
            ['relaxT_30_2_N15_T2', 30],
            ['relaxT_50_3_N15_T2', 50],
            ['relaxT_70_4_N15_T2', 70],
            ['relaxT_90_5_N15_T2', 90],
            ['relaxT_110_6_N15_T2', 110],
            ]

import chemex


# default nucleus is 15N
nucleus = '15N'

from sputil import name_to_spectrum
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

sp_list = list(map(lambda sp: name_to_spectrum(sp[0], s), spec_list))

# set reference peak and position
ref_peak = s.selected_peaks()[0]
pos = ref_peak.position

ref_spec = sp_list[0]
ref_noise = ref_spec.noise
region = ref_spec.region
dim = ref_spec.nuclei.index(nucleus)
npoint = ref_spec.data_size[dim]

# find the lowest signal from other spectra
heights = list(map(lambda sp: sp.data_height(ref_peak.position),
                                sp_list[1:]))
dip_spec = sp_list[heights.index(min(heights)) + 1][0]
dip_noise - dip_spec.noise

ppm_per_pt = ref_spec.spectrum_width[dim] / (npoint - 1)
xdata = np.array(list(map(lambda x: region[1][dim] - ppm_per_pt * x,
                                    range(npoint))))
if dim == 0:
  I0data = np.array(list(map(lambda x:
                  ref_spec.data_height((xdata[x], pos[1])), range(npoint))))
  Idata = np.array(list(map(lambda x:
                  dip_spec.data_height((xdata[x], pos[1])), range(npoint))))
else:
  I0data = np.array(list(map(lambda x:
                  ref_spec.data_height((pos[0], xdata[x])), range(npoint))))
  Idata = np.array(list(map(lambda x:
                  dip_spec.data_height((pos[0], xdata[x])), range(npoint))))
ref_scale = sum(
ydata = np.array(list(map(lambda x:max(0.0, min(1.0, Idata[x] / I0data[x])), range(npoint))))


#scale = sum(
#             values * self.data["intensity"] / self.data["error"] ** 2
#         ) / sum((values / self.data["error"]) ** 2)

# plotting
xlabel = nucleus + ' (ppm)'
ylabel = 'Intensity Ratio (I/I0)'
title = nucleus + ' CEST profile'
if ref_peak.is_assigned:
  title += ' for ' + ref_peak.resonances()[0].group.name

plt.figure()
plt.plot(xdata, ydata, 'bo', markersize=5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.pause(0.1)
plt.show(block=False)
