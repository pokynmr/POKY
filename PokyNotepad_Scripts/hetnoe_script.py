#
# This is an example script to create a hetNOE chart.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Heteronuclear NOE Chart')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project

# Parameters
ref_spec_name = s.show_spectrumselectiondialog('Reference spectrum', 0).strip()
noe_spec_name = s.show_spectrumselectiondialog('NOE spectrum', 0).strip()
only_selected = False # Use all assigned peaks if False
error_bar = True # Draw error bar if True

# Processing start
from sputil import name_to_spectrum, sort_peaks_by_assignment

ref_spec = name_to_spectrum(ref_spec_name, s)
noe_spec = name_to_spectrum(noe_spec_name, s)

if ref_spec == None:
  print('Reference spectrum name %s is wrong.' % ref_spec_name)
  raise SystemExit
if noe_spec == None:
  print('NOE spectrum name %s is wrong.' % ref_spec_name)
  raise SystemExit

# collect peaks and sort by number
temp_peaks = []
for p in ref_spec.peak_list():
  if only_selected and p.selected == 0:
    continue
  if None not in p.resonances() and p.is_assigned == 1:
    temp_peaks.append(p)
peaks = sort_peaks_by_assignment(temp_peaks, False)

# making x, y arrays
import numpy as np
xdata = np.array(list(map(lambda x: x.resonances()[0].group.number,
                                        peaks)))
ref_height = np.array(list(map(lambda x: x.data_height, peaks)))
noe_height = np.array(list(map(lambda x: noe_spec.data_height(x.position),
                                        peaks)))
ydata = np.array(list(map(lambda x: noe_height[x] / ref_height[x],
                                        range(len(peaks)))))
# error estimation
## pick 30 random position
##if error_bar:
#ref_rand, noe_rand = [], []
#for i in range(30):
#  randx=np.random.randint(0, ref_spec.data_size[0])
#  randy=np.random.randint(0, ref_spec.data_size[1])
#  ref_rand.append(ref_spec.data_height((randx, randy)))
#  noe_rand.append(noe_spec.data_height((randx, randy)))
rho_ref = ref_spec.noise
rho_noe = noe_spec.noise
# https://pubs.acs.org/doi/pdf/10.1021/bi00185a040
edata = np.array(list(map(lambda x: ((rho_ref/ref_height[x])**2 +
        (rho_noe/noe_height[x])**2)**.5, range(len(xdata)))))

# plotting
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

import matplotlib.pyplot as plt
xlabel = 'Residue Number'
ylabel = '{1H}-15N NOE'
plt.figure()
if error_bar:
  plt.errorbar(xdata, ydata, edata, fmt='bo', markersize=5)
else:
  plt.plot(xdata, ydata, 'bo', markersize=5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.ylim((0,1))
plt.title('Poky hetNOE Chart')
plt.pause(0.1)
plt.show(block=False)

# print out
print('%-12s %-12s %-12s %-12s %-12s' % ('Assignment', 'I(0)', 'I', 'I/I(0)', 'Err'))
for i in range(len(xdata)):
  print('%-12s %-12.3f %-12.3f %-12.3f %-12.3f' %
    (peaks[i].assignment, ref_height[i], noe_height[i], ydata[i], edata[i]))
