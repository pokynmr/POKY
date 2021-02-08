#
# This is an example script to create a bar chart from selected spectrum.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Peak Intensity Bar Chart')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# Set to False below if all peaks to be considered regardless of assignment
only_assigned = True

# Set to True below if only selected peaks to be considered
only_selected = False

import __main__
s = __main__.main_session
proj = s.project
vname = s.show_viewselectiondialog('Select a view', 0)

from sputil import name_to_view
v = name_to_view(vname, s)

if v == None:
  print('Select a view first.')
  sys.exit()

sp = v.spectrum

from sputil import sort_peaks_by_assignment
peaks = list(map(lambda p: p, sp.peak_list()))

if only_assigned:
  temp_peaks = []
  for p in peaks:
    if only_selected and p.selected == 0:
      continue
    if None not in p.resonances():
      temp_peaks.append(p)
  peaks = sort_peaks_by_assignment(temp_peaks, False)

x, y = [], []

if only_assigned:
  for p in peaks:
    if p.is_assigned == 0:
      continue
    x.append(p.resonances()[0].group.number)
    y.append(p.data_height)
else:
  for p in peaks:
    if only_selected and p.selected == 0:
      continue
    x.append(len(x)+1)
    y.append(p.data_height)

from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.bar(x, y)
ax.tick_params(axis='both', labelbottom=True, bottom=True,
          top=False, labelleft=True, left=True, right=False)
ax.set_ylabel('Data Height')
if only_assigned:
  ax.set_xlabel('Residue Number')
else:
  ax.set_xlabel('Arbitrary Number')
plt.pause(0.1)
plt.show(block=False)
