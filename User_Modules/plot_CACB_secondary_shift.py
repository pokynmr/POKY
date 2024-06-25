#
# This is an example script to plot dCA-dCB seconday chemical shifts.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Secondary CA/CB Shift Plot')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# Reference to J. Biomol. NMR 1995.5 67-81
refCSdict = {# CA    CB    C       H      N       AA
      'A':(  52.5,  19.1,  177.8,  8.24,  123.8),  # Ala
      'Cr':( 58.2,  28.0,  174.6,  8.32,  118.8),  # Cys_r
      'Co':( 55.4,  41.1,  174.6,  8.43,  118.6),  # Cys_o
      'D':(  54.2,  41.1,  176.3,  8.34,  120.4),  # Asp
      'E':(  56.6,  29.9,  176.6,  8.42,  120.2),  # Glu
      'F':(  57.7,  39.6,  175.8,  8.30,  120.3),  # Phe
      'G':(  45.1,  0.00,  174.9,  8.33,  108.8),  # Gly
      'H':(  55.0,  29.0,  174.1,  8.42,  118.2),  # His
      'I':(  61.1,  38.8,  176.4,  8.00,  119.9),  # Ile
      'K':(  56.2,  33.1,  176.6,  8.29,  120.4),  # Lys
      'L':(  55.1,  42.4,  177.6,  8.16,  121.8),  # Leu
      'M':(  55.4,  32.9,  176.3,  8.28,  119.6),  # Met
      'N':(  53.1,  38.9,  175.2,  8.40,  118.7),  # Asn
      'P':(  63.3,  32.1,  177.3,  0.00,   0.00),  # Pro
      'Q':(  55.7,  29.4,  176.0,  8.32,  119.8),  # Gln
      'R':(  56.0,  30.9,  176.3,  8.23,  120.5),  # Arg
      'S':(  58.3,  63.8,  174.6,  8.31,  115.7),  # Ser
      'T':(  61.8,  69.8,  174.7,  8.15,  113.6),  # Thr
      'V':(  62.2,  32.9,  176.3,  8.03,  119.2),  # Val
      'W':(  57.5,  29.6,  176.1,  8.25,  121.3),  # Trp
      'Y':(  57.9,  38.8,  175.9,  8.12,  120.3)   # Tyr
    }

import __main__
s = __main__.main_session
clist = s.project.condition_list()
if len(clist) > 1:
  cname = s.show_conditionselectiondialog('Select a condition to evaluate.', 0)
  from sputil import name_to_condition
  c = name_to_condition(cname, s)
else:
  c = clist[0]

x_list_h, y_list_h = [], []
x_list_e, y_list_e = [], []

data_dict = {}
nmin, nmax = 10**9, -10**9
for resn in c.resonance_list():
  a = resn.group.symbol
  n = resn.group.number
  cs = resn.frequency
  try:
    idx = ['CA', 'CB'].index(resn.atom.name)
  except:
    continue
  if cs == 0.0:
    continue
  if not str(n) in data_dict and not str(n)+'r' in data_dict:
    if a != 'C':
      data_dict[str(n)] = [0, 0]
    else:
      data_dict[str(n)+'r'] = [0, 0]
      data_dict[str(n)+'o'] = [0, 0]
  try:
    if a != 'C':
      row = refCSdict[resn.group.symbol]
      data_dict[str(n)][idx] = cs - row[idx]
    else:
      data_dict[str(n)+'r'][idx] = cs - refCSdict['Cr'][idx]
      data_dict[str(n)+'o'][idx] = cs - refCSdict['Co'][idx]
  except:
    continue
  nmin = min(nmin, n)
  nmax = max(nmax, n)

for n in range(nmin, nmax + 1):
  if str(n) in data_dict:
    diff = data_dict[str(n)][0] - data_dict[str(n)][1]
  elif str(n)+'r' in data_dict:
    diffr = data_dict[str(n)+'r'][0] - data_dict[str(n)+'r'][1]
    diffo = data_dict[str(n)+'o'][0] - data_dict[str(n)+'o'][1]
    diff = min(diffr, diffo)
  else:
    continue
  if diff > 0:
    x_list_h.append(n)
    y_list_h.append(diff)
  else:
    x_list_e.append(n)
    y_list_e.append(diff)

# plotting
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
xlabel = 'Residue Number'
ylabel = 'd(dCA-dCB) (ppm)'
plt.figure()

plt.bar(x_list_h, y_list_h, color = 'green')
plt.bar(x_list_e, y_list_e, color = 'blue')

plt.xlabel(xlabel)
plt.ylabel(ylabel)
ylim = max(np.abs(y_list_h+y_list_e)) * 1.1
plt.xlim((nmin, nmax))
plt.ylim((-1 * ylim, ylim))
plt.title('Poky Secondary CA/CB Shift Plot')
plt.pause(0.1)
plt.show(block=False)
