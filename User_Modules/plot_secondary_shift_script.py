#
# This is an example script to plot seconday chemical shifts.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY Secondary Shift Plot')
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

btn_list = ('CA', 'CB', 'C', 'H', 'N', 'Cancel')
idx = s.show_selectionexdialog('Atom type', 'Select an atom type.', btn_list)
if idx in [-1, len(btn_list)-1]:
  raise SystemExit

x_list, y_list, err_list = [], [], []
for resn in c.resonance_list():
  a = resn.group.symbol
  cs = resn.frequency
  if cs == 0.0 or resn.atom.name != btn_list[idx]:
    continue
  try:
    if a != 'C':
      row = refCSdict[resn.group.symbol]
      if row[idx] == 0.0:
        continue
      dcs = cs - row[idx]
    else: # CYS will be evaluated to the closer value from r. vs o.
      row = refCSdict[resn.group.symbol + 'r']
      dcs = cs - row[idx]
      row2 = refCSdict[resn.group.symbol + 'o']
      dcs2 = cs - row2[idx]
      if abs(dcs2) < abs(dcs): 
        dcs = dcs2
  except:
    continue

  x_list.append(resn.group.number)
  y_list.append(dcs)
  err_list.append(resn.deviation)
  
# plotting
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

import matplotlib.pyplot as plt
xlabel = 'Residue Number'
ylabel = 'd(CSobs - CSref)'
plt.figure()

plt.errorbar(x_list, y_list, err_list, fmt='bo', markersize=5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
ylim = max(max(y_list), abs(min(y_list))) * 1.1
plt.ylim((-1 * ylim, ylim))
plt.title('Poky Secondary Shift Plot: ' + btn_list[idx])
plt.pause(0.1)
plt.show(block=False)
