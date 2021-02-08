#
# This is an example script to plot talos ss or s2 outputs.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

# only select predS2.tab or predSS.tab
#

import __main__
s = __main__.main_session

tabfile = s.open_filedialog('Select a TALOS-N output file.',
                            'TALOS-N file (*.tab);; Any (*)', '')

if tabfile == '':
  raise SystemExit

mode = ''
if tabfile.find('S2.tab') != -1:
  mode = 'S2'
elif tabfile.find('SS.tab') != -1:
  mode = 'SS'
else:
  print('Please choose S2 or SS output file.')
  raise SystemExit

x, y, y2 = [], [], []

f = open(tabfile, 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
  line = lines[i]
  if line.find('VARS') == 0:
    break

for j in range(i, len(lines)):
  sp = lines[j].strip().split()
  if len(sp) < 5:
    continue
  try:
    x.append(int(sp[0]))
    if mode == 'S2':
      y.append(float(sp[-1]))
    elif mode == 'SS':
      y.append(float(sp[4]))
      y2.append(-1. * float(sp[5]))
      print(lines[j])
  except:
    pass

# plotting
import matplotlib.pyplot as plt
xlabel = 'Residue Number'
if mode == 'S2':
  title = 'Protein Flexibility by RCI-S2'
  ylabel = 'RCI-S2 Order Parameter'
elif mode == 'SS':
  title = 'Protein Secondary Structure Prediction by TALOS-N'
  ylabel = 'Propensity'

plt.figure()
barlist = plt.bar(x, y)
plt.title(title)

for i in range(len(x)):
  if mode == 'S2':
    g = y[i]**2
    r = 1.0 - y[i]**2
    b = 0.0
  elif mode == 'SS':
    r = 1 - y[i]
    g = 1.0
    b = 1 - y[i]
  barlist[i].set_color((r,g,b))
if mode == 'SS':
  barlist = plt.bar(x, y2)
  for i in range(len(x)):
    r = 1 + y2[i]
    b = 1.0
    g = 1 + y2[i]
    barlist[i].set_color((r,g,b))

if mode == 'S2':
  plt.xlim(x[0], x[-1])
  plt.ylim(0, 1.0)
elif mode == 'SS':
  plt.xlim(x[0], x[-1])
  plt.ylim(-1.0, 1.0)

plt.pause(0.1)
plt.show(block=False)
