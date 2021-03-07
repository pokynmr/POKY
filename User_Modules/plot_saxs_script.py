#
# This is an example script to plot SAXS data
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

input_file = s.open_filedialog('Select a SAXS file', 
                    'SAXS (*.dat);; Any (*)', '')
if input_file == '':
  raise SystemExit

f = open(input_file, 'r')
lines = f.readlines()
f.close()

xdata, ydata, edata = [], [], []

for line in lines:
  sp = line.strip().split()
  if len(sp) < 3:
    continue
  if sp[0][0] == '#':
    continue
  try:
    x, y, err = float(sp[0]), float(sp[1]), float(sp[2])
    xdata += [x,]
    ydata += [y,]
    edata += [err,]
  except:
    continue

import matplotlib.pyplot as plt

title = 'POKY SAXS PLOT'
xlabel = 'q(A^-1)'
ylabel = 'I(q)'
plt.figure()
plt.errorbar(xdata, ydata, yerr=edata, fmt='bo', markersize=3)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.pause(0.1)
plt.show(block=False)

# Kratky plot
title = 'POKY SAXS KRATKY PLOT'
xlabel = 'q(A^-1)'
ylabel = 'q*q*I(q)'
plt.figure()
ydata2 = list(map(lambda i: xdata[i]**2*ydata[i], range(len(ydata))))
plt.scatter(xdata, ydata2, c='red')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.pause(0.1)
plt.show(block=False)

# Guinier plot
import numpy as np
from scipy.optimize import curve_fit

title = 'POKY SAXS GUINIER PLOT'
xlabel = 'q*q(A^-2)'
ylabel = 'ln I(q)'
plt.figure()
xdata3 = list(map(lambda x: x*x, xdata))
ydata3 = list(map(lambda y: np.log(y), ydata))
plt.scatter(xdata3, ydata3, c='green')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.pause(0.1)
plt.show(block=False)

def guinier_func(x, a, b):
  return a - 1/3 * x * b

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

ans = s.show_inputdialog('Approximation', 'Type q*q range to approximate',
                    '0,0.006')

min_qq, max_qq = ans.split(',')
min_qq, max_qq = float(min_qq), float(max_qq)

idx_min = find_nearest(xdata3, min_qq)
idx_max = find_nearest(xdata3, max_qq)

popt, pcov = curve_fit(guinier_func, 
              xdata3[idx_min:idx_max+1], ydata3[idx_min:idx_max+1])
print(popt)
print('ln(I0): ' + str(popt[0]))
print('I0: ' + str(np.exp(popt[0])))
print('Rg^2: ' + str(popt[1]))
print('Rg: ' + str(popt[1]*0.5))
