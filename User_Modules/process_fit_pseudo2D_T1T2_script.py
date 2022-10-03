#
# This is an example script for processing and fitting T1/T2 pseudo-2D data.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#     and select a directory that has fid file in.

import nmrglue as ng
import numpy as np

import os
import __main__
s = __main__.session

brukdir = s.open_directorydialog('Select the directory that has ser file.', '')
if brukdir == '':
  raise SystemError

ans = s.show_inputdialog('Considered range', 
  'Set range to consider (e.g. 8.5-13.5)', '8.5-13.5')

try:
  ans_list = ans.split('-')
  range_min, range_max = float(ans_list[0].strip()), float(ans_list[1].strip())
except:
  range_min, range_max = -1000, 1000
  s.show_message('Warning', 'Invalid input. A whole spectrum will be used.')

fidfile = os.path.join(brukdir, 'ser')
vdlistfile = os.path.join(brukdir, 'vdlist')

if not os.path.exists(fidfile):
  s.show_message('Error', 'ser file not found.')  
  raise SystemError

if not os.path.exists(vdlistfile):
  s.show_message('Error', 'vdlist file not found.')  
  raise SystemError

# read vdlist
f = open(vdlistfile, 'r')
vdlines = f.readlines()
f.close()
vd_list = []
for i in range(len(vdlines)):
  if vdlines[i].strip().endswith('ms'):
    vd_list.append(float(vdlines[i].replace('ms','').strip())/1000.)
  elif vdlines[i].strip().endswith('ns'):
    vd_list.append(float(vdlines[i].replace('ns','').strip())/1000000.)  
  elif vdlines[i].strip().endswith('ps'):
    vd_list.append(float(vdlines[i].replace('ps','').strip())/1000000.)  
  elif vdlines[i].strip().endswith('s'): 
    vd_list.append(float(vdlines[i].replace('s','').strip()))
  else:
    vd_list.append(float(vdlines[i].strip()))

# data processing
dic, data = ng.fileio.bruker.read(brukdir)
data = ng.bruker.remove_digital_filter(dic, data)
data = ng.proc_base.zf_auto(data)
data = ng.proc_base.rev(data)
data = ng.proc_base.fft(data)
data = ng.proc_autophase.autops(data, 'acme')
udic = ng.bruker.guess_udic(dic,data)
uc = ng.fileiobase.uc_from_udic(udic)

sum_list = []
xdata = uc.ppm_scale()
for i in range(len(vd_list)):
  idc = np.where(np.logical_and(
      xdata > range_min, xdata < range_max))
  sum_list.append(np.sum(data[i,idc].real))

# calculate T1/T2
def exp_func(x, a, b):
  return a * np.exp(-1 * x / b)
from scipy.optimize import curve_fit

# Calculation function definition
# This uses a whole integration
def calcT(vd_list, sum_list):
  # fitting
  popt, pcov = curve_fit(exp_func, vd_list, sum_list,
          maxfev=1000000000,
          bounds=([max(sum_list), 0], [10*max(sum_list), 10000]))
  point_sd = list(np.sqrt(np.diag(pcov)))
  
  print('T-decay: ' + str(popt[1]))
  print('Deviation: ' + str(point_sd[1]))
  
calcT(vd_list, sum_list)