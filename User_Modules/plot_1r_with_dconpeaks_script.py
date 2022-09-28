#
# This is an example script for plotting processed Bruker 1D data with fit.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#     and select a directory that has 1r and dconpeaks.txt file in.
#     Then, you can save .csv file with the data

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

import os
import __main__
s = __main__.session

brukdir = s.open_directorydialog('Select the directory that has 1r file.', '')
if brukdir == '':
  raise SystemError
pfile = os.path.join(brukdir, '1r')
if not os.path.exists(pfile):
  s.show_message('Error', '1r file not found.')  
  raise SystemError
dic, data = ng.fileio.bruker.read_pdata(brukdir, scale_data=True)

# Read dconpeaks.txt if exists
dcfile = os.path.join(brukdir, 'dconpeaks.txt')
if not os.path.exists(dcfile):
  s.show_message('Error', 'dconpeaks.txt file not found.')
  raise SystemError
f = open(dcfile, 'r')
dclines = f.readlines()
f.close()
fit_type_list = ['Mixed Lorentzian and Gaussian', 
                  'Lorentzian', 'Gaussian']
for line in dclines:
  if line.find('Fit type:') != -1:
    fit_type = fit_type_list.index(line.split(':')[1].strip())
    break

# gaussian
def gaussian(x, amp, freq, lw):
  return amp * np.exp(-4*np.log(2)*(x - freq) ** 2. / lw)

# lorentzian
def lorentzian(x, amp, freq, lw):
  return amp * (lw/2.) ** 2. / ((x - freq) ** 2. + (lw/2.) ** 2.)

# plot the spectrum
plt.figure()

# decorate axes
plt.xlabel("Chemical Shift (ppm)")
plt.ylabel("Amplitude")
plt.title("Poky 1D Spectrum")
plt.gca().invert_xaxis()

udic = ng.bruker.guess_udic(dic, data)
uc = ng.fileiobase.uc_from_udic(udic)
xdata = uc.ppm_scale()
vdata = np.array(xdata)
vdata = np.vstack((vdata, data))

plt.plot(xdata, data, 'k-')
# draw fits
for i in range(6, len(dclines)-2):
  sp_list = dclines[i].split()
  if len(sp_list) == 0:
    continue
  try:
    nfit = int(sp_list[0])
  except:
    continue
  sp_list = dclines[i+2].split()
  freq, lw, amp = float(sp_list[0]), float(sp_list[2]), float(sp_list[4])
  gdata = gaussian(xdata, amp, freq, lw)
  ldata = lorentzian(xdata, amp, freq, lw)
  if fit_type == 0: # Pseudo-voift
    lmix = float(sp_list[6]) / 100.
    gmix = 1. - lmix
    ydata = ldata*lmix + gdata*gmix
  elif fit_type == 1: # Lorentzian
    ydata = ldata
  elif fit_type == 1: # Gaussian
    ydata = gdata
  vdata = np.vstack((vdata, ydata))
  
  plt.plot(xdata, ydata, 'r-')  
  
plt.show(block=False)

fname = s.save_filedialog('Save as', 'CSV (*.csv);; Any (*)', '')

if fname != '':
  text = 'PPM,DATA'
  for i in range(len(vdata[2:])):
    text += ',PEAK'+str(i+1)
  text += '\n'
  for i in range(len(vdata[0])):
    for j in range(len(vdata)):
      text += str(vdata[j, i]) + ','
    text = text[:-1]
    text += '\n'
  f = open(fname, 'w')
  f.write(text)
  f.close()
