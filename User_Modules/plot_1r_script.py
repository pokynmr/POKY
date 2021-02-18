#
# This is an example script for plotting processed Bruker 1D data.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#     and select a directory that has 1r file in.

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
dic, data = ng.fileio.bruker.read_pdata(brukdir)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
udic = ng.bruker.guess_udic(dic,data)
uc = ng.fileiobase.uc_from_udic(udic)
ax.plot(uc.ppm_scale(), data, 'k-')

# decorate axes
ax.set_xlabel("Chemical Shift (ppm)")
ax.set_ylabel("Amplitude")
ax.set_title("Poky 1D Spectrum")
ax.invert_xaxis()

plt.show(block=False)
