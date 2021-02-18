#
# This is an example script for plotting processed NMRPipe 1D data.
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

pipefile = s.open_filedialog('Select the 1D NMRPipe file.', 
              'NMRPipe .ft1 (*.ft1);; Any (*)', '')
if pipefile == '':
  raise SystemError
dic, data = ng.pipe.read(pipefile)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
uc = ng.pipe.make_uc(dic,data)
ax.plot(uc.ppm_scale(), data, 'k-')

# decorate axes
ax.set_xlabel("Chemical Shift (ppm)")
ax.set_ylabel("Amplitude")
ax.set_title("Poky 1D Spectrum")
ax.invert_xaxis()

plt.show(block=False)
