#
# This is an example script for processing and plotting Bruker 1D data.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#     and select a directory that has fid file in.

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np

import os
import __main__
s = __main__.session

brukdir = s.open_directorydialog('Select the directory that has fid file.', '')
if brukdir == '':
  raise SystemError
fidfile = os.path.join(brukdir, 'fid')
if not os.path.exists(fidfile):
  s.show_message('Error', 'fid file not found.')  
  raise SystemError
dic, data = ng.fileio.bruker.read(brukdir)
data = ng.bruker.remove_digital_filter(dic, data)
data = ng.proc_base.zf_auto(data)
data = ng.proc_base.rev(data)
data = ng.proc_base.fft(data)
data = ng.proc_autophase.autops(data, 'acme')

ft1file = s.save_filedialog('Save as...', 'NMRPipe .ft1 (*.ft1);; Any (*)', brukdir)
if ft1file != '':
  C = ng.convert.converter()
  C.from_bruker(dic, data)
  pdic, pdata = C.to_pipe()
  ng.pipe.write(ft1file, pdic, pdata, True)
          
if s.show_message_yes_no('Plot', 'Do you want to plot the data?'):
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
