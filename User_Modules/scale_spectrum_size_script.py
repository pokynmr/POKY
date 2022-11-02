#
# This is an example script for scaling the spectrum size.
# Spectra will be chosen from the running Poky session.
#
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Nov. 1, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
import scipy.ndimage
from sputil import name_to_spectrum

# POKY libraries
import __main__
s = __main__.main_session

specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

orig_size_list = list(map(str, sp.data_size))
orig_size = ','.join(orig_size_list)
new_size_list = list(map(lambda x: str(int(x / 2)), sp.data_size))
new_size = ','.join(new_size_list)
ans = s.show_inputdialog('New size', 'Current size: ' + orig_size + \
                  '\n\nNew size: ', new_size)
new_size_list = list(map(lambda x: int(x.strip()), ans.split(',')))
if len(new_size_list) != len(sp.data_size):
  raise SystemError

import os.path

new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

dic, data = ng.sparky.read(sp.data_path)
zoom_list = []
for i in range(sp.dimension):
  dic['w' + str(i+1)]['npoints'] = new_size_list[i]
  zoom_list.append(new_size_list[i] / sp.data_size[i])

new_data = scipy.ndimage.zoom(data, zoom_list, mode='nearest')

ng.sparky.write(new_path, dic, new_data, overwrite=True)

if s.show_message_yes_no('Load data', 
              'Do you want to load the scaled spectrum?'):
  s.open_spectrum(new_path)
