#
# This is an example script for scaling the spectral width.
# Spectra will be chosen from the running Poky session.
#
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jun. 19, 2024
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import nmrglue as ng
from sputil import name_to_spectrum

# POKY libraries
import __main__
s = __main__.main_session

specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

w_list = list(map(lambda x: f'w{x+1}', range(sp.dimension)))
ndim = s.show_selectionexdialog('Scale spectral width', 
                         'Choose the dimension to scale', w_list)
if ndim == -1:
  raise SystemError

try:
  ans = float(s.show_inputdialog('Scale', 'Scale factor: ', '0.5'))
except:
  print('Wrong scale factor.')

import os.path

new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

dic, data = ng.sparky.read(sp.data_path)
dic['w' + str(ndim+1)]['spectral_width'] *= ans
ng.sparky.write(new_path, dic, data, overwrite=True)

if s.show_message_yes_no('Load data', 
              'Do you want to load the scaled spectrum?'):
  s.open_spectrum(new_path)
