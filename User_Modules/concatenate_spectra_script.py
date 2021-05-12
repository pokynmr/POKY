#
# This is an example script for concatenate two spectra.
# Two spectra should share the same ni, spectral width, etc.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: May 12, 2021
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
from sputil import name_to_spectrum
# POKY libraries
import __main__
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select spectrum one', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError


specname2 = s.show_spectrumselectiondialog('Select spectrum two', 0)
sp2 = name_to_spectrum(specname2, s)

if sp2 == None:
  raise SystemError

import os.path
new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

try:
  factor = float(s.show_inputdialog('Scale factor',
            'Multiply by ...', '1'))
except:
  s.show_message('Error', 'Value must be numeric.')
  raise SystemError

proj = s.project

dic, data = ng.sparky.read_lowmem(sp.data_path)
dic2, data2 = ng.sparky.read_lowmem(sp2.data_path)

data3 = np.array(data) + factor*np.array(data2)

ng.sparky.write_lowmem(new_path, dic, data3, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the concatenated spectrum?'):
  s.open_spectrum(new_path)

