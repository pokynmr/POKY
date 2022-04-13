#
# This is an example script for scaling the data height of the 
# selected spectrum in the project.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: April 13, 2022
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
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

import os.path
new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

try:
  factor = float(s.show_inputdialog('Scaling',
            'Multiply by ', '1000000'))
except:
  s.show_message('Error', 'Value must be numeric.')
  raise SystemError
  
dic, data = ng.sparky.read_lowmem(sp.data_path)
new_data = np.multiply(data, factor)

ng.sparky.write_lowmem(new_path, dic, new_data, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(new_path)

