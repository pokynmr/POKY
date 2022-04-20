#
# This is an example script for scaling the data height of the 
# selected spectrum in the project.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: April 20, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# POKY 01/14/22e or higher is required to run this script.
#

import numpy as np
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

transform_engine = s.show_selectionexdialog('Engine selection', 
    'Engine for the transformation: ', 
    ('UCSFTOOL', 'NMRGLUE', 'CANCEL'))
if transform_engine == 2:
  raise SystemError
elif transform_engine == 0:
  #### THIS IS UCSFTOOL VERSION
  import ucsftool3
  ut = ucsftool3.ucsfTool()
  ut.ucsf_open(sp.data_path)
  ut.write_transform(new_path, 'mult', factor, 0, overwrite=1)
  ut.ucsf_close()
elif transform_engine == 1:
  #### THIS IS NMRGLUE VERSION
  import nmrglue as ng
  dic, data = ng.sparky.read_lowmem(sp.data_path)
  new_data = np.multiply(data, factor)
  ng.sparky.write_lowmem(new_path, dic, new_data, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(new_path)

