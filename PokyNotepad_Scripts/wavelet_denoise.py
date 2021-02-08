#
# This is an example script for applying wavelet denoising to the
# selected spectrum in the project.
#
# It can take pretty long depending on the data size.
# Make sure you have saved all the data before using this.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: January 27, 2021
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
from skimage.restoration import denoise_wavelet
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
  factor = float(s.show_inputdialog('Sigma scale',
            'Sigma scale to the estimated noise ' + str(sp.noise), '1'))
except:
  s.show_message('Error', 'Value must be numeric.')
  raise SystemError

proj = s.project
sigma = abs(sp.noise * factor)


dic, data = ng.sparky.read_lowmem(sp.data_path)
data_denoise = denoise_wavelet(data, sigma, method='BayesShrink', mode='soft',
                wavelet_levels=3, wavelet='sym8', rescale_sigma='True')

ng.sparky.write_lowmem(new_path, dic, data_denoise, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(new_path)

