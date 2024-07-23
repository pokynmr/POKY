#
# This is an example script for applying TSVD denoising on the selected spectrum.
#
# It can take pretty long depending on the data size.
# Make sure you have saved all the data before using this.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Last update: July 23, 2024 
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
from sklearn.decomposition import TruncatedSVD
from sputil import name_to_spectrum
# POKY libraries
import __main__
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

import os.path

denoise_path = s.save_filedialog('Save as', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if denoise_path == '':
  raise SystemError

try:
  n_comp = int(s.show_inputdialog('The number of components',
        f'Type the number of components to retain (e.g. {sp.data_size[-1]-1})',
        str(sp.data_size[-1]-1)))
except:
  s.show_message('Error', 
                 'Please type a number.')
  raise SystemError

dic, data = ng.sparky.read(sp.data_path)

svd = TruncatedSVD(n_components=n_comp)
denoised_data = svd.fit_transform(data)
reconstructed_data = svd.inverse_transform(denoised_data)

ng.sparky.write(denoise_path, dic, reconstructed_data, overwrite=True)

if s.show_message_yes_no('Load data', 
                         'Do you want to load the processed spectra?'):
  s.open_spectrum(denoise_path)
