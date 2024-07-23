#
# This is an example script for applying PCA denoising on the selected spectrum.
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
from sklearn.decomposition import PCA
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
  prop_variance = float(s.show_inputdialog('A proportion of variance',
                      'Specify a proportion of variance to retain (e.g. 0.99)',
                      '0.99'))
except:
  s.show_message('Error', 
                 'Please type a number between 0 and 1.')
  raise SystemError

if prop_variance <=0 and prop_variance > 1:
  s.show_message('Error', 
                 'Please type a number between 0 and 1.')
  raise SystemError

dic, data = ng.sparky.read(sp.data_path)

pca = PCA(n_components=prop_variance)
pca_data = pca.fit_transform(data)
reconstruct_data = pca.inverse_transform(pca_data)  

ng.sparky.write(denoise_path, dic, reconstruct_data, overwrite=True)

if s.show_message_yes_no('Load data', 
                         'Do you want to load the processed spectra?'):
  s.open_spectrum(denoise_path)
