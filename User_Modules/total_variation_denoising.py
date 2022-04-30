#
# This is an example script for applying total-variation denoising using 
# split-Bregman optimization on the selected spectrum in the project.
# https://en.wikipedia.org/wiki/Total_variation_denoising
#
# It can take pretty long depending on the data size.
# Make sure you have saved all the data before using this.
#
# Modified from the "wavelet denoising" script by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Developed by Mehdi Rahimi, Ph.D. (mehdi.rahimi@ucdenver.edu)
#
# Last update: February 19, 2021 
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
from skimage.restoration import denoise_tv_bregman
from sputil import name_to_spectrum
# POKY libraries
import __main__
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

import os.path
noisy_path = s.save_filedialog('Provide a name for saving the noisy spectrum', 'UCSF (*.ucsf);; Any (*)',
                                os.path.dirname(sp.data_path))
if noisy_path == '':
  raise SystemError


denoise_path = s.save_filedialog('Provide a name for saving the de-noised spectrum', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if denoise_path == '':
  raise SystemError
  
  
dic, data = ng.sparky.read_lowmem(sp.data_path)

# Creating a noisy data by adding a Gaussian noise to a selected spectrum
noisy_data = data + np.random.normal(scale=sp.noise*80, size=data.shape)

ng.sparky.write_lowmem(noisy_path, dic, noisy_data, overwrite=True)


# De-noising the noisy data using the tv bregman method
data_denoise = denoise_tv_bregman(noisy_data, 1e-5, max_iter=100, eps=0.001, isotropic=True)

ng.sparky.write_lowmem(denoise_path, dic, data_denoise, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectra?'):
  s.open_spectrum(noisy_path)
  s.open_spectrum(denoise_path)


