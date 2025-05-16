#
# This is an example script for normalizing the spectrum against the reference.
# Spectra will be chosen from the running Poky session.
#
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: May 15, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import os.path
import numpy as np
import nmrglue as ng
from sputil import name_to_spectrum

# POKY libraries
import __main__
s = __main__.main_session

specname = s.show_spectrumselectiondialog('Select a reference spectrum', 0)
ref_sp = name_to_spectrum(specname, s)

if ref_sp == None:
  raise SystemError
ref_nuc = ''.join(ref_sp.nuclei)

specname = s.show_spectrumselectiondialog('Select a target spectrum to normalize', 0)
tar_sp = name_to_spectrum(specname, s)

if tar_sp == None:
  raise SystemError
tar_nuc = ''.join(tar_sp.nuclei)

if ref_nuc != tar_nuc:
  s.show_message('Error', f'Nuclei mismatch. {ref_nuc} != {tar_nuc}')
  raise SystemError

new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(tar_sp.data_path))
if new_path == '':
  raise SystemError

mode = s.show_selectionexdialog('Normalization', 'Select a mode', 
            ('Scale range', 'Standardize using mean and deviation', 'Cancel'))
if mode in [-1, 2]:
  raise SystemError


def scale_to_range(array_to_normalize, reference_array):
  min_ref = np.min(reference_array)
  max_ref = np.max(reference_array)
  min_norm = np.min(array_to_normalize)
  max_norm = np.max(array_to_normalize)

  if max_norm == min_norm:
    return np.full_like(array_to_normalize, (min_ref + max_ref) / 2)

  scaled_array = min_ref + (array_to_normalize - min_norm) * \
                (max_ref - min_ref) / (max_norm - min_norm)
  return scaled_array

def standardize_against(array_to_normalize, reference_array):
  mean_ref = np.mean(reference_array)
  std_ref = np.std(reference_array)

  if std_ref == 0:
    return np.full_like(array_to_normalize, mean_ref)

  standardized_array = (array_to_normalize - mean_ref) / std_ref
  return standardized_array

ref_dic, ref_data = ng.sparky.read(ref_sp.data_path)
tar_dic, tar_data = ng.sparky.read(tar_sp.data_path)

if mode == 0:
  data = scale_to_range(tar_data, ref_data)
elif mode == 1:
  data = standardize_against(tar_data, ref_data)

ng.sparky.write(new_path, tar_dic, data, overwrite=True)

if s.show_message_yes_no('Load data', 
              'Do you want to load the scaled spectrum?'):
  s.open_spectrum(new_path)