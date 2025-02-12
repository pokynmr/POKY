#
# This is an example script for concatenating two spectra.
# Spectra will be chosen from the running Poky session.
#
# This version uses the first spectrum as template, and add intensities 
#   from the second spectrum while concatenate_spectra_script.py requires
#   two spectrum parameters identifical. It means this is more flexible,
#   but slower. Still dimension order must be consistent. Be patient.
#
# One example use is to make a HNCACB + HNCA. 
# Because HNCACB is less sensitive, sometime a signal in HNCA doesn't exist
#   in HNCACB. If you don't want to have two spectra opened together,
#   you can use this script to generate a combined version.
#
# Also, we use noise level to match the scale difference automatically.
# If you wish to use different scale, you can modify "factor" down there.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jan. 12, 2025
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

specname = s.show_spectrumselectiondialog('Select a template spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

specname2 = s.show_spectrumselectiondialog('Select a spectrum to add', 0)
sp2 = name_to_spectrum(specname2, s)

if sp2 == None:
  raise SystemError

# check if dimension orders are the same between two spectra.
if sp.nuclei != sp2.nuclei:
  s.show_message('Error', 'Dimension orders are not the same.')
  raise SystemError

add_subtract = s.show_selectiondialog('Add or subtract', 
                                      'Do you want to add or subtract?',
                                      ('Add', 'Subtract', 'Cancel'))
if add_subtract in [-1, 2]:
  raise SystemError

# I think we don't need this for now.
#scale_offset = np.array(sp.scale_offset) - np.array(sp2.scale_offset)

import os.path
new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

# We will match scale factor by comparing noise level.
# You can manually adjust if you want something else.
factor = sp.noise / sp2.noise
if add_subtract == 1:
  factor = factor * -1.

proj = s.project

dic, data = ng.sparky.read(sp.data_path)
uc_list = []
for i in range(sp.dimension):
  uc = ng.sparky.make_uc(dic, data, i)
  ppm_scale = uc.ppm_scale()
  ppm_0, ppm_1 = uc.ppm_limits()
  uc_list.append( (uc, ppm_scale, ppm_0, ppm_1))

def get_ppm(uc_list, pts, dim):
  return uc_list[dim][0].ppm(pts)

data2 = np.zeros(data.shape)

# get relevant data heights from the second spectrum 
for i in range(data.shape[0]):
  for j in range(data.shape[1]):
    x_ppm = get_ppm(uc_list, i, 0)
    y_ppm = get_ppm(uc_list, j, 1)
    if sp.dimension == 2:
      # 2D
      hts = sp2.data_height( (x_ppm + sp.scale_offset[0], 
                    y_ppm + sp.scale_offset[1]) )
      data2[i, j] = hts
    else:
      for k in range(data.shape[2]):
        z_ppm = get_ppm(uc_list, k, 2)
        if sp.dimension == 3:
          # 3D
          hts = sp2.data_height( (x_ppm + sp.scale_offset[0], 
                    y_ppm + sp.scale_offset[1],
                    z_ppm + sp.scale_offset[2]) )
          data2[i, j, k] = hts
        else:
          # 4D
          for l in range(data.shape[3]):
            a_ppm = get_ppm(uc_list, l, 3)
            hts = sp2.data_height( (x_ppm + sp.scale_offset[0], 
                    y_ppm + sp.scale_offset[1],
                    z_ppm + sp.scale_offset[2],
                    a_ppm + sp.scale_offset[3]) )
            data2[i, j, k, l] = hts

# combine
data3 = np.array(data) + factor * data2
ng.sparky.write(new_path, dic, data3, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the concatenated spectrum?'):
  s.open_spectrum(new_path)
