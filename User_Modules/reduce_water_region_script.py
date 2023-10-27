#
# This is an example script for reduce water region intensities of 
# selected spectrum in the project.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: October 27, 2023
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
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

dim_list = []
for i in range(len(sp.nuclei)):
  if sp.nuclei[i] == '1H':
    dim_list.append(f'w{i+1}')

if len(dim_list) == 0:
  raise SystemError

if len(dim_list) > 1:
  ans = s.show_selectionexdialog('Dimension', 'Dimension to apply: ',
                          dim_list)
  if ans == -1:
    raise SystemError
else:
  ans = 0

dimension = int(dim_list[ans][-1]) - 1

import os.path
new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

dic, data = ng.fileio.sparky.read(sp.data_path)

try:
  factor = float(s.show_inputdialog('Max Reduction',
              'Max Division: ', '10'))
except:
  s.show_message('Error', 'Value must be numeric.')
  raise SystemError

try:
  ans = s.show_inputdialog('Water signal region',
            'Water signal region to reduce', '3.7 - 5.7')
  a, b = ans.split('-')
  x1 = float(a.strip())
  x3 = float(b.strip())  
  x2 = (x1 + x3) / 2
  y1, y3, y2 = 1, 1, factor
except:
  s.show_message('Error', 'Values must be numeric.')
  raise SystemError

def quadratic(a, b, c, x):
  return a * x**2 + b * x + c

A = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
b = np.array([y1, y2, y3])
coef = np.linalg.solve(A, b)


uc = ng.sparky.make_uc(dic, data, dimension)
ppm_scale = uc.ppm_scale()
ppm_scale += sp.scale_offset[dimension]

if sp.dimension == 2:
  if dimension == 0:
    for i in range(data.shape[dimension]):
      data[i,:] /= max(1, quadratic(coef[0], coef[1], coef[2], 
                                    ppm_scale[i]))
  elif dimension == 1:
    for i in range(data.shape[dimension]):
      data[:,i] /= max(1, quadratic(coef[0], coef[1], coef[2], 
                                    ppm_scale[i]))
else:
  if dimension == 0:
    for i in range(data.shape[dimension]):
      data[i,:,:] /= max(1, quadratic(coef[0], coef[1], coef[2], 
                                    ppm_scale[i]))
  elif dimension == 1:
    for i in range(data.shape[dimension]):
      data[:,i,:] /= max(1, quadratic(coef[0], coef[1], coef[2], 
                                    ppm_scale[i]))
  elif dimension == 2:
    for i in range(data.shape[dimension]):
      data[:,:,:i] /= max(1, quadratic(coef[0], coef[1], coef[2], 
                                    ppm_scale[i]))
  
ng.sparky.write(new_path, dic, data, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(new_path)

