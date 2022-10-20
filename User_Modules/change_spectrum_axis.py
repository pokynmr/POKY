#
# This is an example script for changing dimension order of a spectrum.
# A spectrum will be chosen from the running Poky session.
# When processing, please do not work on spectral view. It may crash.
# Be patient on running this script.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Oct. 20, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import os
from sputil import name_to_spectrum

# POKY libraries
import __main__
s = __main__.main_session

specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

new_path = s.save_filedialog('New spectrum name', 'UCSF (*.ucsf);; Any (*)',
                                  os.path.dirname(sp.data_path))
if new_path == '':
  raise SystemError

# Ask dimensions to change.  
current_dim = ', '.join(sp.nuclei)


import ucsftool3
# ucsftool's swap is only between two dimensions.
# If more than two dimensions need to be swapped, we need to repeat.
ut = ucsftool3.ucsfTool()

# mapper. 2D is straight-forward. 4D needs to swap one by one.
# 2D
if sp.dimension == 2:
  ut.ucsf_open(sp.data_path)
  ut.write_swapped_axis(new_path, 1, 2, overwrite=1)
  ut.ucsf_close()
# 4D
elif sp.dimension == 4:
  ans = s.show_inputdialog('Swap Axis', 
    f'Current: {current_dim}\n' + \
    'Please specify two dimensions to swap (e.g. 1, 2)', '1, 2')
  dims = ans.split(',')
  try:
    dim1 = int(dims[0].strip())
    dim2 = int(dims[1].strip())
    if dim1 not in [1, 2, 3, 4] or dim2 not in [1, 2, 3, 4] or dim1 == dim2:
      print('Invalid selection of dimensions.')
      raise SystemError
    ut.ucsf_open(sp.data_path)
    ut.write_swapped_axis(new_path, dim1, dim2, overwrite=1)
    ut.ucsf_close()
  except:
    print('Invalid selection of dimensions.')
    raise SystemError
# 3D
elif sp.dimension == 3:
  mapper_3D = ('132', '213', '231', '312', '321', 'Cancel')
  mapping = s.show_selectionexdialog('Axis Order', 
    f'Current: {current_dim}\nSelect new axis order:', mapper_3D)
  if mapping == len(mapper_3D)-1:
    raise SystemError  
  new_order = mapper_3D[mapping]
  ut.ucsf_open(sp.data_path)
  # 132: 2 <-> 3
  if new_order == '132':
    ut.write_swapped_axis(new_path, 2, 3, overwrite=1)
  # 213: 1 <-> 2
  elif new_order == '213':
    ut.write_swapped_axis(new_path, 1, 2, overwrite=1)
  # 321: 1 <-> 3
  elif new_order == '321':
    ut.write_swapped_axis(new_path, 1, 3, overwrite=1)
  else:
    import tempfile
    import random
    tmp_file = os.path.join(tempfile.gettempdir(), 
          'tmp_' + str(random.randrange(10**6)) + '.ucsf')
    # 231: 1 <-> 3, 1 <-> 2
    if new_order == '231':
      ut.write_swapped_axis(tmp_file, 1, 3, overwrite=1)
      ut.ucsf_close()
      ut.ucsf_open(tmp_file)
      ut.write_swapped_axis(new_path, 1, 2, overwrite=1)
    # 312: 1 <-> 3, 2 <-> 3
    elif new_order == '312':
      ut.write_swapped_axis(tmp_file, 1, 3, overwrite=1)
      ut.ucsf_close()
      ut.ucsf_open(tmp_file)
      ut.write_swapped_axis(new_path, 2, 3, overwrite=1)
  ut.ucsf_close()

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(new_path)
