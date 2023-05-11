#
# This is an example script to take 2D slices from a 3D spectrum.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#
# Runs on BUILD 02/13/2023g or newer
#

import __main__
s = __main__.main_session
spec_name = s.show_spectrumselectiondialog('Select a spectrum', 0)
if spec_name == '':
  raise SystemExit      

from sputil import name_to_spectrum
spec = name_to_spectrum(spec_name, s)
if spec == None:
  s.show_message('Error', 'Spectrum is not set properly.')
  raise SystemExit

if spec.dimension != 3:
  s.show_message('Error', 'Please select 3D spectrum.')
  raise SystemExit

idx = s.show_selectionexdialog('Dimension', 
                    'Select the dimension to slice.', 
                    ('w1', 'w2', 'w3', 'Cancel'))
if idx in [3, -1]:
  raise SystemExit

import os
spec_dir = os.path.dirname(os.path.abspath(spec.data_path))
out_dir = s.open_directorydialog('Select directory to save', spec_dir)
if out_dir == '':
  raise SystemExit

import ucsftool3
ut = ucsftool3.ucsfTool()
ut.ucsf_open(spec.data_path)
try:
  ut.write_planes('slice', out_folder=out_dir, dim=idx+1, overwrite=1)
except:
  os.chdir(out_dir)
  ut.write_planes('slice', dim=idx+1, overwrite=1)
ut.ucsf_close()

s.show_message('Finished.', 'Finished.')