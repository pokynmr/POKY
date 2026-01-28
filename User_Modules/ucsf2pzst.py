#
# This is an example script to convert ucsf to compressed pzst file.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

print('\n\n\n------------------------------------------------------')
print('POKY UCSF2PZST conversion')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import os
import __main__
s = __main__.main_session

try:
  import poky_pzst
except:
  s.show_message('Error', 'Please upgrade your POKY to use this feature.')
  raise SystemError

import nmrglue as ng
import numpy as np

ucsfpath = s.open_filedialog('Select a UCSF file', 'UCSF file (*.ucsf)', '')
if ucsfpath == '':
  raise SystemError

pzstpath = ucsfpath[:-4] + 'pzst'
pzstpath = s.save_filedialog('Save as', 'PZST file (*.pzst)', pzstpath)
if pzstpath == '':
  raise SystemError

# ask bit depth
ans = s.show_selectionexdialog('Choose the bit depth for the export. ',
        'A larger bit depth results in more detail and a larger file size.',
        ('8', '16', '32', 'Cancel'))
if ans in [-1, 3]:
  raise SystemError
bit_depth = 2**(3+ans)

# ask noise scale
ans = s.show_selectionexdialog('Choose the noise scale for the export. ',
    'A larger scale results in lossier and a smaller file size. (Default: 1)',
    ('0', '0.5', '1', '2', '3', 'Cancel'))
if ans in [-1, 5]:
  raise SystemError
noise_scale = [0, 0.5, 1, 2, 3][ans]

dic, data = ng.sparky.read(ucsfpath)
flat_data = data.flatten()
random_indices = np.random.choice(flat_data.size, size=30, replace=False)
noise_estimate = np.median(np.abs(flat_data[random_indices]))
print(f"Estimated Noise: {noise_estimate}")

print(f'--- Compressing {ucsfpath} ---')
poky_pzst.ucsf2pzst(ucsfpath, pzstpath, bit_depth=bit_depth,
                    noise_level=noise_estimate, 
                    noise_threshold=noise_estimate * noise_scale)

s.show_message('Done.', f'{pzstpath} created.')