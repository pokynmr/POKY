#
# This is an example script to convert pzst to ucsf file.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

print('\n\n\n------------------------------------------------------')
print('POKY PZST2UCSF conversion')
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

pzstpath = s.open_filedialog('Select a pzst file', 'PZST file (*.pzst)', '')
if pzstpath == '':
  raise SystemError

ucsfpath = pzstpath[:-4] + 'ucsf'
ucsfpath = s.save_filedialog('Save as', 'UCSF file (*.ucsf)', ucsfpath)
if ucsfpath == '':
  raise SystemError

ans = s.show_inputdialog('Smooth sigma', 
      'Type smooth sigma between 0 and 1. (Default: 0.7)',
      '0.7')
smooth_sigma = float(ans)
if not 0.0 <= smooth_sigma <= 1:
  raise SystemError

poky_pzst.pzst2ucsf(pzstpath, ucsfpath, smooth_data=True, 
                    smooth_sigma=smooth_sigma)
if s.show_message_yes_no('Load data', 
                         'Do you want to load the converted spectrum?'):
  s.open_spectrum(ucsfpath)