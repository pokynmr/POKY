#
# This is an example script to export peak intensities in csv.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY Export Peak Intensities')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

spec_names = s.show_spectrumselectiondialog('Select spectra', 1)
spec_name_list = spec_names.split('\t')
if len(spec_name_list) == 0:
  raise SystemError

import os
from sputil import name_to_spectrum

spec_list = list(map(lambda x: name_to_spectrum(x, s), spec_name_list))
if None in spec_list:
  raise SystemError

csv_dir = s.open_directorydialog('Export into', '')

for spec in spec_list:
  print(f'Process {spec.name}')
  peaks = spec.peak_list()
  content = ''
  for peak in peaks:
    strfreqs = ','.join([str(x) for x in peak.frequency])
    line = f'{peak.assignment},{strfreqs},{peak.data_height}\n'
    content += line
  csv_path = os.path.join(csv_dir, f'{spec.name}.csv')
  print(f'Writing {csv_path}')
  f = open(csv_path, 'w')
  f.write(content)
  f.close()
print('Done.')