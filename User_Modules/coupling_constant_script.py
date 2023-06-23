#
# This is an example script to calculate coupling constants in a spectrum.
# You must have a pair of same assignment in the spectrum.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#

import os
import __main__
s = __main__.main_session
spec_name = s.show_spectrumselectiondialog('Select a spectrum', 0)
if spec_name == '':
  raise SystemExit    

# spectrum to be simulated
from sputil import name_to_spectrum
spec = name_to_spectrum(spec_name, s)
if spec == None:
  s.show_message('Error', 'Spectrum is not set properly.')
  raise SystemExit

peak_dict = {}
for p in spec.peak_list():
  if p.is_assigned != 1:
    continue
  try:
    peak_dict[p.assignment].append(p)
  except:
    peak_dict[p.assignment] = [p,]

if spec.dimension == 2:
  content = '      Assignment        dw1(Hz)    dw2(Hz)\n'
elif spec.dimension == 3:
  content = '      Assignment        dw1(Hz)    dw2(Hz)       dw3(Hz)\n'

print(content)
for key in peak_dict.keys():
  if len(peak_dict[key]) != 2:
    continue
  p, p2 = peak_dict[key]
  dw_list = list(map(lambda dim: abs(p.frequency[dim] - p2.frequency[dim]), 
                      range(spec.dimension)))
  line = '%18s' % (p.assignment)
  for dim in range(spec.dimension):
    line += ' %10.3f' % (dw_list[dim] * spec.hz_per_ppm[dim])
  print(line)
  content += line
  content += '\n'

s.set_clipboard(content)
if s.show_message_yes_no('Finished',
      'Coupling constant data is in your clipboard. ' + \
      'Do you want to open in your spreadsheet program?'):
  lines = content.split('\n')
  csv_content = ''
  for line in lines:
    csv_content += ','.join(line.split())
  import tempfile
  import random
  tmp_file = os.path.join(tempfile.gettempdir(), 
              'tmp_' + str(random.randrange(10**6)) + '.csv')
  f = open(tmp_file, 'w')
  f.write(content)
  f.close()

  s.open_external(tmp_file)