#
# This is an example script to change sequence numbers of NDPPLOT INI file.
# This script only replaces X numbers not the labels!
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session
input_ini = s.open_filedialog('Select a NDP-PLOT INI file',
                  'NDPPLOT INI (*.ini);; Any (*)', '')
if input_ini == '':
  raise SystemExit

output_ini = s.save_filedialog('Save as',
                  'NDPPLOT INI (*.ini);; Any (*)', '')
if output_ini == '':
  raise SystemExit

n = s.show_inputdialog('Renumber', 'What will be the first residue number?', '1')
user_first = int(n)

f = open(input_ini,'r')
lines = f.readlines()
f.close()

input_first = 1
for line in lines:
  if line.find('XMIN=') != -1:
    input_first = int(line.strip().split('=')[1])
    break
d = user_first-input_first

text = ''
for line in lines:
  if line.find('XMIN=') != -1 or line.find('XMAX=') != -1:
    s1 = line.strip().split('=')
    tmp = '%s=%d' % (s1[0], int(s1[1])+d)
    text += tmp
    text += '\n'
    continue
  if line.find('COLOR_') != -1:
    s1 = line.strip().split('=')
    s2 = s1[1].split(',')
    s3 = '%s=%d,%d' % (s1[0], int(s2[0])+d, int(s2[1])+d)
    for seg in s2[2:]:
      s3 = s3 + ',' + seg
    text += s3
    text += '\n'
    continue
  if line.find('DATA_') != -1:
    s1 = line.strip().split('=')
    s2 = s1[1].split(',')
    s3 = '%s=%d' % (s1[0], int(s2[0])+d)
    for seg in s2[1:]:
      s3 = s3 + ',' + seg
    text += s3
    text += '\n'
    continue
  text += line

f = open(output_ini, 'w')
f.write(text)
f.close()

s.show_message('Finished', 'Finished.')