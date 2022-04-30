#
# This is an example script for making a POKY resonance list
# from a tabular file.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: April 30, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

##############################
# USER PARAMETERS TO SET
# Columns in the tabular file
# Column starts from zero (0)
col_nseq = 1    # RESIDUE NUMBER
col_aaa = 2     # THREE LETTER A.A.
col_a = -1      # ONE LETTER A.A.
col_atm = 3     # ATOM
col_cs = 6      # CHEMICAL SHIFT
#
##############################

# POKY libraries
import __main__
s = __main__.main_session
from wlutil import AAA_dict

nuc_dict = {'H': '1H', 'C': '13C', 'N': '15N', 'F': '19F', 'P': '31P'}

tab_file = s.open_filedialog('Select a tabular file', 
                'Any (*)', '')

if tab_file == '':
  raise SystemError

import os
from sys import platform
tab_path = os.path.dirname(tab_file)

out_file = s.save_filedialog('Save as', 
                'shifts.list (*.list);; Any (*)', '')


if out_file == '':
  raise SystemError

f = open(tab_file, 'r')
lines = f.readlines()
f.close()

text = ''
for line in lines:
  try:
    seg_list = line.strip().split()
    if col_a == -1:
        a = AAA_dict[seg_list[col_aaa]]
    else:
        a = seg_list[col_a]
    nseq = seg_list[col_nseq]
    atm = seg_list[col_atm]
    nuc = nuc_dict[atm[0]]
    cs = seg_list[col_cs]
    data = '%s%s %5s %4s %9s   0.0   0\n' % (a, nseq, atm, nuc, cs)
    text += data
  except:
    print('WARNING: ' + line)
    continue

f = open(out_file, 'w')
f.write(text)
f.close()

s.show_message('Finished.', 'Load the output file in the resonance list window.')

out_path = os.path.dirname(out_file)
if 'linux' in platform:
  os.system(f'xdg-open {out_path}')
elif platform == 'darwin':
  os.system(f'open {out_path}') # OS X
elif platform == 'win32' or platform == 'cygwin':
  os.startfile(out_path) # Windows
if os.path.getsize(out_file) > 0:
  s.open_notepad(out_file)
  s.command_characters('rl')
