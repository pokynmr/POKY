#
# This is an example script for assigning a chain ID to PDB.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jan. 30, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

chainID = 'A'

import __main__
s = __main__.main_session
in_file = s.open_filedialog('message', 'Any (*);; PDB file (*.pdb)', '')
if in_file == '':
  raise SystemError

out_file = s.save_filedialog('message', 'Any (*);; PDB file (*.pdb)', '')
if out_file == '':
  raise SystemError

f = open(in_file, 'r')
lines = f.readlines()
f.close()

content = ''

for line in lines:
  if not line.startswith('ATOM'):
    content += line
    continue
  line = line[:21] + chainID + line[22:]
  content += line
f = open(out_file, 'w')
f.write(content)
f.close()

s.show_message('Saved.', 'Assignment completed in ' + out_file)