#
# This is an example script to change sequence number in the PDB(x) file
# This script runs on POKY BUILD 02/06/21e or newer
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

input_file = s.open_filedialog('Select a PDB file', 
                    'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')
if input_file == '':
  raise SystemExit

output_file = s.save_filedialog('Save as', 
                    'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')

offset = s.show_inputdialog('Offset', 'What is the offset?', '0')
try:
  offset = int(offset)
except:
  raise SystemExit

from pymol import cmd
cmd.delete('for_renumber')
cmd.load(input_file, 'for_renumber')
cmd.alter('for_renumber', 'resi=str(int(resi)+%d)' % (offset))
cmd.save(output_file, 'for_renumber', state=0)
cmd.delete('for_renumber')

s.show_message('Finished', 'Finished.')