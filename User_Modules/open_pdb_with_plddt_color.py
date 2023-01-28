#
# This is an example script for opening a PDB with pLDDT coloring.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jan. 28, 2023
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import re
import os

# POKY libraries
import __main__
s = __main__.main_session
      
in_name = s.open_filedialog('Open PDB', 'PDB (*.pdb);; mmCIF (*.mmcif);; Any (*)', '')
if in_name == '':
  raise SystemError

f = open(in_name, 'r')
lines = f.readlines()
f.close()
lines = list(filter(lambda x: x.find('ATOM') == 0 and len(x) > 65, lines))
b_list = list(map(lambda x: x[60:66].strip(), lines))
b_list = list(map(lambda x: float(x), b_list))
if max(b_list) < 1:
  cutoffs = ['0.90', '0.70', '0.50']
else:
  cutoffs = ['90', '70', '50']

objname = os.path.splitext(os.path.basename(in_name))[0]  
c_list = ['0x0053d7', '0x57caf9', '0xffdb12', '0xff7e45']
cmd = f'load {in_name}; '
cmd += f'color {c_list[0]}, ({objname}) and (b >{cutoffs[0]} or b ={cutoffs[0]}); '
cmd += f'color {c_list[1]}, ({objname}) and ((b <{cutoffs[0]} and b >{cutoffs[1]}) or (b ={cutoffs[1]})); '
cmd += f'color {c_list[2]}, ({objname}) and ((b <{cutoffs[1]} and b >{cutoffs[2]}) or (b ={cutoffs[2]})); '
cmd += f'color {c_list[2]}, ({objname}) and ((b <{cutoffs[2]} and b >0.0 ) or (b =0.0))'
s.set_clipboard(f'{cmd}')
s.show_message('PyMOL', 'Ctrl+V into PyMOL commandline.')
s.open_pymol('')