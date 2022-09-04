#
# This is an example script to add hydrogen atoms to PDB.
# This will remove heavy atoms if not N, C, O, S
# This script runs on POKY BUILD 02/06/21e or newer
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

in_file = s.open_filedialog('Select the PDB ensemble file', 
                            'Any (*);; mmCIF (*.cif);; PDB (*.pdb)', '')
if in_file == '':
  raise SystemExit

out_file = s.save_filedialog('Save as', 'PDB (*.pdb);; Any (*)', '')
if out_file == '':
  raise SystemExit

from pymol import cmd
cmd.load(in_file, 'for_hydrogen')
cmd.do('remove not name N*+C*+O*+S*')
cmd.do('h_add for_hydrogen')
cmd.do('save ' + out_file)

# rename protons
from myseq import AAA_dict, A_dict, aaa2a

f = open(out_file, 'r')
lines = f.readlines()
f.close()

H12_H13 = 2
wlines, wlines2 = [], []
iFirstSeq, iCurSeq = -9999, -9999
tenSeq = '' 
for i in range(len(lines)):
  line = lines[i]
  if len(line) < 66: continue
  if line[0:4] != 'ATOM': continue
  heavyatom = line[12:16].strip()
  if iFirstSeq == -9999:
    iFirstSeq = int(line[22:26].strip())

  iSeq = int(line[22:26].strip())
  if (iCurSeq != iSeq) and len(tenSeq) < 10:
    tenSeq = tenSeq + aaa2a(line[17:20])
    iCurSeq = iSeq

  if heavyatom[0] == 'H': 
    continue

  line = line[0:12] + ('%4s' % (heavyatom.strip())) + line[16:]
  line = line[0:6] + ('%5d' % (len(wlines) + 1)) + line[11:]
  wlines.append(line)
      
  # find matching protons
  x = float(line[30:38].strip())
  y = float(line[38:46].strip())
  z = float(line[46:54].strip())
    
  iMatch = 0
  match_list = []
    
  for j in range(i+1,len(lines)):
    line2 = lines[j]
    if len(line2) < 66: break
    if line2[0:4] != 'ATOM': break
    if line[22:26].strip() != line2[22:26].strip():
      break
    if line2[12:16].strip()[0] != 'H':
      continue
      
    # distance under 1.25
    x2 = float(line2[30:38].strip())
    y2 = float(line2[38:46].strip())
    z2 = float(line2[46:54].strip())
    if (x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2) < 1.25*1.25:
      match_list.append(line2)

  for j in range(len(match_list)):
    line2 = match_list[j]
    protonatom = 'H' + heavyatom[1:]
    if len(match_list) == 2:
      protonatom = '%s%d' % (protonatom, j+H12_H13)
    elif len(match_list) != 1:
      protonatom = '%s%d' % (protonatom, j+1)      
    line2 = line2[0:12] + ('%4s' % (protonatom.strip())) + line2[16:]
    line2 = line2.replace('CYSS', 'CYS ')
    line2 = line2[0:6] + ('%5d' % (len(wlines) + 1)) + line2[11:]
    wlines.append(line2)
    
f = open(out_file,'w')
for line in wlines:
  f.write(line)    
f.close()
print('DONE')
