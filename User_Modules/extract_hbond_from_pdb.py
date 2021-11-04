#
# Extract hydrogen bond constraints using a PDB file. 
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Last Update: Nov 4, 2021
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('Extract hydrogen bond constraints using a PDB file.')
print('created by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

#################
import os, sys

pdb_file = s.open_filedialog('Select a file to read', 
                            'PDB file (*.pdb)', '')
if pdb_file == '':
  raise SystemExit

upl_file = s.save_filedialog('Select a file to write', 
      'DIANA upper limit file (*.upl)', os.path.dirname(pdb_file))
if upl_file == '':
  raise SystemExit
lol_file = upl_file[:-3] + 'lol'

if not os.path.exists(pdb_file):
  s.show_message('Error', '%s not found.' % (pdb_file))
  sys.exit()

f = open(pdb_file, 'r')
pdb_lines = f.readlines()
f.close()

"""
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
"""
atom_list = []
seq_list = []
for i in range(1000): # N H O
  atom_list.append([[9999,9999,9999], [9999,9999,9999], [9999,9999,9999]])
  seq_list.append('')

atm_dic = {'N':0, 'H':1, 'HN':1, 'O':2}

for line in pdb_lines:
  if line.find('ENDMDL') == 0:
    break # only model 1 is supported.
  if line.find('ATOM') != 0 or len(line) < 55: 
    continue
  atm = line[12:16].strip()
  if not atm in ['N', 'H', 'HN', 'O']: 
    continue
  seq, x, y, z = int(line[22:26].strip()), float(line[30:38].strip()), \
    float(line[38:46].strip()), float(line[46:54].strip())
  atom_list[seq][atm_dic[atm]][0], atom_list[seq][atm_dic[atm]][1], \
    atom_list[seq][atm_dic[atm]][2], seq_list[seq] = x, y, z, \
                                          line[17:20].strip()
    
content = ''    
for i in range(len(atom_list)):
  n,h,o = atom_list[i]
  for j in range(len(atom_list)):
    if abs(i-j) < 4: continue
    n2,h2,o2 = atom_list[j]
    if [9999,9999,9999] in [n, h, o, n2, h2, o2]:
      continue
    no_dist = (n[0]-o2[0])**2+(n[1]-o2[1])**2+(n[2]-o2[2])**2
    ho_dist = (h[0]-o2[0])**2+(h[1]-o2[1])**2+(h[2]-o2[2])**2
    if no_dist > 3.4**2 or ho_dist > 2.3**2: continue
    if no_dist < 2.6**2 or ho_dist < 1.7**2: continue
    #hLine = 'assign ( resid  %3d and name HN ) ( resid %3d and name O ) 2.0 0.3 0.3' % (i, j)
    #nLine = 'assign ( resid  %3d and name  N ) ( resid %3d and name O ) 3.0 0.3 0.3' % (i, j)
    hLine = '%3d %3s H %3d %3s O 2.2' % (i, seq_list[i], j, seq_list[j])
    nLine = '%3d %3s N %3d %3s O 3.3' % (i, seq_list[i], j, seq_list[j])
    content = content + hLine + '\n'
    content = content + nLine + '\n'
f = open(upl_file, 'w')
f.write(content)
f.close()

if s.show_message_yes_no('Finished', 
      'Done. Do you want to open in Poky Notepad?'):
  s.open_notepad(upl_file)