#
# This is an example script to renumber the following files. 
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Supported files: pdb, upl, aco, rdc, prot, str, seq
#
# Last Update: June 14, 2023
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY renumbering tool')
print('created by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

#################
import os, sys

ext_list = ('pdb', 'upl', 'tbl', 'aco', 'rdc', 'prot', 'str', 'seq')
desc_list = ('PDB coordinate file', 'DIANA upper limit file', 'XPLOR TBL file',
            'DIANA angle constraint file', 'DIANA residual dipolar coupling file',
            'XEASY PROT chemical shift file', 'BMRB NMR-STAR file', 'SEQ file')
extfilter = 'Any (*)'

for i in range(len(ext_list)):
  extfilter = extfilter + f';; {desc_list[i]} (*.{ext_list[i]})'

input_file = s.open_filedialog('Select a file to read', 
                            extfilter, '')
if input_file == '':
  raise SystemExit
file_ext = input_file.split('.')[-1]

try:
  ext_idx = ext_list.index(file_ext)
except:
  s.show_message('Error', 'Not supported file type.')
  raise SystemExit

output_file = s.save_filedialog('Select a file to write', 
  f'{desc_list[ext_idx]} (*.{ext_list[ext_idx]})', os.path.dirname(input_file))
if output_file == '':
  raise SystemExit

d = s.show_inputdialog('The Shift', 'Type the differential number', '1')

try:
  diff = int(d)
except:
  raise SystemExit

if diff == 0 or abs(diff) > 999:
  raise SystemExit

f = open(input_file, 'r')
lines = f.readlines()
f.close()

content = ''

# PDB file
if file_ext == 'pdb':
  for line in lines:
    try:
      seqnum = line[22:26]
      nseqnum = int(seqnum.strip())
      nseqnum2 = nseqnum+diff
      seqnum2 = '%4d' % (nseqnum2)
      line2 = line[0:22] + seqnum2 + line[26:]
      content = content + line2
    except:
      content = content + line
elif file_ext == 'aco':
  for line in lines:
    try:
      seqnum = line[0:4]
      nseqnum = int(seqnum.strip())
      nseqnum2 = nseqnum+diff
      seqnum2 = '%4d' % (nseqnum2)
      line2 = seqnum2 + line[4:]
      content = content + line2
    except:
      content = content + line
elif file_ext == 'upl':
  for line in lines:
    try:
      slist = line.strip().split()
      seqnum = slist[0]
      seqnum2 = slist[3]
      nseqnum = int(seqnum)
      nseqnum2 = int(seqnum2)
      nseqnum3 = nseqnum+diff
      nseqnum4 = nseqnum2+diff      
      seqnum3 = '%4d' % (nseqnum3)
      seqnum4 = '%4d' % (nseqnum4)      
      line2 = '%-5s %-5s %-5s %-5s %-5s %-5s %-s\n' % \
          (seqnum3, slist[1], slist[2], seqnum4, slist[4], slist[5], slist[6])
      content = content + line2
    except:
      content = content + line
elif file_ext == 'prot':
  for line in lines:
    try:
      slist = line.strip().split()
      seqnum = slist[-1]
      nseqnum = int(seqnum)
      nseqnum2 = nseqnum + diff      
      seqnum2 = '%4d' % (nseqnum2)
      line2 = '%-7s %-7s %-7s %-7s %-s\n' % \
          (slist[0], slist[1], slist[2], slist[3], seqnum2)
      content = content + line2
    except:
      content = content + line
elif file_ext == 'rdc':
  bFirst = False
  for line in lines:
    sp = line.split()
    if not bFirst: 
      content = content + line
      if line.find('First') != -1: 
        bFirst = True
      continue
    else:  
      try:
        nseqnum = int(sp[0])
        nseqnum2 = int(sp[3])
        nseqnum3 = nseqnum+diff
        nseqnum4 = nseqnum2+diff      
        #seqnum3 = '%6d' % (nseqnum3)
        #seqnum4 = '%6d' % (nseqnum4)      
        #line2 = seqnum3 + line[6:17] + seqnum4 + line[23:]
        line2 = ''
        for i in range(len(sp)):
          if i == 0:
            line2 = line2 + '%6d' % (nseqnum3)
          elif i == 3: 
            line2 = line2 + '%6d' % (nseqnum4)
          else:
            line2 = line2 + '%8s' % (sp[i])
        content = content + line2 + '\n'
      except:
        content = content + line
elif file_ext == 'tbl':
  for line in lines:
    new_line = ''
    skip = -1
    for i in range(len(line)):
      if skip != -1 and skip > i:
        continue
      if line[i:i+len('resid')] == 'resid':
        tmp = ''
        for j in range(i+6, i+16):
          if line[j] == 'a' or line[j] == 0:
            break
          tmp += line[j]
        new_line += 'resid'
        new_line += ' %4d' % (int(tmp.strip()) + diff)
        skip = j-1
      else:
        new_line += line[i]
    content = content + new_line + '\n'
elif file_ext == 'str':
  for line in lines:
    try:
      seqnum = line[5:9]
      nseqnum = int(seqnum.strip())
      nseqnum2 = nseqnum+diff
      seqnum2 = '%4d' % (nseqnum2)
      line2 = line[0:5] + seqnum2 + line[9:]
      if len(line.strip().split()) < 3:
        content = content + line2
        continue
    except:
      pass
    try:
      seqnum = line[8:11]
      seqnum2 = line[13:16]      
      nseqnum = int(seqnum.strip())
      nseqnum2 = int(seqnum2.strip())      
      nseqnum3 = nseqnum+diff
      nseqnum4 = nseqnum2+diff      
      seqnum3 = '%3d' % (nseqnum3)
      seqnum4 = '%3d' % (nseqnum4)      
      line2 = line[0:8] + seqnum3 + line[11:13] + seqnum4 + line[16:]
      content = content + line2
    except:
      content = content + line
elif file_ext == 'seq':
  for line in lines:
    try:
      sp_list = line.strip().split()
      nseq = int(sp_list[1])+diff
      line2 = '%6s %6d\n' % (sp_list[0], nseq)
      content = content + line2
    except:
      content = content + line

f = open(output_file, 'w')
f.write(content)
f.close()

if s.show_message_yes_no('Finished', 
      'Done. Do you want to open in Poky Notepad?'):
  s.open_notepad(output_file)
