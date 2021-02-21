
#
# This is an example script to convert between DIANA UPL and XPLOR TBL
# TBL file must be single line single restraint.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

from atomnames import s12s32
from myseq import A_dict, AAA_dict

import re

def one_or_three(seq_file):
  f = open(seq_file, 'r')
  tmp = f.read()
  f.close()
  tmp2 = tmp.upper()
  regex = re.compile('[^A-Z]')
  seq = regex.sub('', tmp2)
  if len(seq) % 3 != 0:    return 3, seq
  for i in range(int(len(seq) / 3)):
    subseq = seq[i*3:i*3+3]
    try:
      AAA_dict[subseq]
    except:
      return 1, seq
  return 3, seq

def extract_seq_file(seq_file):
  seq_mode, seq = one_or_three(seq_file)
  a_list = list(map(lambda x: '', range(1000)))
  if seq_mode == 1:
    for i in range(len(seq)): a_list[i] = seq[i]
    return a_list
  regex = re.compile('[^A-Z+0-9+ ]')
  f = open(seq_file, 'r')
  lines = f.readlines()
  f.close()
  iCur = 1
  for line in lines:
    if line.strip() == '': continue
    if line.strip()[0] == '#': continue
    tmp = regex.sub('', line.upper()).split()
    if len(tmp) == 2:
      try:
        iCur = int(tmp[0])
        a_list[iCur] = AAA_dict[tmp[1]]
        iCur = iCur + 1
      except:
        iCur = int(tmp[1])
        a_list[iCur] = AAA_dict[tmp[0]]
        iCur = iCur + 1
        continue
    if len(tmp) == 1:
      a_list[iCur] = AAA_dict[tmp[0]]
      iCur = iCur + 1
      continue

  return a_list
  
def translate12to23(aaa, atm):
  if atm == 'HN': return 'H'
  if atm == 'CO': return 'C'
  if len(aaa) == 3:     a = AAA_dict[aaa]
  else:                 a = aaa
  for grp in s12s32:
    if (grp[0] == a) and (grp[1] == atm):  return grp[2]
  return atm

def translate23to12(aaa, atm):
  if atm == 'H': return 'HN'
  if atm == 'C': return 'CO'
  if len(aaa) == 3:     a = AAA_dict[aaa]
  else:                 a = aaa
  for grp in s12s32:
    if (grp[0] == a) and (grp[2] == atm):  return grp[1]
  return atm

def xplortbl2dyanaupl(tbl_line, a_list):
  tmp3 = tbl_line.strip().replace('assign', '').replace('resid', '')
  tmp2 = tmp3.replace('name','').replace('and','').replace('(','').replace(')','')
  if tmp2 == '': return ''
  if tmp2[0] == '!' or tmp2[0] == '#' or tmp2[0] == '/': return ''
  tmp = tmp2.replace('#','').replace('*','')

  num_list = []
  atm_list = []
  sp = tmp.split()
  if len(sp) != 7: return ''
  for s in sp:
    try:
      num_list.append(float(s))
    except:
      atm_list.append(s)
      
  if len(num_list) != 5 or len(atm_list) != 2: return ''
  n1 = int(num_list[0])
  n2 = int(num_list[1])
  a1 = translate12to23(a_list[n1], atm_list[0])
  a2 = translate12to23(a_list[n2], atm_list[1])
  s1 = A_dict[a_list[n1]]
  s2 = A_dict[a_list[n2]]
  dist = num_list[2] + num_list[4]
  upl_line = '%4d %4s %4s %4d %4s %4s %8.2f' % (n1, s1, a1, n2, s2, a2, dist)
  return upl_line

def dyanaupl2xplortbl(upl_line):
  tmp = upl_line.strip().replace('CQ','C').replace('QQ','H').replace('Q','H')
  if tmp == '': return ''
  if tmp[0] == '!' or tmp[0] == '#' or tmp[0] == '/': return ''
  sp = tmp.split()
  if len(sp) < 7: return ''
  try:
    n1 = int(sp[0])
    n2 = int(sp[3])
    s1 = sp[1]
    s2 = sp[4]
    a1 = translate23to12(s1, sp[2])
    a2 = translate23to12(s2, sp[5])
    if len(a1) > 1 and a1 != 'HN': a1 = a1 + '*'
    if len(a2) > 1 and a2 != 'HN': a2 = a2 + '*'
    dist = float(sp[6])
    if a1[0] == 'H' or a2[0] == 'H':    lbn = dist-1.8
    elif a1[0] == 'C' and a2[0] == 'C': lbn = dist-4.0
    else: lbn = dist-2.7

    tbl_line = \
        'assign ( resid %3d and name %4s ) ( resid %3d and name %4s ) %6.2f %6.2f   0.00' % \
        (n1, a1, n2, a2, dist, lbn)
    return tbl_line
  except:
    return ''

def tbl2upl(tbl_file, seq_file):
  f = open(tbl_file, 'r')
  lines = f.readlines()
  f.close()
  text = ''
  a_list = extract_seq_file(seq_file)
  for line in lines:
    upl = xplortbl2dyanaupl(line, a_list)
    if upl != '': 
      text += upl
      text += '\n'
  return text

def upl2tbl(upl_file):
  f = open(upl_file, 'r')
  lines = f.readlines()
  f.close()
  text = ''
  for line in lines:
    tbl = dyanaupl2xplortbl(line)
    if tbl != '': 
      text += tbl
      text += '\n'
  return text

import os.path
import __main__
s = __main__.session
dist_file = s.open_filedialog('Select a distance restraint file.', 
    'DIANA UPL (*.upl);; XPLOR TBL (*.tbl)', '')
if dist_file == '':
  raise SystemExit

p = os.path.dirname(dist_file)
if dist_file.endswith('.upl'):
  target_file = s.save_filedialog('Save as...', 
    'XPLOR TBL (*.tbl);; Any (*)', p)
  if target_file == '':
    raise SystemExit
  text = upl2tbl(dist_file)
elif dist_file.endswith('.tbl'):
  seq_file = s.open_filedialog('Select a sequence file.', 
    'Sequence file (*.seq *.fasta *.txt);; Any (*)', p)
  target_file = s.save_filedialog('Save as...', 
    'DIANA UPL (*.upl);; Any (*)', p)
  if target_file == '':
    raise SystemExit  
  text = tbl2upl(dist_file, seq_file)
else:
  raise SystemExit

f = open(target_file, 'w')
f.write(text)
f.close()
print(text)