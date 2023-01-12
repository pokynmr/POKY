
#
# This is an example script to generate POKY resonance file from XEASY.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

# -----------------------------------------------------------------------------  
# POKY libraries
import __main__
s = __main__.main_session
from myseq import AAA_dict, A_dict, aaa2a, a2aaa

nuc_map = {'C': '13C', 'N': '15N', 'H': '1H', 'Q': '1H', 'M': '1H'}
# -----------------------------------------------------------------------------  
def SavePokyResonance(poky_file, external_cs_list):   
  line_queue = []
  for rec in external_cs_list:
    nseq,  a, atm, cs = rec
    try:
      nuc = nuc_map[atm[0]]
    except:
      continue
    ga = '%s%d' % (a, nseq)
    szLine = '%6s %5s %4s %7.3f   0.0     0\n' % (ga, atm, nuc, cs)
    line_queue.append(szLine)
  
  f = open(poky_file, 'w')
  f.write(''.join(line_queue))
  f.close()   

seq_file = s.open_filedialog('Select the SEQ file', 
                            'SEQ file (*.seq);; Any (*)', '')
if seq_file == '':
  raise SystemError
prot_file = s.open_filedialog('Select the XEASY prot file', 
                            'SEQ file (*.prot);; Any (*)', '')
if prot_file == '':
  raise SystemError
poky_file = s.save_filedialog('Save as', 
                  'POKY resonance list (*.list);; Any (*)', '')
if poky_file == '':
  raise SystemError

f = open(seq_file, 'r')
seq_lines = f.readlines()
f.close()

seq_list = list(map(lambda x: [x, 'XXX'], range(1000)))
  
for line in seq_lines:
  splited = line.strip().split()
  if len(splited) < 2: continue
  if splited[0][0] == '#': continue
  try:
    seq_idx = int(splited[1])
  except:
    continue
  seq_list[seq_idx][1] = splited[0]    

f = open(prot_file, 'r')
cs_lines = f.readlines()
f.close()
  
cs_list = []
for line in cs_lines:
  splited = line.strip().split()
  if len(splited) < 5: continue
  if splited[0][0] == '#': continue
  try:
    if float(splited[1]) > 500: continue
    record = [int(splited[4]), aaa2a(seq_list[int(splited[4])][1]), 
              splited[3], float(splited[1])]
    cs_list.append(record)
  except:
    continue
  
SavePokyResonance(poky_file, cs_list)
s.open_notepad(poky_file)

s.show_message('Done', 'You can read the generated file from Resonances tab.')
