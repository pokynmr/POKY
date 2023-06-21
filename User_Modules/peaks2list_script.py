#
# Convert a XEASY peak list to a POKY peak list.
#
# It needs three inputs. 
#   Inputs: .seq, .prot, .peaks
#   Output: .list
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: June 21, 2023
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import __main__
s = __main__.main_session

from myseq import AAA_dict

seq_path = s.open_filedialog('Select a sequence file', 
                             'SEQ file (*.seq)', '')
if seq_path == '':
  raise SystemError

prot_path = s.open_filedialog('Select a prot file', 
                             'PROT file (*.prot)', '')
if prot_path == '':
  raise SystemError

peaks_path = s.open_filedialog('Select a peaks file', 
                             'PEAKS file (*.peaks)', '')
if peaks_path == '':
  raise SystemError

list_path = s.save_filedialog('Save a list file', 
                             'LIST file (*.list)', '')
if list_path == '':
  raise SystemError

seq_shift = int(s.show_inputdialog('Sequence shift', 
                                   'Sequence number shift', '0'))

# READ FILES
# SEQ
f = open(seq_path, 'r')
lines = f.readlines()
f.close()
seq_dict = {}
for line in lines:
  sp_list = line.strip().split()
  try:
    a = AAA_dict[sp_list[0]]
    nseq = int(sp_list[1])
    seq_dict[nseq] = a
  except:
    continue

# PROT
f = open(prot_path, 'r')
lines = f.readlines()
f.close()
prot_dict = {}
for line in lines:
  sp_list = line.strip().split()
  try:
    nidx = int(sp_list[0])
    cs = float(sp_list[1])
    atm = sp_list[3]
    nseq = int(sp_list[4])
    prot_dict[nidx] = [cs, atm, nseq]
  except:
    continue

# PEAKS
f = open(peaks_path, 'r')
lines = f.readlines()
f.close()

ndim = int(lines[0].strip().split()[-1])

if ndim == 2:
  content = 'Assignments        w1        w2\n\n'
elif ndim == 3:
  content = 'Assignments        w1        w2        w3\n\n'

for line in lines:
  sp_list = line.strip().split()
  try:
    n = int(sp_list[0])
    w1, w2 = float(sp_list[1]), float(sp_list[2])
    i1, i2 = int(sp_list[7+ndim]), int(sp_list[8+ndim])
    
    asg1, asg2, asg3 = '?', '?', '?'
    if i1 != 0:
      cs1, atm1, nseq1 = prot_dict[i1]
      a1 = seq_dict[nseq1]
      asg1 = f'{a1}{nseq1+seq_shift}{atm1}'
    if i2 != 0:
      cs2, atm2, nseq2 = prot_dict[i2]
      a2 = seq_dict[nseq2]
      asg2 = f'{a2}{nseq2+seq_shift}{atm2}'

    if ndim == 3:
      w3 = float(sp_list[3])
      i3 = int(sp_list[9+ndim])
      if i3 != 0:
        cs3, atm3, nseq3 = prot_dict[i3]
        a3 = seq_dict[nseq3]
        asg3 = f'{a3}{nseq3+seq_shift}{atm3}'
      asg = f'{asg1}-{asg2}-{asg3}'
      content += '%24s %7.3f %7.3f %7.3f\n' % (asg, w1, w2, w3)
    else:
      asg = f'{asg1}-{asg2}'
      content += '%24s %7.3f %7.3f\n' % (asg, w1, w2)
  except:
    continue


f = open(list_path, 'w')
f.write(content)
f.close()

print(content)

s.show_message('Finished.', 
    f'{list_path} has been made. Use "rp" to read the file in.')