#
# This is an example script to count short-, medium- and long-range restraints.
# Each line should represent a single restraint!
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.session

distfile = s.open_filedialog('Select distance restraint file.', 
    'DIANA UPL (*.upl);; XPLOR TBL (*.tbl);; Any (*)', '')

if distfile == '':
  raise SystemExit

f = open(distfile, 'r')
lines = f.readlines()
f.close()

shortrange = 0
mediumrange = 0
longrange = 0

for line in lines:
  temp = line.strip()
  temp = temp.replace('(', ' ')
  temp = temp.replace(')', ' ')
  if temp == '':
    continue
  if temp[0] in ['!', '#', '*', '/']:
    continue

  sp = temp.split()
  idx_list = []

  for seg in sp:
    try:
      nseg = int(seg)
      idx_list.append(nseg)
    except:
      pass
    if len(idx_list) == 2:
      if abs(idx_list[0]-idx_list[1]) <= 1:
        shortrange += 1
      elif abs(idx_list[0]-idx_list[1]) < 5:
        mediumrange += 1
      else:
        longrange += 1
      break
msg = 'Short range |i-j| <= 1 : ' + str(shortrange)
msg = msg + '\nMedium range 1 < |i-j| < 5 : ' + str(mediumrange)
msg = msg + '\nLong range |i-j| >= 5 : ' + str(longrange)  
s.show_message('Distance Restraints', msg)