#
# This is an example script to convert NMRVIEW XPK to POKY/SPARKY LIST.
# However, this does not make assignment labels. Only peaks!
# That is because I do not know NMRVIEW at all.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.session

xpkfile = s.open_filedialog('Select NMRVIEW XPK file.', 
    'NMRVIEW XPK (*.xpk);; Any (*)', '')
if xpkfile == '':
  raise SystemExit    

import os.path
listfile = s.save_filedialog('Save as...', 
    'POKY LIST (*.list);; Any (*)', os.path.dirname(xpkfile))
if listfile == '':
  raise SystemExit    
  
f = open(xpkfile, 'r')
lines = f.readlines()
f.close()

# check if nmrview file
# if not, just print same file to let it be used as it is
if len(lines) < 10:
  for line in lines:
    print(line.rstrip())
  raise SystemExit
    
iXpk = 0
for i in range(5):
  line = lines[i]
  if line.find('.nv') != -1 or line.find('.xpk') != -1 or \
      line.find('label dataset') != -1:
    iXpk = 1
    break

if iXpk == 0:
  for line in lines:
    print(line.rstrip())
  raise SystemExit

# find column description line
iDesc = -1
for i in range(10):
  line = lines[i]
  slist = line.strip().split()
  if len(slist) < 8: continue
  # we expect .P to be in.
  if line.find('.P') != -1:
    iDesc = i
    break

if iDesc == -1:    
  for line in lines:
    print(line.rstrip())
  raise SystemExit
    
# check columns
col_list = []
int_clm = -1
for i in range(len(slist)):
  if slist[i].find('.P') != -1:
    col_list.append(i+1)
    continue
  if slist[i] == 'int':
    int_clm = i+1
    continue

if len(col_list) == 2:
  header = '    %10s     %5s     %5s    ' % ('Assignment', 'w1', 'w2')
elif len(col_list) == 3:
  header = '    %10s     %5s     %5s     %5s    ' % \
            ('Assignment', 'w1', 'w2', 'w3')
elif len(col_list) == 4:
  header = '    %10s     %5s     %5s     %5s     %5s    ' % \
            ('Assignment', 'w1', 'w2', 'w3', 'w4')
else:    
  for line in lines:
    print(line.rstrip())
  sys.exit()

text = ''
scaled = 1
if int_clm != -1:
  header = header + '   Data Height'
  iHigh = 0
  try:
    for i in range(iDesc+1, min(20, len(lines))):
      slist = lines[i].strip().split()    
      if abs(float(slist[int_clm])) > 10000:
        iHigh = iHigh + 1
  except:
    pass        
text += header
text += '\n\n'
    
if iHigh > 5: scaled = 0
    
for i in range(iDesc+1, len(lines)):
  slist = lines[i].strip().split()
  if len(slist) < 10: continue
  pline = '%14s' % ('?-?-?')
  
  try:
    for clm in col_list:  
      pline = pline + ' %9.3f' % (float(slist[clm]))
    if int_clm != -1:
      if scaled == 1:
        pline = pline + ' %17.3f' % (float(slist[int_clm]) * 1000000.0)
      else:
        pline = pline + ' %17.3f' % (float(slist[int_clm]))      
    text += pline
    text += '\n'
  except:
    continue

print(text)
f = open(listfile, 'w')
f.write(text)
f.close()
