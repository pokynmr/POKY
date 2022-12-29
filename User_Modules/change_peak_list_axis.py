#
# This is an example script for changing dimension order of a peak list.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Dec. 29, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import os

# POKY libraries
import __main__
s = __main__.main_session

from sputil import complete_assignment

infile = s.open_filedialog('Select a peak list', 
            'POKY peak list (*.list);; Any (*)', '')

if infile == '':
  raise SystemError

f = open(infile, 'r')
lines = f.readlines()
f.close()

ndim = 0
for line in lines[2:]:
  if line.find('-') != -1:
    ndim = len(line.split()[0].split('-'))
    break
if ndim == 0:
  print('The peak list looks incorrect.')
  raise SystemError

def parse_data(lines):
  data_list, skip_list = [], []
  for line in lines:
    sp_list = line.split()
    if len(sp_list) < ndim + 1:
      print('Skipped: ' + line.rstrip())
      skip_list.append(line.rstrip())
      continue
    try:
      w_list = sp_list[1:ndim+1]
      wf_list = list(map(lambda w: float(w), sp_list[1:ndim+1]))
    except:
      print('Skipped: ' + line.rstrip())
      skip_list.append(line.rstrip())
      continue
    others = ' '.join(sp_list[ndim+1:])
    a_list = complete_assignment(sp_list[0]).split('-')
    if len(a_list) != ndim:
      print('Skipped: ' + line.rstrip())
      skip_list.append(line.rstrip())
      continue
    data_list.append([a_list, w_list, wf_list, others])
  return data_list, skip_list
data_list, skip_list = parse_data(lines)

if len(data_list) == 0:
  print('The peak list does not have any valid peak line.')
  raise SystemError

outfile = s.save_filedialog('Save as', 'POKY peak list (*.list);; Any (*)',
                                  os.path.dirname(infile))
if outfile == '':
  raise SystemError

def change_axis(x1, x2, data):
  for d in data:
    a_list, w_list, wf_list, others = d
    a_list[x1], a_list[x2] = a_list[x2], a_list[x1]
    w_list[x1], w_list[x2] = w_list[x2], w_list[x1]
    wf_list[x1], wf_list[x2] = wf_list[x2], wf_list[x1]
  return data

def build_lines(data, content):
  for d in data:
    a_list, w_list, wf_list, others = d
    lbl = '-'.join(a_list)
    w = ' '.join(list(map(lambda x: '%9s' % x, w_list)))
    line = '%20s    %s    %s\n' % (lbl, w, others)
    content += line
  return content

# Ask dimensions to change.  
content = '\n'.join(skip_list) + '\n'

# 2D
if ndim == 2:
  data_list = change_axis(0, 1, data_list)
  content = build_lines(data_list, content)

# 3D
elif ndim == 3:
  mapper_3D = ('132', '213', '231', '312', '321', 'Cancel')
  mapping = s.show_selectionexdialog('Axis Order', 
    f'Select new axis order:', mapper_3D)
  if mapping == len(mapper_3D)-1:
    raise SystemError  
  new_order = mapper_3D[mapping]
  # 132: 2 <-> 3
  if new_order == '132':
    data_list = change_axis(1, 2, data_list)
  # 213: 1 <-> 2
  elif new_order == '213':
    data_list = change_axis(0, 1, data_list)
  # 321: 1 <-> 3
  elif new_order == '321':
    data_list = change_axis(0, 2, data_list)
  # 231: 1 <-> 3, 1 <-> 2
  elif new_order == '231':
    data_list = change_axis(0, 2, data_list)
    data_list = change_axis(0, 1, data_list)
  # 312: 1 <-> 3, 2 <-> 3
  elif new_order == '312':
    data_list = change_axis(0, 2, data_list)
    data_list = change_axis(1, 2, data_list)
  content = build_lines(data_list, content)

# 4D
elif ndim == 4:
  ans = s.show_inputdialog('Swap Axis',
    'Please specify two dimensions to swap (e.g. 1, 2)', '1, 2')
  dims = ans.split(',')
  try:
    dim1 = int(dims[0].strip())
    dim2 = int(dims[1].strip())
    if dim1 not in [1, 2, 3, 4] or dim2 not in [1, 2, 3, 4] or dim1 == dim2:
      print('Invalid selection of dimensions.')
      raise SystemError
    data_list = change_axis(dim1-1, dim2-1, data_list)
    content = build_lines(data_list, content)
  except:
    print('Invalid selection of dimensions.')
    raise SystemError

f = open(outfile, 'w')
f.write(content)
f.close()

s.show_message('Finished', 
  'New file is generated. You can read the file in the "rp" window.')
