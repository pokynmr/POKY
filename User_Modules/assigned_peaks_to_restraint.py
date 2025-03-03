#
# This is an example script for creating restraints from assigned peaks.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Mar. 3, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

# POKY libraries
import __main__
s = __main__.main_session
from sputil import name_to_spectrum

specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

new_path = s.save_filedialog('Save as', 'UPL file (*.upl);; Any (*)', '')
if new_path == '':
  raise SystemError

from myseq import A_dict
peaks = sp.peak_list()
assigned_peaks = list(filter(lambda x: x.is_assigned, peaks))
heights = list(map(lambda x: abs(x.data_height), assigned_peaks))
heights.sort(reverse=True)

cutoffs = [heights[int(len(heights)/3)], 
           heights[int(2*len(heights)/3)], 
           heights[-1]]
distances = [4.9, 6.4, 7.7]

ans = s.show_inputdialog('Height cutoffs for bins',
  f'Specify height cutoffs for bins (e.g. {cutoffs[0]:.2f}, ' + \
  f'{cutoffs[1]:.2f}, {cutoffs[2]:.2f})',
  f'{cutoffs[0]:.2f}, {cutoffs[1]:.2f}, {cutoffs[2]:.2f}')
if ans == '':
  raise SystemError
ans_list = ans.split(',')
cutoffs = list(map(lambda x: float(x.strip()), ans_list))

ans = s.show_inputdialog('Distances for bins',
        f'Specify upper limit distances for bins (e.g. {distances[0]:.2f}, ' + \
        f'{distances[1]:.2f}, {distances[2]:.2f})',
        f'{distances[0]:.2f}, {distances[1]:.2f}, {distances[2]:.2f}')
if ans == '':
  raise SystemError
ans_list = ans.split(',')
distances = list(map(lambda x: float(x.strip()), ans_list))

if sp.dimension == 2:
  i1, i2 = 0, 1
else:
  if sp.dimension == 3:
    dim_list = [[0, 1], [0, 2], [1, 2]]
  elif sp.dimension == 4:
    dim_list = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
  dim_list2 = list(map(lambda x: f'w{x[0]+1},w{x[1]+1}', dim_list))
  ans = s.show_selectionexdialog('Select dimensions to use', 
                                 dim_list2 + ['Cancel',])
  if ans in [-1, len(ans)-1]:
    raise SystemError
  i1, i2 = dim_list[ans]

content, nupl = '', 0
for peak in assigned_peaks:
  res = peak.resonances()
  if res[i1] == res[i2]: # diagonal skipping
    continue
  for i, cutoff in enumerate(cutoffs):
    if abs(peak.data_height) > cutoff:
      line = '%6d %4s %4s %6d %4s %4s %8.2f\n' % (res[i1].group.number, 
          A_dict[res[i1].group.name[0]], res[i1].atom.name, res[i2].group.number, 
          A_dict[res[i2].group.name[0]], res[i2].atom.name, distances[i])
      content += line
      nupl += 1
      break

if content == '':
  s.show_message('Error', 'No available restraints to create.')
  raise SystemError

f = open(new_path, 'w')
f.write(content)
f.close()

s.show_message('Finished.', 
               f'{nupl} restraints saved to {new_path}')