#
# POKY PCA on Assigned Peaks (POKY POAP)
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Last Update: May 21, 2022 
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY PCA on Assigned Peaks (POKY POAP)')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

######################################################################
# PARAMETERS
# 1. Use GUI to select spectra, and plot PC1 vs. PC2.
spec_names = s.show_spectrumselectiondialog('Select spectra', 1)
specname_list = spec_names.split('\t')

# Adjusting gyromagnetic ratio
NtoH = 0.2
CtoH = 0.25
######################################################################

# Make a pseudo chemical shift from two different dimensions
CombineNuclei = s.show_message_yes_no('Nuclei consideration', 
  'Do you want to consider two dimensions as one?')

from sputil import name_to_spectrum, parse_assignment
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt 
from sklearn.decomposition import PCA

from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

sp_list = list(map(lambda sp: name_to_spectrum(sp, s), specname_list))
if len(sp_list) < 3:
  print('At least use more than two spectra.')
  raise SystemExit

if None in sp_list:
  none_name = sp_list[sp_list.index(None)][0]
  print('Experiment (' + none_name + ') does not exist in your spectrum list.')
  raise SystemExit

def preprocess(data1d):
  # mean center
  avg = np.average(data1d)
  mc_data = np.subtract(data1d, avg)

  # pareto scale
  std = np.std(mc_data)
  data1d = np.divide(mc_data, np.sqrt(std))
  return data1d

res_list = []
for i in range(len(sp_list)):
  sp = sp_list[i]
  pl = sp.peak_list()
  scale = [1, 1]
  if '15N' in sp.nuclei:
    scale[sp.nuclei.index('15N')] = NtoH
  if '13C' in sp.nuclei:
    scale[sp.nuclei.index('13C')] = CtoH
  if i == 1:
    res_dict = dict(res_list)
    print(res_dict)
  for p in pl:
    if p.is_assigned != 1:
      continue
    # not assigned in all spectra definitely.
    if p.resonances()[0].peak_count < len(sp_list):
      continue
    if p.resonances()[1].peak_count < len(sp_list):
      continue
    asgn = parse_assignment(p.assignment)
    if asgn == None:
      continue

    if CombineNuclei:
      if i == 0: # make dictionary
        res_list.append([asgn[0][0], [[p.frequency[0] * scale[0], p.frequency[1] * scale[1]],]])
      else:
        res_dict[asgn[0][0]].append([p.frequency[0] * scale[0], p.frequency[1] * scale[1]])
    else:  
      if i == 0: # make dictionary
        res_list.append([''.join(asgn[0]), [p.frequency[0] * scale[0],]])
        res_list.append([''.join(asgn[1]), [p.frequency[1] * scale[1],]])
      else:
        res_dict[''.join(asgn[0])].append(p.frequency[0] * scale[0])
        res_dict[''.join(asgn[1])].append(p.frequency[1] * scale[1])

keys = list(res_dict.keys())

if CombineNuclei:
  # measure euclidean distances from the very first one.
  for i in range(len(keys)):
    res = keys[i]
    data_list = res_dict[res]
    dist = [0, ]
    refc1, refc2 = data_list[0] # reference point
    for j in range(1, len(data_list)):
      c1, c2 = data_list[j]
      dist.append(((refc1-c1)**2 + (refc2-c2)**2)**.5)
    res_dict[res] = dist # re-assign with distances

# preprocess
for i in range(len(keys)):
  res = keys[i]
  scaled = preprocess(res_dict[res])
  #res_dict[res] = scaled
  if i == 0:
    data_stack = scaled
  else:
    data_stack = np.vstack((data_stack, scaled))  

pca = PCA(n_components=2)
converted_data = pca.fit_transform(data_stack)

plt.figure()

# PC1 vs. PC2
plt.scatter(converted_data[:, 0], converted_data[:, 1])
plt.xlabel('PC1') , plt.ylabel('PC2')
print('----------------------------------------------')
print('%-20s %-15s %-15s' % \
  ('RESONANCE', 'PC1', 'PC2'))
print('----------------------------------------------')
for i in range(len(keys)):
  plt.annotate(keys[i], converted_data[i])
  print('%20s %15.3f %15.3f' % \
    (keys[i], converted_data[i,0], converted_data[i,1]))

plt.show(block=False)
