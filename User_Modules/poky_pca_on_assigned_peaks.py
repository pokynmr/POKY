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
# Only selected peaks will be considered if there are any.
# If there aren't all peaks will be considered.
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
# 0. Variance cutoff.
rho_cutoff = 0.0001
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

scale_method = s.show_selectionexdialog('Scaler', 'Scale method: ', 
                                ('Pareto', 'Unit', 'Raw'))

# If there are some selected peaks, we only used them.
# To do so, we will put the spectrum to the reference level (0)
if len(s.selected_peaks()) > 2:
  sp_list.remove(s.selected_peaks()[0].spectrum)
  sp_list = [s.selected_peaks()[0].spectrum,] + sp_list  

def preprocess(data):
  # normalize
  for i in range(len(data)):
    dmin, dmax = min(data[i]), max(data[i])
    if dmax == dmin:
      continue
    data[i] = np.divide(np.subtract(data[i], dmin), (dmax-dmin))

  for i in range(len(data[0])):
    d = data[:,i]
    # mean center
    avg = np.average(d)
    mc_data = np.subtract(d, avg)
    
    # filter
    std = np.std(mc_data)
    if std < rho_cutoff**.5:
      data[:,i] = np.multiply(mc_data, 0)
      continue

    # standard
    #mmax, mmin = np.max(mc_data), np.min(mc_data)
    #mc_data = np.divide(mc_data, (mmax - mmin))
    # apply scale
    if scale_method == 0: # pareto
      data1d = np.divide(mc_data, np.sqrt(std))
    elif scale_method == 1: # unit
      data1d = np.divide(mc_data, std)
    else: # raw
      data1d = mc_data
    data[:,i] = data1d
  
    # remove zero column  
  data = data[:,~(data==0).all(0)]
  return data
  
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
  for p in pl:
    if p.is_assigned != 1:
      continue
    if len(s.selected_peaks()) > 2 and i == 0 and p.selected == 0:
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
        try:
          res_dict[asgn[0][0]].append([p.frequency[0] * scale[0], p.frequency[1] * scale[1]])
        except:
          continue
    else:  
      if i == 0: # make dictionary
        res_list.append([''.join(asgn[0]), [p.frequency[0] * scale[0],]])
        res_list.append([''.join(asgn[1]), [p.frequency[1] * scale[1],]])
      else:
        try:
          res_dict[''.join(asgn[0])].append(p.frequency[0] * scale[0])
          res_dict[''.join(asgn[1])].append(p.frequency[1] * scale[1])
        except:
          continue

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

for i in range(len(keys)):
  res = keys[i]
  #scaled = preprocess(res_dict[res])
  #res_dict[res] = scaled
  if i == 0:
    data_stack = np.array(res_dict[res]) #scaled
  else:
    data_stack = np.vstack((data_stack, np.array(res_dict[res])))# scaled))  

# preprocess
data_stack = preprocess(data_stack)

pca = PCA(n_components=2)
converted_data = pca.fit_transform(data_stack)
loadings = pca.components_.T
print('Loadings')
print(loadings)
plt.figure()

# PC1 vs. PC2
plt.scatter(converted_data[:, 0], converted_data[:, 1])
for i in range(len(loadings)):
  l = loadings[i]
  plt.arrow(0, 0, l[0], l[1],color = 'r',alpha = 0.5)
  plt.text(l[0]* 1.15, l[1] * 1.15, 
          "Var"+str(i+1), color = 'g', ha = 'center', va = 'center')
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
