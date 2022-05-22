#
# POKY PCA on Raw Spectra (POKY PORS)
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Last Update: May 21, 2022 
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# Warning! This can be very slow!
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY PCA on Raw Spectra (POKY PORS)')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

######################################################################
# PARAMETERS
# 0. Variance cutoff
rho_cutoff = 0.6
# 1. Use GUI to select spectra, and plot PC1 vs. PC2.
spec_names = s.show_spectrumselectiondialog('Select spectra', 1)
specname_list = spec_names.split('\t')

# 2. Use below, and plot user define X vs. PC2. 
#spec_list = [ # spec name,    X value
#            ['protein_ligand_1_0',  0],
#            ['protein_ligand_1_1',  1],
#            ['protein_ligand_1_2',  2],
#            ['protein_ligand_1_4',  4],
#            ['protein_ligand_1_8',  8],
#            ]
#x_unit = 'Molar Ratio'
#specname_list = list(map(lambda x: x[0], spec_list))
#x_list = list(map(lambda x: x[1], spec_list))
######################################################################

from sputil import name_to_spectrum
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

ref_peak = None
if len(s.selected_peaks()) == 1:
  if s.show_message_yes_no('Peak Selected', 
      'Is this selected peak the reference standard?'):
    ref_peak = s.selected_peaks()[0]


def preprocess(data):
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

for i in range(0, len(sp_list)):
  sp = sp_list[i]
  dic, data = ng.sparky.read_lowmem(sp.data_path)
  
  if i == 0:
    data_stack = np.array(data).flatten()
  else:
    data_stack = np.vstack((data_stack, np.array(data).flatten()))

# apply reference if exists
  if ref_peak != None:
    ref_hts = sp.data_height(ref_peak.frequency)
    data = np.divide(data, ref_hts)
  
print(data_stack.shape)
data_stack = preprocess(data_stack)
print(data_stack.shape)

pca = PCA(n_components=2)
converted_data = pca.fit_transform(data_stack)

plt.figure()

# PC1 vs. PC2
if 'spec_list' not in locals():
  plt.scatter(converted_data[:, 0], converted_data[:, 1])
  plt.xlabel('PC1') , plt.ylabel('PC2')
  print('----------------------------------------------')
  print('%-20s %-15s %-15s' % \
    ('SPECTRUM', 'PC1', 'PC2'))
  print('----------------------------------------------')
  for i in range(len(specname_list)):
    plt.annotate(specname_list[i], converted_data[i])
    print('%20s %15.3f %15.3f' % \
      (specname_list[i], converted_data[i,0], converted_data[i,1]))

# X vs. PC1
else:
  plt.scatter(x_list, converted_data[:, 0])
  plt.xlabel(x_unit) , plt.ylabel('PC1')
  print('-------------------------------------------------------------')
  print('%-20s %-15s %-15s %-15s' % \
    ('SPECTRUM', x_unit, 'PC1', 'PC2'))
  print('-------------------------------------------------------------')
  for i in range(len(specname_list)):
    plt.annotate(specname_list[i], (x_list[i], converted_data[i,0]))
    print('%20s %15.3f %15.3f %15.3f' % \
      (specname_list[i], x_list[i], converted_data[i,0], converted_data[i,1]))

plt.show(block=False)
