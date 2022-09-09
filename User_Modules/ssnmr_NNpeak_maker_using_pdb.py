#
# This is an example script to simulate ssNMR peaks with assignments.
# It enables the "predict-and-confirm" assignment method for ssNMR data
# Labels with "_s" mean they are from prediction (pacsy, bmrb, shiftx2)
# Labels without "_s" mean they are from existing resonances
# "_s" tags can be removed by "ut" (untag _s) or "cu" (untag _s and center)
# 
# This version uses PDB and distances. Supports CC and NN now.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#
# Future development may include non-short mixing time data
#


import __main__
s = __main__.main_session


spec_name = s.show_spectrumselectiondialog('Select a spectrum', 0)
if spec_name == '':
  raise SystemExit             
  
idx = s.show_selectionexdialog('Reference', 
                    'Select prediction type.', 
                    ('PACSY', 'BMRB', 'SHIFTX2', 'Tranfer Only', 'Cancel'))
if idx in [4, -1]:
  raise SystemExit

if idx == 1: # BMRB
  reference_shift = \
    s.show_inputdialog('BMRB', 'Type BMRB entry number or drop the file.', 
                      '15380') # this bmrb code is GB1 from Rienstra group
elif idx == 2: # SHIFTX
  reference_shift = \
    s.show_inputdialog('SHIFTX2', 
                  'Paste SHIFTX2 CSV URL path or drop the file.', 
                  'http://www.shiftx2.ca/wkdir/1611508737/2KQ4.cs')
        # this shiftx2 url is my prediction with 2KQ4 PDB (also GB1)
elif idx == 3: # Transfer existing resonances only
  reference_shift = 'transfer' 
else: # use pacsy chemical shift database prediction (default)
  reference_shift = 'pacsy'

# pdb
pdb_file = s.open_filedialog('Select a PDB file.',
                             'PDB file (*.pdb)', '')

if pdb_file == '' or not pdb_file.endswith('.pdb'):
  print("You didn't select a PDB file")
  raise SystemExit

model_number = s.show_inputdialog('Model number',
                                  'Select the model number in PDB',
                                  '1')

# distance cutoff
ans = \
    s.show_inputdialog('Cutoff', 'Distance cutoff (e.g. H_5.5, C_7.7 or N_7.7)', 
                      'H_5.5') 
try:
  cutoffAtom = ans.split('_')[0]
  cutoffDist = float(ans.split('_')[1])
except:
  raise SystemExit

# use existing assignments
use_resonances = True

# spectrum to be simulated
from sputil import name_to_spectrum
spec = name_to_spectrum(spec_name, s)
if spec == None:
  s.show_message('Error', 'Spectrum is not set properly.')
  raise SystemExit

#avail_exp_list = ['2D-CC', '2D-NN']
#idx = s.show_selectionexdialog('Experiment', 
#                          'Select the experiment to simulate', 
#                          tuple(avail_exp_list + ['Cancel',]))
#if idx in [len(avail_exp_list), -1]:
#  raise SystemExit
#spec_type = avail_exp_list[idx]
spec_type = '2D-NN'

from myseq import AAA_dict, A_dict, ReadSequence

# sequence, dict: 'M1CA': 63.2
from nmrstar import fetch_bmrb_to_dict

# dict: 'M1CA': 63.2
from shiftx2 import fetch_shiftx2_to_dict

# dict: 'M1CA': 63.2
from sputil import parse_poky_shift_to_dict

import os, sys
import numpy as np
from math import sqrt
from pacsysimulate import get_pacsy_freq

oneline_a = 'ARNDCQEGHILKMFPSTWYV'

# check sequence file
seq_list = ReadSequence(s) # =[(1, 'A'), ]
if len(seq_list) == 0:
  s.show_message('Error', 'Sequence is not set properly.')
  raise SystemExit

sequence = ''
first_seqidx = seq_list[0][0]
for x in seq_list:
  sequence += x[1]

print('Spectrum Type: ' + spec_type)
print('Sequence: ' + sequence)

proj = s.project

if reference_shift.endswith('.list'):
  cs_dict = parse_poky_shift_to_dict(reference_shift)
elif reference_shift.endswith('.cs') or reference_shift.endswith('.csv'):
  cs_dict = fetch_shiftx2_to_dict(reference_shift)
elif reference_shift.endswith('.str') or reference_shift.endswith('.star') \
  or reference_shift.isnumeric():
  _, cs_dict = fetch_bmrb_to_dict(reference_shift)
else:
  cs_dict = None # pacsy or transfer only

# function for getting shift from dict
def get_freq_from_dict(cs_dict, ga):
  try:
    return cs_dict[ga]
  except:
    return None
# function for preparing for each dimension
def get_prepared(j, sequence, res_prof, idx):
  a = sequence[j+res_prof[idx][0]]
  if j+res_prof[idx][0] < 0: a = None
  if j+res_prof[idx][0] > len(sequence)-1: a = None
  i = first_seqidx+j+res_prof[idx][0]
  atm = atom_prof[idx]
  if j+res_prof[idx][0]-1 < 0: ap = None
  else:  ap = sequence[j+res_prof[idx][0]-1]
  if j+res_prof[idx][0]+1 > len(sequence)-1: an = None
  else: an = sequence[j+res_prof[idx][0]+1]
  return a, i, atm, ap, an

def get_user_freq(g, a):
  for c in proj.condition_list():
    r = c.find_resonance(g, a)
    if r != None:
      return r.frequency
  return None

def distance3D(atom1, atom2):
  """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                  (26.913, 26.639, -3.51))
      returns the distance
  """
  return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

def readPDB(pdb, modelnumber):
  """
  Parameters
  ----------
  pdb : address to a .pdb file
  modelnumber : the model number in the pdb (starting from 1)
  Returns
  -------
  pdb_list: a list of each [atom's coordinates, aa, atom]
  """
  pdbLines, modelList = [], []
  tempLines = open(pdb, 'r').readlines()
  
  # clean
  for line in tempLines:
    if line[0:4] in ['MODE', 'ATOM', 'ENDM']:
      pdbLines.append(line)

  # fill modelList
  for line in pdbLines:
    if line[0:5] == 'MODEL':
      modelList.append([])
    if line[0:4] != 'ATOM':
      continue
    if line[12:16].strip()[0] not in ['C', 'N', 'H', 'O']:
        continue
    aaa = line[17:20].strip()
    atm = line[12:16].strip()
    nSeq = int(line[23:26].strip())
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())

    # in case MODEL not in PDB.
    if len(modelList) == 0:
      modelList.append([])
    try:
      modelList[-1].append( [nSeq, x, y, z, aaa, AAA_dict[aaa], atm] )
    except:
      continue

  return modelList[modelnumber - 1]

def createDistanceMatrix(pdb_list, cutoffAtom, cutoffDist):
  distMat = np.zeros( (len(pdb_list), len(pdb_list)) )
  keyList = []
  for i in range(len(pdb_list)):
    nSeq, x, y, z, aaa, a, atm = pdb_list[i]
    keyList.append(a+str(nSeq)+atm)
    for j in range(i+1, len(pdb_list)):
      nSeq2, x2, y2, z2, aaa2, a2, atm2 = pdb_list[j]
      dist = distance3D( (x, y, z), (x2, y2, z2))
      distMat[i, j] = distMat[j, i] = dist
  return distMat, keyList

# 2D-NN
pt_mode = 'N'
if spec_type == '2D-NN':
  if pt_mode.find('NX') == -1:
    c_list = ['N',]
  else: 
    c_list = ['N','NH','NH1','NH2','ND','ND1','ND2','NE','NE1','NE2','NZ']
# 2D-CC- not meant to implement in this script.
elif spec_type == '2D-CC':
  c_list = ['C', 'CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 'CE', 
      'CE1', 'CE2', 'CE3', 'CH', 'CH1', 'CH2', 'CZ']

# create distance matrix
# read PDB
pdb_list = readPDB(pdb_file, int(model_number))
print(pdb_list)
distMat, keyList = createDistanceMatrix(pdb_list, cutoffAtom, cutoffDist)
peak_list = []

for i in range(len(keyList)):
  key = keyList[i]
  nSeq, x, y, z, aaa, a, atm = pdb_list[i]
  
  if atm != cutoffAtom: continue

  # dim1
  ga1 = '%s%dN' % (a, nSeq)  
  pcs1 = get_pacsy_freq(a, 'N', None, None)
  cs1 = get_freq_from_dict(cs_dict, ga1)
  ucs1 = get_user_freq(ga1[:-1], 'N')
    
  if cs1 == None and pcs1 < 500: cs1 = pcs1
  if ucs1 != None and use_resonances: cs1 = ucs1
  else: ga1 += '_s'

  if cs1 == None: cs1 = pcs1            
  if cs1 == None or cs1 > 500: continue
  
  for j in range(i+1, len(keyList)):
    key2 = keyList[j]
    nSeq2, x2, y2, z2, aaa2, a2, atm2 = pdb_list[j]
    if atm2 != cutoffAtom: continue
    if distMat[i, j] > cutoffDist: continue
            
    ga2 = '%s%dN' % (a2, nSeq2)  
    pcs2 = get_pacsy_freq(a, 'N', None, None)
    cs2 = get_freq_from_dict(cs_dict, ga2)
    ucs2 = get_user_freq(ga2[:-1], 'N')
      
    if cs2 == None and pcs2 < 500: cs2 = pcs2
    if ucs2 != None and use_resonances: cs2 = ucs2
    else: ga2 += '_s'

    if cs2 == None: cs2 = pcs2            
    if cs2 == None or cs2 > 500: continue
        
    if reference_shift == 'transfer' and None in [ucs1, ucs2]:
      continue
    if [[ga1, ga2], [cs1, cs2]] not in peak_list:           
      peak_list.append( [[ga1, ga2], [cs1, cs2]] )
      
# place assigned peaks
from sputil import split_group_atom

for peak in peak_list:
  gas, freqs = peak[0], peak[1]
  gas2, freqs2 = list(gas), list(freqs)
  peak = spec.place_peak(freqs2)
  # attach assignment labels
  for i in range(spec.dimension):
    g, a = split_group_atom(gas2[i])
    peak.assign(i, g, a)
  peak.show_assignment_label()

  # swap
  for a in range(spec.dimension):
    gas2[1-a], freqs2[1-a] = gas[a], freqs[a]
  peak2 = spec.place_peak(freqs2)
  # attach assignment labels
  for i in range(spec.dimension):
    g, a = split_group_atom(gas2[i])
    peak2.assign(i, g, a)
  peak2.show_assignment_label()
