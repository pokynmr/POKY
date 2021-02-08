#
# This is an example script to simulate ssNMR peaks with assignments.
# It enables the "predict-and-confirm" assignment method for ssNMR data
# Labels with "_s" mean they are from prediction (pacsy, bmrb, shiftx2)
# Labels without "_s" mean they are from existing resonances
# "_s" tags can be removed by "ut" (untag _s) or "cu" (untag _s and center) 
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
                    ('PACSY', 'BMRB', 'SHIFTX2', 'Cancel'))
if idx in [3, -1]:
  raise SystemExit

if idx == 1: # BMRB
  reference_shift = \
    s.show_inputdialog('BMRB', 'Type BMRB entry nunmber or drop the file.', 
                      '15380') # this bmrb code is GB1 from Rienstra group
elif idx == 2: # SHIFTX
  reference_shift = \
    s.show_inputdialog('SHIFTX2', 
                  'Paste SHIFTX2 CSV URL path or drop the file.', 
                  'http://www.shiftx2.ca/wkdir/1611508737/2KQ4.cs')
        # this shiftx2 url is my prediction with 2KQ4 PDB (also GB1)
else: # use pacsy chemical shift database prediction (default)
  reference_shift = 'pacsy' 

# use existing assignments
use_resonances = True

# spectrum to be simulated
from sputil import name_to_spectrum
spec = name_to_spectrum(spec_name, s)
if spec == None:
  s.show_message('Error', 'Spectrum is not set properly.')
  raise SystemExit

exp_list = ['2D-NCA', '2D-NCACB', '2D-NCO', '2D-CC', 'NCACB', 'NCACX', 
  'NCACO', 'NCOCA', 'NCOCACB', 'NCOCX', 'CANCO', 'CAN(CO)CX']
if spec.dimension == 2:
  avail_exp_list = exp_list[0:4]
elif spec.dimension == 3:
  avail_exp_list = exp_list[4:]
else:
  s.show_message('Error', '4D spectrum is not set supported.')
  raise SystemExit

idx = s.show_selectionexdialog('Experiment', 
                          'Select the experiment to simulate', 
                          tuple(avail_exp_list + ['Cancel',]))
if idx in [len(avail_exp_list), -1]:
  raise SystemExit

spec_type = avail_exp_list[idx]

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}     

# sequence, dict: 'M1CA': 63.2
from nmrstar import fetch_bmrb_to_dict

# dict: 'M1CA': 63.2
from shiftx2 import fetch_shiftx2_to_dict

# dict: 'M1CA': 63.2
from sputil import parse_poky_shift_to_dict

"""
'2D-NCA', ...     % N(i)-CA(i)
'2D-NCACB', ...   % N(i)-CA(+,i)/CB(-,i)
'2D-NCO', ...     % N(i)-CO(i-1)
'2D-CC', ...      % CX/O(i)-CX/O(i)
'NCACB', ...      % N(i)-CA(i)-CA/CB(i)
'NCACX', ...      % N(i)-CA(i)-CX/O(i)
'NCACO', ...      % N(i)-CA(i)-CO(i)
'NCOCA', ...      % N(i)-CO(i-1)-CA(i-1)
'NCOCACB', ...      % N(i)-CO(i-1)-CA/CB(i-1)
'NCOCX', ...      % N(i)-CO(i-1)-CX/O(i-1)
'CANCO', ...      % N(i)-CA(i)-CO(i-1)
'CAN(CO)CX' ...   % N(i)-CA(i)-CX/O(i-1)
"""

exp_atom_prof = [['N','CA'],  ['N','CACB'], ['N','C'],    
                ['CX','CX'], ['N','CA','CACB'], ['N','CA','CX'],  
                ['N','CA','C'],   ['N','C','CA'],     ['N','C','CACB'],   
                ['N','C','CX'],     ['N','CA','C'],    ['N','CA','CX']]
exp_res_prof =  [[[0,],[0,]], [[0,],[0,]],  [[0,],[-1,]], 
                [[0,],[0,]], [[0,],[0,],[0,]],  [[0,],[0,],[0,]], 
                [[0,],[0,],[0,]], [[0,],[-1,],[-1,]], [[0,],[-1,],[-1,]], 
                [[0,],[-1,],[-1,]], [[0,],[0,],[-1,]], [[0,],[0,],[-1,]]]

try:
  expidx = exp_list.index(spec_type)
except:
  s.show_message('Error', 'Experiment is not set properly.')
  raise SystemExit

import os, sys

oneline_a = 'ARNDCQEGHILKMFPSTWYV'

# check sequence file
from myseq import ReadSequence
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
  cs_dict = None # pacsy

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

from pacsysimulate import get_pacsy_freq

# N CA and N CO-1 based
exp = exp_list[expidx]
atom_prof = exp_atom_prof[expidx]
res_prof = exp_res_prof[expidx]
peak_list, peak_list2 = [], []
if atom_prof[0] == 'N':
  if len(atom_prof) < 3 and exp != '2D-NCACB':
    c_list = []
  elif atom_prof[2] == 'CA':
    c_list = ['CA',]
  elif atom_prof[2] == 'CACB' or exp == '2D-NCACB':  
    c_list = ['CA','CB']
  elif atom_prof[2] == 'C':
    c_list = ['C']
  else:  # CX
    c_list = ['C', 'CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 
              'CE', 'CE1', 'CE2', 'CE3', 'CH', 'CH1', 'CH2', 'CZ']
  
  for j in range(len(sequence)):
    a1, i1, atm1, a1p, a1n = get_prepared(j, sequence, res_prof, 0)
    a2, i2, atm2, a2p, a2n = get_prepared(j, sequence, res_prof, 1)
    if exp == '2D-NCACB': atm2 = 'CA'
    g1 = '%s%d' % (a1, i1)
    ga1 = '%s%d%s' % (a1, i1, atm1)
    g2 = '%s%d' % (a2, i2)
    ga2 = '%s%d%s' % (a2, i2, atm2)

    pcs1 = get_pacsy_freq(a1, atm1, a1p, a1n)
    pcs2 = get_pacsy_freq(a2, atm2, a2p, a2n)

    cs1 = get_freq_from_dict(cs_dict, ga1)
    cs2 = get_freq_from_dict(cs_dict, ga2)

    ucs1 = get_user_freq(g1, atm1)
    ucs2 = get_user_freq(g2, atm2)
    
    if cs1 == None and pcs1 < 500: cs1 = pcs1
    if cs2 == None and pcs2 < 500: cs2 = pcs2
    # user resonance assignment
    if ucs1 != None and use_resonances: cs1 = ucs1
    else: ga1 += '_s'
    if ucs2 != None and use_resonances: cs2 = ucs2
    else: ga2 += '_s'
    

    shift_list = [cs1, cs2]
    if None in shift_list or cs1 > 500 or cs2 > 500: 
      continue

    # 2D experiment
    if len(atom_prof) == 2:
      peak_list.append( [[ga1, ga2], shift_list] )

    # 3D experiment
    if len(c_list) != 0:
      if exp == '2D-NCACB':
        a3, i3, a3p, a3n = a2, i2, a2p, a2n
      else:
        a3, i3, atm3, a3p, a3n = get_prepared(j, sequence, res_prof, 2)
      for cx in c_list:
        g3 = '%s%d' % (a3, i3)
        ga3 = '%s%d%s' % (a3, i3, cx)
        pcs3 = get_pacsy_freq(a3, cx, a3p, a3n)
        cs3 = get_freq_from_dict(cs_dict, ga3)
        ucs3 = get_user_freq(g3, cx)
        if cs3 == None: cs3 = pcs3
        # user resonance assignment
        if ucs3 != None and use_resonances: cs3 = ucs3
        else: ga3 += '_s'
        if cs3 == None or cs3 > 500:
          continue
        shift_list = [cs1, cs2, cs3]
        if exp == '2D-NCACB':
          peak_list.append( [[ga1, ga3], [shift_list[0], shift_list[2]]] ) 
        else:
          peak_list.append( [[ga1, ga2, ga3], shift_list] )
else:
  # 2D-CC
  c_list = ['C', 'CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 'CE', 
        'CE1', 'CE2', 'CE3', 'CH', 'CH1', 'CH2', 'CZ']
  for j in range(len(sequence)):
    a = sequence[j]
    i1 = first_seqidx+j
    g1 = '%s%d' % (a, i1)

    if j-1 < 0: ap = None
    else: ap = sequence[j-1]
    if j+1 > len(sequence)-1: an = None
    else: an = sequence[j+1]

    # dim 1
    for cx in c_list:
      ga1 = '%s%d%s' % (a, i1, cx)  
      pcs1 = get_pacsy_freq(a, cx, ap, an)
      cs1 = get_freq_from_dict(cs_dict, ga1)
      ucs1 = get_user_freq(ga1, cx)
      
      if cs1 == None and pcs1 < 500: cs1 = pcs1
      if ucs1 != None and use_resonances: cs1 = ucs1
      else: ga1 += '_s'

      if cs1 == None: cs1 = pcs1            
      if cs1 == None or cs1 > 500: continue            
      
      # dim 2
      for cx2 in c_list:
        ga2 = '%s%d%s' % (a, i1, cx2)        
        pcs2 = get_pacsy_freq(a, cx2, ap, an)
        cs2 = get_freq_from_dict(cs_dict, ga2)
        ucs2 = get_user_freq(ga2, cx2)

        if cs2 == None and pcs2 < 500: cs2 = pcs2
        if ucs2 != None and use_resonances: cs2 = ucs2
        else: ga2 += '_s'
        
        if cs2 == None or cs2 > 500: continue           
        peak_list.append( [[ga1, ga2], [cs1, cs2]] )

# determine axis order
# In our prediction created here,
#   first dimension in peak_list is always N
#   second dimension is covalently attached carbon
#   third dimension is through-space
# We will use average and deviation to determine.
def detect_dimension(peak_list, spec):
  if len(peak_list) == 0: return None
  if len(peak_list[0][0]) != spec.dimension: return None

  dim_list, avg_list, dif_list, sorted_list, idx_list = [], [], [], [], []
  rmin, rmax = spec.region[0], spec.region[1]
  ravg, rwid = [], []

  for i in range(spec.dimension):
    dim_list.append([i, i]) # default: 0->0, 1->1, 2->2
    avg_list.append(0)
    dif_list.append([])
    sorted_list.append([])
    ravg.append( (rmin[i] + rmax[i])/2.0 )
    rwid.append( rmax[i] - rmin[i] )
  
  # the way it detects dimensions
  # 1. Use nuclei specific frequency
  # 2. Intensity at peak position- not really sure if that helps
  
  # nuclei check
  #  region = ((.4, -.5), (11.6, 10.7)) # (min-ppm, max-ppm)
  iLineCount = 0
  for peak in peak_list:
    frequency = peak[1]
    for j in range(len(frequency)):
      avg_list[j] = avg_list[j] + frequency[j]
    iLineCount += 1
  if iLineCount == 0:
    return dict(dim_list) # not a good file
    
  for i in range(spec.dimension):
    avg_list[i] = avg_list[i] / float(iLineCount)
    
  # just find the best.
  # digits
  for i in range(spec.dimension):
    # dimensions
    for j in range(spec.dimension):
      dif_list[i].append(abs(avg_list[i]-ravg[j])/rwid[j])
    sorted_list[i] = sorted_list[i] + sorted(dif_list[i])
    idx_list.append(dif_list[i].index(sorted_list[i][0]))
  if len(list(set(idx_list))) != len(idx_list): 
    return dict(dim_list)    # not distinguishable
  
  for i in range(spec.dimension):
    dim_list[i][1] = idx_list[i] # distinguishable
  return dict(dim_list)

# place assigned peaks
from sputil import split_group_atom, alias_onto_spectrum
from pyutil import subtract_tuples
dim_dict = detect_dimension(peak_list, spec)
if dim_dict == None:
  s.show_message('Error', 
            'Spectrum does not look the correct experiment to use.')
  raise SystemExit

for peak in peak_list:
  gas, freqs = peak[0], peak[1]
  gas2, freqs2 = list(gas), list(freqs)
  for a in range(spec.dimension):
    b = dim_dict[a]
    gas2[b], freqs2[b] = gas[a], freqs[a]
  peak = spec.place_peak(freqs2)

  # alias if edge of the spectrum
  freqs3 = alias_onto_spectrum(freqs2, spec)
  if freqs2 != freqs3:
    peak.position = freqs3
    peak.alias = subtract_tuples(freqs2, freqs3)

  # attach assignment labels
  for i in range(spec.dimension):
    g, a = split_group_atom(gas2[i])
    peak.assign(i, g, a)
  peak.show_assignment_label()
