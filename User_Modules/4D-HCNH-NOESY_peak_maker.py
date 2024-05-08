#
# This is an example script to simulate 4D-HCNH NOESY peaks with assignments.
# It enables the "predict-and-confirm" assignment method for 4D-HCNH NOESY.
# Labels with "_s" mean they are from prediction (pacsy, bmrb, shiftx2)
# Labels without "_s" mean they are from existing resonances
# "_s" tags can be removed by "ut" (untag _s) or "cu" (untag _s and center) 
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#     Thomas Evangelidis's suggested edits (axis order) incorporated (05/08/24) 
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#
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
                      '6457') # this bmrb code is GB1 from Rienstra group
elif idx == 2: # SHIFTX
  reference_shift = \
    s.show_inputdialog('SHIFTX2', 
                  'Paste SHIFTX2 CSV URL path or drop the file.', 
                  'http://www.shiftx2.ca/wkdir/1611508737/1d3z.cs')
elif idx == 3: # Transfer existing resonances only
  reference_shift = 'transfer' 
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

if not (spec.nuclei == ('HC', 'C', 'N', 'HN') or \
        ''.join(sorted(spec.nuclei)) == '13C15N1H1H'):
  s.show_message('Error', 'Please choose a 4D-HCNH spectrum.')
  raise SystemExit

# determine axis order
try:
  C_idx = [i for i, n in enumerate(spec.nuclei) if 'H' not in n and 'C' in n][0]
  N_idx = [i for i, n in enumerate(spec.nuclei) if 'H' not in n and 'N' in n][0]
  H_idx = [i for i, n in enumerate(spec.nuclei) if 'H' in n][0]
except:
  s.show_message('Error', 'Please choose a 4D-HCNH spectrum.')
  raise SystemExit

h_idx = [x for x in [0, 1, 2, 3] if x not in [C_idx, N_idx, H_idx]][0]
if abs(H_idx-C_idx)+abs(h_idx-N_idx) < abs(h_idx-C_idx)+abs(H_idx-N_idx):
  HC_idx, HN_idx = H_idx, h_idx
else:
  HN_idx, HC_idx = H_idx, h_idx

pt_list = ['Intra', 'Sequential', 'Through-Space']
idx = s.show_selectionexdialog('Peak Type', 
                        'Select the peak type to simulate', 
                        tuple(pt_list + ['Cancel',]))
if idx in [len(pt_list), -1]:
  raise SystemExit
pt_mode = pt_list[idx]  

# through-space, pdb is needed
if pt_mode == 'Through-Space':
  pdb_file = s.open_filedialog('Select a PDB file', 
                    'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')
  if pdb_file == '':
    raise SystemExit
  try:
    cutoff = float(s.show_inputdialog('Cutoff', 
                                    'What is the cutoff in Angstrom?', '5'))
  except:
    s.show_message('Error', 'Please type a number.')
    raise SystemExit
  from pymol import cmd
  cmd.delete('hcnh')
  cmd.load(pdb_file, 'hcnh')

# sequence, dict: 'M1CA': 63.2
from nmrstar import fetch_bmrb_to_dict

# dict: 'M1CA': 63.2
from shiftx2 import fetch_shiftx2_to_dict

# dict: 'M1CA': 63.2
from sputil import parse_poky_shift_to_dict

# check sequence file
from myseq import ReadSequence
from atomnames import protein_atoms_by_group_no_pseudo, \
                      protein_attached_heavy_atoms, \
                      protein_pseudo_atoms

from sputil import split_group_atom, alias_onto_spectrum
from pyutil import subtract_tuples
from pacsysimulate import get_pacsy_freq
from math import sqrt

NXHX_dict, CXHX_dict = {}, {}
oneline_a = 'ARNDCQEGHILKMFPSTWYV'
# patch
protein_attached_heavy_atoms['K']['HZ'] = 'NZ'
protein_attached_heavy_atoms['K']['MZ'] = 'NZ'
for a in oneline_a:
  atm_list = protein_atoms_by_group_no_pseudo[a]
  attach_list = protein_attached_heavy_atoms[a]
  
  HN_list = list(filter(lambda x: x[0] == 'H' and attach_list[x][0] == 'N', 
                        atm_list))
  HC_list = list(filter(lambda x: x[0] == 'H' and attach_list[x][0] == 'C', 
                        atm_list))
  
  NXHX_list = list(map(lambda x: [attach_list[x], x], HN_list))
  CXHX_list = list(map(lambda x: [attach_list[x], x], HC_list))
  NXHX_dict[a] = NXHX_list
  CXHX_dict[a] = CXHX_list

seq_list = ReadSequence(s) # =[(1, 'A'), ]
if len(seq_list) == 0:
  s.show_message('Error', 'Sequence is not set properly.')
  raise SystemExit

sequence = ''
first_seqidx = seq_list[0][0]
for x in seq_list:
  sequence += x[1]

print('Spectrum Type: 4D-HCNH-NOESY')
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

def get_user_freq(g, a):
  for c in proj.condition_list():
    r = c.find_resonance(g, a)
    if r != None:
      return r.frequency
  return None

def GetShift(i, atm):
  a = sequence[i]
  if i == 0:
    ap = ''
  else:
    ap = sequence[i-1]
  if i == len(sequence) - 1:
    an = ''
  else:
    an = sequence[i+1]
  
  n = first_seqidx+i
  g = '%s%d' % (a, n)
  ga = g + atm
  
  pcs = get_pacsy_freq(a, atm, ap, an)
  cs = get_freq_from_dict(cs_dict, ga)
  ucs = get_user_freq(g, atm)

  if cs == None and pcs < 500 and idx == 0: 
    cs = pcs
  
  if ucs != None and use_resonances: 
    cs = ucs
  else: 
    ga += '_s'
  
  if reference_shift == 'transfer' and ucs == None:
    return ga, 9999
  if cs == None:
    cs = 9999
  return ga, cs

def GetDistance(i1, atm1, i2, atm2):
  try:
    atm1_list = protein_pseudo_atoms[sequence[i1-first_seqidx]][atm1]
  except:
    atm1_list = [atm1,]
  try:
    atm2_list = protein_pseudo_atoms[sequence[i2-first_seqidx]][atm2]
  except:
    atm2_list = [atm2,]

  dist = 10**10
  for at1 in atm1_list:
    for at2 in atm2_list:
      try:
        tmp = cmd.get_distance(f'(/hcnh///{i1}/{at1})', 
                          f'(/hcnh///{i2}/{at2})')
      except:
        continue
      dist = min(dist, tmp)
  return dist

# propagate NX & CX combinations
NXHX_list = list(map(lambda x: NXHX_dict[x], sequence)) # [[['N', 'H'], ...]]
CXHX_list = list(map(lambda x: CXHX_dict[x], sequence))

NXHX_dict, CXHX_dict = {}, {}
# Generate Shift dictionary for combinations
for i in range(len(sequence)):
  NXHX_combination = NXHX_list[i]
  for nx, hnx in NXHX_combination:
    nx_ga, nx_cs = GetShift(i, nx)
    hnx_ga, hnx_cs = GetShift(i, hnx)
    NXHX_dict[(i, nx, hnx)] = (nx_ga, nx_cs, hnx_ga, hnx_cs)
  
  CXHX_combination = CXHX_list[i]
  for cx, hcx in CXHX_combination:
    cx_ga, cx_cs = GetShift(i, cx)
    hcx_ga, hcx_cs = GetShift(i, hcx)
    CXHX_dict[(i, cx, hcx)] = (cx_ga, cx_cs, hcx_ga, hcx_cs)
      

if pt_mode == 'Through-Space':
  # Generate distance matrix between NXHX and CXHX
  DIST_dict = {}
  for i, nx, hnx in NXHX_dict.keys():
    nx_ga, nx_cs, hnx_ga, hnx_cs = NXHX_dict[(i, nx, hnx)]
    for j, cx, hcx in CXHX_dict.keys():
      cx_ga, cx_cs, hcx_ga, hcx_cs = CXHX_dict[(j, cx, hcx)]
      dist = GetDistance(i + first_seqidx, hnx, j + first_seqidx, hcx)
      DIST_dict[i, hnx, j, hcx] = dist

peak_list = []
for i in range(len(sequence)):
  NXHX_combination = NXHX_list[i]
  for j in range(len(sequence)):
    if pt_mode == 'Intra' and i != j:
      continue
    elif pt_mode == 'Sequential' and abs(i-j) > 1:
      continue
    
    CXHX_combination = CXHX_list[j]
    for nx, hnx in NXHX_combination:
      nx_ga, nx_cs, hnx_ga, hnx_cs = NXHX_dict[(i, nx, hnx)]
      if nx_cs > 500 or hnx_cs > 500:
        continue
   
      for cx, hcx in CXHX_combination:
        cx_ga, cx_cs, hcx_ga, hcx_cs = CXHX_dict[(j, cx, hcx)]
        if cx_cs > 500 or hcx_cs > 500:
          continue
        
        # Check distance is under cutoff if through-space mode
        if pt_mode == 'Through-Space':
          dist = DIST_dict[i, hnx, j, hcx]
          if dist > cutoff:
            continue
        
        peak_list.append( [[nx_ga, hnx_ga, cx_ga, hcx_ga], 
                           [nx_cs, hnx_cs, cx_cs, hcx_cs]] )

axis_order = [N_idx, HN_idx, C_idx, HC_idx]
for gas, freqs in peak_list:
  gas2, freqs2 = list(gas), list(freqs)
  print(f"Creating peak with label {gas2} as position {freqs2}")
  for i in range(spec.dimension):
    # Use N_idx, HN_idx, C_idx, HC_idx 
    gas2[axis_order[i]], freqs2[axis_order[i]] = \
      gas[i], freqs[i]
    
  peak = spec.place_peak(freqs2)

  # alias if edge of the spectrum
  freqs3 = alias_onto_spectrum(freqs2, spec)
  if tuple(freqs2) != freqs3:
    peak.position = freqs3
    peak.alias = subtract_tuples(freqs2, freqs3)

  # attach assignment labels
  for i in range(spec.dimension):
    g, a = split_group_atom(gas2[i])
    peak.assign(i, g, a)
  peak.show_assignment_label()

if pt_mode == 'Through-Space':
  cmd.delete('hcnh')
