#
# This is an example script to simulate nucleic acid peaks with assignments.
# It enables the "predict-and-confirm" assignment method for DNA/RNA data
# Labels with "_s" mean they are from prediction (bmrb statistics)
# Labels without "_s" mean they are from existing resonances
# "_s" tags can be removed by "ut" (untag _s) or "cu" (untag _s and center) 
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Required POKY BUILD 02/06/21d or newer
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#

# Change two lines below
first_seq = 1
sequence = 'AUGCAUGCAUGC'

import __main__
s = __main__.main_session
proj = s.project

spec_name = s.show_spectrumselectiondialog('Select a spectrum', 0)
if spec_name == '':
  raise SystemExit             
  
# use existing assignments
use_resonances = True

# DNA or RNA?
deoxy = (sequence.find('T') != 1)

# Clean sequence
seq = ''
for x in sequence:
  if 'AUGTC'.find(x) != -1:
    seq += x

# spectrum to be simulated
from sputil import name_to_spectrum
spec = name_to_spectrum(spec_name, s)
if spec == None:
  s.show_message('Error', 'Spectrum is not set properly.')
  raise SystemExit

exp_list = ['N-HSQC', 'C-HSQC']

idx = s.show_selectionexdialog('Experiment', 
                          'Select the experiment to simulate', 
                          tuple(exp_list + ['Cancel',]))
if idx in [len(exp_list), -1]:
  raise SystemExit

spec_type = exp_list[idx]

from atomnames import dna_rna_atoms_by_group, dna_rna_attached_heavy_atoms
from na_shiftstats import na_atom_statistics
from sputil import place_peak_in_spectral_region

def get_user_freq(g, a):
  for c in proj.condition_list():
    r = c.find_resonance(g, a)
    if r != None:
      return r.frequency
  return None

if idx == 0: # N-HSQC
  nuc = 'N'
elif idx == 1: # C-HSQC
  nuc = 'C'

count = 0
for i in range(len(seq)):
  seqid = i + first_seq
  x = seq[i]
  if deoxy:
    dx = 'D' + x
  else:
    dx = x
  dg = dx + str(seqid)
  g = x + str(seqid)

  group_dict = dna_rna_attached_heavy_atoms[x]
  Hlist = list(group_dict.keys())
  for hatm in Hlist:
    heavyatm = group_dict[hatm]
    if heavyatm[0] != nuc:
      continue

    Hfreq = get_user_freq(g, hatm)
    NCfreq = get_user_freq(g, heavyatm) 

    try:
      if Hfreq == None:
        Hfreq = na_atom_statistics(dx, hatm).average_shift
        hatm += '_s'
      if NCfreq == None:
        NCfreq = na_atom_statistics(dx, heavyatm).average_shift
        heavyatm += '_s'
    except:
      continue
    
    if Hfreq > 500 or NCfreq > 500:
      continue

    if spec.nuclei[0] == '1H':
      pk = place_peak_in_spectral_region(spec, [Hfreq, NCfreq])
      if pk == None:
        continue
      pk.assign(0, g, hatm)
      pk.assign(1, g, heavyatm)
    elif spec.nuclei[1] == '1H':
      pk = place_peak_in_spectral_region(spec, [NCfreq, Hfreq])
      if pk == None:
        continue
      pk.assign(0, g, heavyatm)
      pk.assign(1, g, hatm)
    else: continue # not N/C-HSQC
    
    pk.show_assignment_label() 
    count+=1

msg = '%d peaks were populated.' % (count)

s.show_message('Simulate assignments onto spectrum', msg)
