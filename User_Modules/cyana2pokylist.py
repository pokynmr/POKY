#
# Convert a Cyana peak list (.peaks) together with a Cyana proton file
# (.prot) and a sequence file (.seq) into a POKY peak list file (.list).
#
# This is a Notepad script counterpart of the "cy" / "Read Cyana peak list"
# module (modules/poky/readcyana.py), but instead of placing peaks onto a
# spectrum it writes an assigned POKY space-delimited peak list that can be
# read back with the "rp" command.
#
#   Inputs: .seq, .prot, .peaks
#   Output: .list
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: September 24, 2026
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import __main__
s = __main__.main_session

from myseq import aaa2a
from cyana import CyanaDictionary

amino = {'ALA' : 'A',  # alanine
   'CYS' : 'C',  # cysteine
   'CYSS' : 'C',  # cysteine
   'ASP' : 'D',  # aspartate
   'ASP-' : 'D',  # aspartate
   'GLU' : 'E',  # glutamate
   'GLU-' : 'E',  # glutamate
   'PHE' : 'F',  # phenylalanine
   'GLY' : 'G',  # glycine
   'HIS' : 'H',  # histidine
   'HIS+' : 'H',  # histidine
   'HEM' : 'H',  # histidine
   'ILE' : 'I',  # isoleucine
   'LYS' : 'K',  # lysine
   'LYS+' : 'K',  # lysine
   'LEU' : 'L',  # leucine
   'ASN' : 'N',  # asparagine
   'MET' : 'M',  # methionine
   'PRO' : 'P',  # proline
   'GLN' : 'Q',  # glutamine
   'ARG' : 'R',  # arginine
   'ARG+' : 'R',  # arginine
   'SER' : 'S',  # serine
   'THR' : 'T',  # threonine
   'VAL' : 'V',  # valine
   'TRP' : 'W',  # tryptophan
   'TYR' : 'Y'}  # tyrosine

seq_path = s.open_filedialog('Select a Cyana sequence file',
                             'SEQ file (*.seq);; Any (*)', '')
if seq_path == '':
  raise SystemError

prot_path = s.open_filedialog('Select a Cyana prot file (optional)',
                             'PROT file (*.prot);; Any (*)', '')
if prot_path == '':
  prot_path = None

peaks_path = s.open_filedialog('Select a Cyana peaks file',
                             'PEAKS file (*.peaks);; Any (*)', '')
if peaks_path == '':
  raise SystemError

list_path = s.save_filedialog('Save a POKY peak list file',
                             'LIST file (*.list);; Any (*)', '')
if list_path == '':
  raise SystemError

sel = s.show_selectiondialog('Peak list dimension',
                             'Select the dimensions of the peak list:',
                             ('2D', '3D', 'Cancel'))
if sel == 2:
  raise SystemError
is3d = (sel == 1)

seq_shift = int(s.show_inputdialog('Sequence shift',
                                   'Sequence number shift', '0'))

# READ SEQ
# seq[nseq] = 3-letter code, keyed by integer residue number.
f = open(seq_path, 'r')
slines = f.readlines()
f.close()

seq = {}
next_idx = 1
for line in slines:
  sf = line.lstrip().split()
  if len(sf) < 1: continue
  if sf[0][0] == '#': continue
  if len(sf) > 1:
    try:
      nseq = int(sf[1])
      aa3 = sf[0]
    except:
      try:
        nseq = int(sf[0])
        aa3 = sf[1]
      except:
        continue
  else:
    aa3 = sf[0]
    nseq = next_idx
  next_idx = nseq + 1
  seq[nseq] = aa3

# READ PROT (optional)
prot_dict = {}
if prot_path is not None:
  f = open(prot_path, 'r')
  plines = f.readlines()
  f.close()
  for line in plines:
    pfields = line.split()
    if len(pfields) < 5: continue
    if pfields[0][0] == '#': continue
    try:
      nidx = int(pfields[0])
      atm = pfields[3]
      nseq = int(pfields[4])
      prot_dict[nidx] = [atm, nseq]
    except:
      continue

cyana_dict = CyanaDictionary()

def resolve(nn_entry):
  # Return (aa1, nseq, atom) for a Cyana assignment entry, or None.
  if '.' in nn_entry:
    parts = nn_entry.split('.')
    if len(parts) != 2: return None
    try:
      rns = int(parts[0]); atm = parts[1]
    except:
      try:
        atm = parts[0]; rns = int(parts[1])
      except:
        return None
  else:
    try:
      nidx = int(nn_entry)
    except:
      return None
    if nidx == 0: return None
    if nidx not in prot_dict: return None
    atm, rns = prot_dict[nidx]
    rns = int(rns)

  if rns not in seq: return None
  try:
    aa3 = seq[rns]
  except:
    return None
  try:
    aa1 = amino[aa3]
  except:
    aa1 = aaa2a(aa3)
  if aa1 == 'X': return None

  # Cyana -> IUPAC nomenclature
  atm = cyana_dict.toIUPAC(aa1, atm)
  return (aa1, rns, atm)

# READ PEAKS
f = open(peaks_path, 'r')
splines = f.readlines()
f.close()

if is3d:
  header = '      Assignment         w1         w2         w3   Data Height\n\n'
else:
  header = '      Assignment         w1         w2   Data Height\n\n'
content = header

for line in splines:
  ls = line.lstrip()
  if len(ls) == 0: continue
  if ls[0] == '#': continue

  fields = line.split()
  try:
    if is3d:
      if len(fields) < 12: continue
      x = float(fields[1]); y = float(fields[2]); z = float(fields[3])
      integ = fields[6] if len(fields) > 6 else ''
      if len(fields) >= 13:
        nn = [fields[10], fields[11], fields[12]]
      else:
        nn = [fields[0], fields[1], fields[2]]
      freqs = [x, y, z]
    else:
      if len(fields) < 10: continue
      x = float(fields[1]); y = float(fields[2])
      integ = fields[5] if len(fields) > 5 else ''
      if len(fields) >= 11:
        nn = [fields[9], fields[10]]
      else:
        nn = [fields[0], fields[1]]
      freqs = [x, y]
  except:
    continue

  asg_parts = []
  for entry in nn:
    res = resolve(entry)
    if res is None:
      asg_parts.append('?')
    else:
      aa1, rns, atm = res
      asg_parts.append('%s%d%s' % (aa1, rns + seq_shift, atm))
  asg = '-'.join(asg_parts)

  try:
    dh = float(integ)
  except:
    dh = 0.0

  line_out = '%20s' % asg + ''.join('%10.3f' % w for w in freqs) + '%12.0f' % dh
  content += line_out + '\n'

f = open(list_path, 'w')
f.write(content)
f.close()

print(content)
s.show_message('Finished.',
    f'{list_path} has been made. Use "rp" to read the file in.')
