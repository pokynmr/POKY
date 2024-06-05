#
# This is an example script for running esmfold with a dimer.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: June 5, 2024
#
# Reference: https://github.com/JinyuanSun/PymolFold
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

######################
# USER PARAMETER START
######################
sequence_1 = '''MDSNTVSSFQVDCFLWHVRK RFADQELGDAPFLDRLRRDQ KSLRGRGSTLGLDIETATRA 
                GKQIVERILKEES'''
sequence_2 = '''MDSNTVSSFQVDCFLWHVRK RFADQELGDAPFLDRLRRDQ KSLRGRGSTLGLDIETATRA 
                GKQIVERILKEES'''
linker_size = 25  # the number of artificial glycine residues
pymol = True
pLDDTcolor = True
######################
# USER PARAMETER END
######################

import requests
import re
import os

# POKY libraries
import __main__
s = __main__.main_session

def clean_seq(seq):
  seq = re.sub("[^A-Z:]", "", seq.replace("/",":").upper())
  seq = re.sub(":+",":", seq)
  seq = re.sub("^[:]+","", seq)
  seq = re.sub("[:]+$","", seq)
  return seq

def query_esmfold(seq1, seq2, linker, outname):
  seq1 = clean_seq(seq1)
  seq2 = clean_seq(seq2)
  seq = seq1 + linker + seq2

  headers = {
    'Content-Type': 'application/x-www-form-urlencoded',
  }

  response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                            headers=headers, data=seq, verify=False)
  new_pdb_string = ''
  pdb_string = response.content.decode('utf-8')
  if pdb_string.startswith("HEADER"):
    # remove artificial linker
    pdb_lines = pdb_string.split('\n')
    for line in pdb_lines:
      if not line.startswith('ATOM '):
        new_pdb_string += line + '\n'
        continue
      nseq = int(line[23:26].strip())
      if nseq <= len(seq1):               # CHAIN A
        new_pdb_string += line + '\n'
        continue
      if nseq > len(seq1) + len(linker):  # CHAIN B
        new_pdb_string += '%s B %3d%s\n' % (line[:20], 
                                       nseq - len(linker) - len(seq1),
                                       line[26:])
        continue
      # linker skipped
    '''
ATOM    336  HA  ASP A  21
    '''    
    with open(outname, "w") as out:
      out.write(new_pdb_string)
    print(f"Results saved to {outname}")
    return True
  if new_pdb_string == '':
    print(pdb_string)
    print('Linker removal failed.')
  else:
    print(new_pdb_string)
  return False
      
outname = s.save_filedialog('Save as', 'PDB (*.pdb);; Any (*)', '')
if outname == '':
  raise SystemError

linker = 'G' * linker_size
if not query_esmfold(sequence_1, sequence_2, linker, outname):
  raise SystemError

if pymol:
  objname = os.path.splitext(os.path.basename(outname))[0]  
  c_list = ['0x0053d7', '0x57caf9', '0xffdb12', '0xff7e45']
  cmd = f'load {outname}; '
  if pLDDTcolor:
    cmd += f'color {c_list[0]}, ({objname}) and (b >0.90 or b =0.90); '
    cmd += f'color {c_list[1]}, ({objname}) and ((b <0.90 and b >0.70) or (b =0.70)); '
    cmd += f'color {c_list[2]}, ({objname}) and ((b <0.70 and b >0.50) or (b =0.50)); '
    cmd += f'color {c_list[2]}, ({objname}) and ((b <0.50 and b >0.0 ) or (b =0.0))'
  s.set_clipboard(f'{cmd}')
  s.show_message('PyMOL', 'Ctrl+V into PyMOL commandline.')
  s.open_pymol('')
