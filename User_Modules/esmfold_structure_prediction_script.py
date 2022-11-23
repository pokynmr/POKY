#
# This is an example script for running esmfold and openning in PyMOL.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Nov. 23, 2022
#
# Reference: https://github.com/JinyuanSun/PymolFold
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# For Multimer Calculation:
#   Use Colab version from POKY Notebooks (two-letter-code PN or File Menu) 
#

######################
# USER PARAMETER START
######################
sequence = '''MQIFVKTLTG KTITLEVEPS DTIENVKAKI QDKEGIPPDQ QRLIFAGKQL
              EDGRTLSDYN IQKESTLHLV LRLRGG'''
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

def query_esmfold(seq, outname):
  seq = re.sub("[^A-Z:]", "", seq.replace("/",":").upper())
  seq = re.sub(":+",":", seq)
  seq = re.sub("^[:]+","", seq)
  seq = re.sub("[:]+$","", seq)

  headers = {
    'Content-Type': 'application/x-www-form-urlencoded',
  }

  response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                            headers=headers, data=seq)
  
  pdb_string = response.content.decode('utf-8')
  if pdb_string.startswith("HEADER"):
    with open(outname, "w") as out:
      out.write(pdb_string)
    print(f"Results saved to {outname}")
    return True
  print(pdb_string)
  return False
      
outname = s.save_filedialog('Save as', 'PDB (*.pdb);; Any (*)', '')
if outname == '':
  raise SystemError

if not query_esmfold(sequence, outname):
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
