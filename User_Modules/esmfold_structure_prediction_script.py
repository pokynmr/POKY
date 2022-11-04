#
# This is an example script for running esmfold and open in PyMOL.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Nov. 3, 2022
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
sequence = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
pymol = True
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
  s.set_clipboard(f'load {outname}')
  s.show_message('PyMOL', 'Ctrl+V into PyMOL commandline.')
  s.open_pymol('')