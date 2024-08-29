#
# This is an example script for simulate a spectrum using sequence.
# It basically integrates ESMFold and UCBShift or SHIFTX2.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Apr. 7, 2023
#
# Reference:  https://github.com/JinyuanSun/PymolFold
#             https://github.com/THGLab/CSpred
#             https://shiftx2.ca
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# In POKY CSPred Client, click Check the results and Simulate butons. 
#

######################
# USER PARAMETER START
######################
sequence = '''MQIFVKTLTG KTITLEVEPS DTIENVKAKI QDKEGIPPDQ QRLIFAGKQL
              EDGRTLSDYN IQKESTLHLV LRLRGG'''
method = 'UCBShift' # or 'SHIFTX2'. SHIFTX2 requires POKY 021323e or newer
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
                            headers=headers, data=seq, verify=False)
  
  pdb_string = response.content.decode('utf-8')
  if pdb_string.startswith("HEADER"):
    with open(outname, "w") as out:
      out.write(pdb_string)
    print(f"Results saved to {outname}")
    return True
  print(pdb_string)
  return False

import random
import tempfile
import time
base = tempfile.gettempdir()
random.seed()
token = time.strftime("%y%m%d_%H%M%S_") + \
            "%06d" % random.randrange(0, 1000000)
tmp_outname = os.path.join(base, token + '.pdb')

if not query_esmfold(sequence, tmp_outname):
  s.show_message('Error', 'ESMFold failed. Check your sequence.')
  raise SystemError

import Poky_CSpred
try:
  if method == 'UCBShift':
    d = Poky_CSpred.show_CSpred(s)
  elif method == 'SHIFTX2':
    d = Poky_CSpred.show_SHIFTXpred(s)
  d.pdb_field.set(tmp_outname)
except:
  s.show_message('Error', 'Update POKY. Your version is old.')
  raise SystemError

if not d.submit_to_server():
  s.show_message('Error', 'POKY CSPred failed. Try later.')
  raise SystemError
 
s.show_message('Submitted', 
  'Check the results and Simulate once done in the Poky CSpred Client.')
