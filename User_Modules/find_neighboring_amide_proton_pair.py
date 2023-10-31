#
# Find a neighboring inter-amide proton pair.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jan. 10, 2023
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

######################
# USER PARAMETER START
######################
sequence = '''MQIFVKTLTG KTITLEVEPS DTIENVKAKI QDKEGIPPDQ QRLIFAGKQL
              EDGRTLSDYN IQKESTLHLV LRLRGG'''
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

from pymol import cmd
cmd.delete('dt_neighbor')
cmd.load(tmp_outname, 'dt_neighbor')
nmodel = cmd.count_states('dt_neighbor')
cmd.do('h_add dt_neighbor & name N')

# get min, max residue number
nmin, nmax = 10**5, -10**5
atom_list = cmd.get_model("dt_neighbor & name CA").atom
for atm in atom_list:
  nmin = min(int(atm.resi), nmin)
  nmax = max(int(atm.resi), nmax)

import numpy as np
h_dist = np.zeros((nmax-nmin+1, nmax-nmin+1))

# calculate inter-amide-protein distance
for i in range(nmin, nmax):
  xyz = cmd.get_coords(f'dt_neighbor & i. {i} & n. H*')
  for j in range(i+1, nmax+1):
    xyz2 = cmd.get_coords(f'dt_neighbor & i. {j} & n. H*')
    try:
      dist = ((xyz[0][0]-xyz2[0][0])**2 + (xyz[0][1]-xyz2[0][1])**2 + \
        (xyz[0][2]-xyz2[0][2])**2)**.5
    except:
      continue
    h_dist[i-nmin, j-nmin] = h_dist[j-nmin, i-nmin] = dist

msg = 'Select the contact distance cutoff in Angstrom (e.g. 5)'
cutoff = float(s.show_inputdialog('Distance cutoff', msg, '5'))

cmd.delete('dt_neighbor')

if cutoff == '':
  raise SystemError

while True:
  msg = 'Type the residue number that you want to find neighbors.\n' +\
      ' Cancel if you want to stop finding.'
  nseq = int(s.show_inputdialog('Residue number', msg, str(nmin)))
  if msg == '':
    break  
  neighbor = ''
  for i in range(nmin, nmax):
    if i == nseq:
      continue
    if h_dist[nseq-nmin, i-nmin] <= cutoff and h_dist[nseq-nmin, i-nmin] > 1:
      neighbor = neighbor + str(i) + ', '
  if neighbor == '':
    s.show_message('Done', 'No neighboring amide proton found.')
  else:
    s.show_message('Done', 'Neighboring amide protons: \n' + neighbor[:-2])
