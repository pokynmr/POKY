#
# This is an example script to convert sequence formats.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

infile = '/path/to/your/sequence/file.seq'
outfile = '/path/to/save/your/sequence/file.fasta'

# three2one or one2three
conversion_mode = 'three2one'

f = open(infile, 'r')
lines = f.readlines()
f.close()

if mode == 'three2one':
  fasta = ''
  for line in lines:
    a = AAA_dict[line.strip().split()[0]]
    fasta += a

  print(fasta)
  f = open(outfile, 'w')
  f.write(fasta)
  f.close()
elif mode == 'one2three':
  f = open(infile, 'w')
  idx = 1
  for line in lines:
    if line.strip() == '' or '#!>'.find(line.strip()[0]) == -1:
      continue
    for a in line:
      try:
        aaa = A_dict[a]
        out = '%4s %d' % (aaa, idx)
        print(out)
        f.write(out+'\n')
        idx += 1
      except:
        pass
  f.close()
