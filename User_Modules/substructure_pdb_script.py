#
# This is an example script to sub-structure from a PDB file.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

input_file = s.open_filedialog('Select a PDB file', 
                    'PDB (*.pdb);; Any (*)', '')
if input_file == '':
  raise SystemExit

output_file = s.save_filedialog('Save as', 
                    'PDB (*.pdb);; Any (*)', '')

f = open(input_file, 'r')
lines = f.readlines()
f.close()

nmin, nmax = 10**9, -10**9
for line in lines:
  if line.startswith('ATOM'):
    nseq = int(line[22:26].strip())
    nmin = min(nseq, nmin)
    nmax = max(nseq, nmax)

subrng = s.show_inputdialog('Range', 'What is the substructure range?', 
                            f'{nmin}-{nmax}')
sr_list = subrng.split('-')
nmin = int(sr_list[0].strip())
nmax = int(sr_list[1].strip())
print((nmin, nmax))

if s.show_message_yes_no('Renumbering', 'Renumber to make it start from 1?'):
  offset = -1 * nmin + 1
else:
  offset = 0

content = ''
inc = 1
for line in lines:
  if line.startswith('ATOM'):
    nseq = int(line[22:26].strip())
    if nseq < nmin or nseq > nmax:
      continue
    new_line = f'{line[:5]}{inc:6d}{line[11:22]}{nseq+offset:4d}{line[26:]}'
    content += new_line
    inc += 1
  elif line.startswith('TER'):
    nseq = int(line[22:26].strip())
    new_line = f'{line[:5]}{inc:6d}{line[11:22]}{nmax+offset:4d}{line[26:]}'
    content += new_line
    inc += 1
  elif line.startswith('MODEL') or line.startswith('ENDMDL'):
    content += line


f = open(output_file, 'w')
f.write(content)
f.close()

s.show_message('Finished', 'Finished.')
