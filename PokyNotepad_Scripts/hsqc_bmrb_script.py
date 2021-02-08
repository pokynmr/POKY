#
# This is an example script to crawl BMRB entry and create N-HSQC peaks.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY BMRB HSQC crawler')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project

import sys
from sputil import name_to_spectrum
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)
if sp == None:
  s.show_message('Error', 'A N-HSQC spectrum needs to be selected.')
  sys.exit(1)

if sp.dimension != 2 or sp.nuclei not in (('15N', '1H'), ('1H', '15N')):
  s.show_message('Error', '%s does not look like a N-HSQC.' % (sp.name))
  sys.exit(1)


# BMRB entry to apply
BMRB_entry = s.show_inputdialog('BMRB entry', 'Type BMRB entry code to use.', '6457')
#BMRB_entry = '6457' # Change to your protein's BMRB entry code

# Offset to apply (ppm). Set if uniformly shifted through w1 or w2
szoffset1 = s.show_inputdialog('Offset', 'Type w1 offset to apply (ppm).', '0')
szoffset2 = s.show_inputdialog('Offset', 'Type w2 offset to apply (ppm).', '0')
offset1 = float(szoffset1)
offset2 = float(szoffset2)


# BMRB SPARKY file is H-N. Check if flipping is necessary.
if sp.nuclei == ('15N', '1H'):
  flip = True
else:
  flip = False

# download from BMRB
import requests
url = 'https://api.bmrb.io/current/entry/%s/simulate_hsqc?format=sparky&filter=all' % (BMRB_entry)
response = requests.get(url)
if response.status_code != 200:
  s.show_message('Error', 'Failed to get data from BMRB:' + response.status_code)
  sys.exit(1)

# parse and apply
n = 0
from sputil import parse_assignment
lines = response.text.split('\n')
for line in lines:
  l = line.strip().split()
  if len(l) != 3:
    continue
  try:
    asgn = parse_assignment(l[0])
  except:
    continue # failed to parse
  if asgn == None or len(asgn) != 2:
    continue
  try:
    if flip:
      p = sp.place_peak((float(l[2])+offset1, float(l[1])+offset2))
      p.assign(0, asgn[1][0], asgn[1][1])
      p.assign(1, asgn[0][0], asgn[0][1])
    else:
      p = sp.place_peak((float(l[1])+offset1, float(l[2])+offset2))
      p.assign(0, asgn[0][0], asgn[0][1])
      p.assign(1, asgn[1][0], asgn[1][1])
    p.show_assignment_label()
    n += 1
  except:
    continue

s.show_message('Info', str(n) + ' peak(s) added.')
