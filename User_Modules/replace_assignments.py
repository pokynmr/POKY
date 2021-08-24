#
# Replace assignments 
# 
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: August 24, 2021
#
# Usage:
#    1. Define mapping dictionary.
#    2. File > Run Python Module
#

# User defined mapping dictionary.
# Need to change.
map_dict = {
  # For example...
  'M1CA': 'M2CB', # specific atom from specific amino acid
  'V4': 'K5',     # all atoms from V4 
}

# No need to change below
# POKY libraries
import __main__
s = __main__.main_session
from sputil import name_to_spectrum
from sputil import parse_assignment, parse_group_name, split_group_atom

specnames = s.show_spectrumselectiondialog(
          'Select spectra to replace assignments', 1)
specname_list = specnames.split('\t')
if len(specname_list) == 0:
  raise SystemError

# Check spec names
spec_list = []
for specname in specname_list:
  sp = name_to_spectrum(specname, s)
  if sp == None:
    print(specname + ' selects no spectrum.')
    raise SystemError
  spec_list.append(sp)

# Check map_dict
map_list = list(map_dict.items())
for mp in map_list:
  k, v = mp
  try:
    g, n = parse_group_name(k)
  except:
    print('Map key error: ' + k)
    raise SystemError
  try:  
    g, n = parse_group_name(v)
  except:
    print('Map value error: ' + v)
    raise SystemError

# Loop spectra
for sp in spec_list:
  for peak in sp.peak_list():
    for i in range(len(peak.resonances())):
      r = peak.resonances()[i]
      if r == None:
        continue
      for mp in map_list:
        k, v = mp
        ga = split_group_atom(k)
        ga2 = split_group_atom(v)
        if ga == None:
          ga = [k, '']
        if ga2 == None:
          ga2 = [v, '']  
        if r.group.name == ga[0]:
          if ga[1] != r.atom.name and ga[1] != '':
            continue
          if ga2[1] == '':
            peak.assign(i, ga2[0], r.atom.name)
          else:
            peak.assign(i, ga2[0], ga2[1])
          break
s.command_characters('dr')
print('Done.')