#
# This is an example script to convert HX1/2 <-> HX2/3 system.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module


# This script does not require any modification to run.


import __main__
s = __main__.main_session
proj = s.project

# Buttons to be displayed. Return value is the index (0, 1, 2, ....)
selection = ['SYBYL', 'DIANA', 'XPLOR', 'BMRB', 'CANCEL']

# Selection dialog with buttons
cur=s.show_selectionexdialog('Hydrogen nomenclature',
  'Select current nomenclature', selection)
# Cancel pressed
if cur in [len(selection)-1, -1]:
  raise SystemExit;

# Selection dialog with buttons
tar=s.show_selectionexdialog('Hydrogen nomenclature',
    'Select target nomenclature', selection)
# Cancel pressed
if tar in [len(selection)-1, -1]:
  raise SystemExit;

# change from here
# this uses Eldon's table (very old)
# at least XPLOR and DIANA should be correct
# however it QX or CQX or %, * are tricky
# use with caution.

from nomenclature import translate

for sp in proj.spectrum_list():
  for p in sp.peak_list():
    for i in range(len(p.resonances())):
      r = p.resonances()[i]
      if r.atom.name == '': continue
      new_atm = translate(r.group.symbol, r.atom.name)
      if new_atm != '':
        print('%s -> %s%s' % (r.name, r.group.name, new_atm))
        p.assign(i, r.group.name, new_atom)
      else:
        if r.atom.name.find('CQ') == 0 and selection[tar] == 'XPLOR':
          new_atm = r.atom.name.replace('CQ', 'C') + '*'
          print('%s -> %s%s' % (r.name, r.group.name, new_atm))
          p.assign(i, r.group.name, new_atom)
        elif r.atom.name.find('Q') == 0 and selection[tar] == 'XPLOR':
          new_atm = r.atom.name.replace('Q', 'H') + '*'
          print('%s -> %s%s' % (r.name, r.group.name, new_atm))
          p.assign(i, r.group.name, new_atom)
        elif r.atom.name[-1] in ['*', '%']:
          new_atom = r.atom.name.replace('*', '').replace('%', '')
          print('%s -> %s%s' % (r.name, r.group.name, new_atm))
          p.assign(i, r.group.name, new_atom)

s.command_characters('dr')
