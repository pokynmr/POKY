#
# This is an example script to make peak labels from assignments.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#

from sputil import name_to_spectrum

# POKY libraries
import __main__
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

ans = s.show_selectionexdialog('Assignment to peak label',
          'Change assignments to labels for',
          ('All assigned peaks', 'Selected peaks', 'Cancel'))
if ans in [-1, 2]:
  raise SystemError

counter = 0
for p in sp.peak_list():
  if p.is_assigned != 1:
    continue
  if ans == 1 and p.selected == 0:
    continue
  lbl = p.assignment
  for i in range(sp.dimension):
    p.assign(i, '', '')
  p.show_label(lbl)
  counter += 1

s.show_message('Done', 
  f'{counter} assignment(s) have been changed to peak label(s).')
