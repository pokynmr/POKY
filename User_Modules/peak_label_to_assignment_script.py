#
# This is an example script to make assignments from peak labels.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#

from sputil import name_to_spectrum, parse_assignment

# POKY libraries
import __main__
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

ans = s.show_selectionexdialog('Peak label to assignment',
          'Change labels to assignments for',
          ('All labeled peaks', 'Selected peaks', 'Cancel'))
if ans in [-1, 2]:
  raise SystemError

counter = 0
for p in sp.peak_list():
  if p.label == None:
    continue
  if ans == 1 and p.selected == 0:
    continue
  lbl = p.label.text
  asgn = parse_assignment(lbl)
  if asgn == None:
    continue
  try:
    for i in range(sp.dimension):
      p.assign(i, asgn[i][0], asgn[i][1])
  except:
    continue
  p.show_assignment_label()
  counter += 1

s.show_message('Done', 
  f'{counter} peak label(s) have been changed to assignment(s).')