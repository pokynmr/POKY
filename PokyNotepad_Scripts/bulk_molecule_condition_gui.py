#
# This is an example script to set molecule and condition name for many spectra.
# This is GUI version
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project

spec_names = s.show_spectrumselectiondialog('Select spectra', 1)
if spec_names == '':
  raise SystemExit
spec_list = spec_names.split('\t')

molecule_name = s.show_inputdialog('Molecule Name',
                                    'Type molecule name to use','')
condition_name = s.show_inputdialog('Condition Name',
                                    'Type condition name to use','')

if not s.show_message_yes_no('Change Molecule/Condition',
  'Do you want to change molecule/condition to %s/%s?' %
  (molecule_name, condition_name)):
  raise SystemExit

# where it really changes
c = 0
for spec in proj.spectrum_list():
  if spec.name in spec_list:
    spec.set_molecule_condition(molecule_name, condition_name)
    c+=1
    print('%s -> %s / %s' % (spec.name, molecule_name, condition_name))
print('%d experiment(s) found and set.\nDone.' % (c))
