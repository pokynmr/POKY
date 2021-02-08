#
# This is an example script to set molecule and condition name for many spectra.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project

# this is a keyword in spectrum names.
find_key = '_T1'

# molecule name to be used. it can be nothing if desired.
molecule_name = 'Nanog_APO'

# condition name to be used. it can be nothing if desired.
condition_name = 'T1'

# where it really changes
c = 0
for spec in proj.spectrum_list():
  if spec.name.find(find_key) != -1:
    spec.set_molecule_condition(molecule_name, condition_name)
    c+=1
    print('%s -> %s / %s' % (spec.name, molecule_name, condition_name))
print('%d experiment(s) found and set.\nDone.' % (c))
