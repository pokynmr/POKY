#
# This is an example script to show and hide views.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

print('\n\n\n------------------------------------------------------')
print('POKY View Visibility')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

import __main__
s = __main__.main_session
proj = s.project

hide_list = [ # list view names here
              # 'HNCACB', 'CBCACONH'
            ]
show_list = [ # list view names here
              # 'HSQC', 'HNCO'
            ]

for v in proj.view_list():
  if v.name in show_list:
    v.is_shown = 1
  elif v.name in hide_list:
    v.is_shown = 0

## If you want to use spectrum names, you can do like this:
# for v in proj.view_list():
#   if v.spectrum.name in show_list:
#     v.is_shown = 1
#   elif v.spectrum.name in hide_list:
#     v.is_shown = 0

## If you want to use condition names, you can do like this:
# for v in proj.view_list():
#   if v.spectrum.condition == None:
#     continue
#   if v.spectrum.condition.name in show_list:
#     v.is_shown = 1
#   elif v.spectrum.condition.name in hide_list:
#     v.is_shown = 0

## If you want to use molecule names, you can do like this:
# for v in proj.view_list():
#   if v.spectrum.molecule == None:
#     continue
#   if v.spectrum.molecule.name in show_list:
#     v.is_shown = 1
#   elif v.spectrum.molecule.name in hide_list:
#     v.is_shown = 0
