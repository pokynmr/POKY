#
# This is an example script to plot contours from selected view.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project

print('\n\n\n------------------------------------------------------')
print('POKY Contour Drawer')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

from plot_strips import draw_spectrum

plotting_methods = ('2D contour', '3D wireframe', '3D surface', '3D scatter', '3D contour', 'Cancel')
dr = s.show_selectionexdialog('Plotting', 'Plotting method', plotting_methods)
if dr in [-1, 5]:
  raise SystemExit

if dr == 0:
  draw_spectrum(s)
else:
  draw_spectrum(s, None, True, plotting_methods[dr][3:])

# Alternatively a user can plot by view name
# Uncomment below to try
#
#from sputil import name_to_view
#v = name_to_view('HNCACB', s) # this should be your view name
#if v:
#  draw_spectrum(s, v)
