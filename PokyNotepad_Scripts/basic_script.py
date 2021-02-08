#
# This is an example script to show basic information from Poky.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project

print('\n\n\n------------------------------------------------------')
print('POKY Project Information')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# Show project information
print('* Project Information')
print('Project Path: ' + proj.save_path)
print('Spectrum Count: ' + str(len(proj.spectrum_list())))
print('View Count: ' + str(len(proj.view_list())))
print('Molecule Count: ' + str(len(proj.molecule_list())))
print('Condition Count: ' + str(len(proj.condition_list())))

# Show spectrum information
print('\n* Spectrum Information')
for spec in proj.spectrum_list():
  print('\n[Spectrum: %s]' % (spec.name))
  print('Data Path: ' + str(spec.data_path))
  print('Save Path: ' + str(spec.save_path))
  print('Description: ' + spec.desc)
  print('Dimension: ' + str(spec.dimension))
  print('Data Size: ' + str(spec.data_size))
  print('Nuclei: ' + str(spec.nuclei))
  print('Sweep Width: ' + str(spec.sweep_width))
  print('Region: ' + str(spec.region))
  print('Offset: ' + str(spec.scale_offset))
  print('Noise: ' + str(spec.noise))
  print('Peak Count: ' + str(len(spec.peak_list())))
  
  # Show peak information
  for peak in spec.peak_list():
    print((peak.assignment, peak.frequency, peak.data_height))

# View information
print('\n* View Information')
for view in proj.view_list():
  print('\n[View: %s]' % (view.name))
  print('Spec Name: ' + str(view.spectrum.name))
  print('Center: ' + str(view.center))
  print('Region: ' + str(view.region))
  print('Visible Depth: ' + str(view.visible_depth))
  print('Axis Order: ' + str(view.axis_order))
  pl = view.positive_levels
  print('Positive Levels: ' + str((pl.lowest, pl.factor, pl.levels, pl.color)))
  nl = view.negative_levels
  print('Negative Levels: ' + str((nl.lowest, nl.factor, nl.levels, nl.color)))

