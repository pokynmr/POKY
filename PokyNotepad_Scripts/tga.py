#
# This is an example script to transfer peaks to all spectra
# and save peak lists
#
# Saving entities:
# ppm, Hz, intensity, S/N
#
# You can customize as you wish.
#
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#         Make sure the source spectrum is selected,


# if this is not set, files will be saved in spectrum path.
import __main__
s = __main__.main_session
output_dir = s.open_directorydialog('Select a directory to save peak lists.', '')
#output_dir = '/path/to/save/your/files'

import os.path
from sputil import parse_assignment, split_group_atom

def save_peaks(spec):
  fname = os.path.basename(spec.data_path)
  a, b = os.path.splitext(fname)
  path = a + '.list'
  if output_dir != '':
    fullpath = os.path.join(output_dir, path)
  else:
    fullpath = os.path.join(os.path.dirname(spec.data_path) +
                        path)

  f = open(fullpath, 'w')
  for peak in spec.peak_list():
    line = '%-13s %7.3f %7.3f %10.3f %10.3f %11d %5d\n' % \
        (peak.assignment, peak.frequency[0], peak.frequency[1],
        peak.frequency[0] * spec.hz_per_ppm[0], 
        peak.frequency[1] * spec.hz_per_ppm[1],
        peak.data_height, int(peak.data_height / spec.noise))
    f.write(line)
  f.close()

def tga(session):
  ref_spec = session.selected_spectrum()
  save_peaks(ref_spec)
  ref_peaks = ref_spec.peak_list()
  for spec in session.project.spectrum_list():
    if ref_spec == spec:
      continue
    for peak in ref_peaks:
      peak2 = spec.place_peak(peak.frequency)
      if peak.is_assigned == 1:
        tasn = parse_assignment(peak.assignment)
        peak2.assign(0, tasn[0][0], tasn[0][1])
        peak2.assign(1, tasn[1][0], tasn[1][1])
        peak2.show_assignment_label()
    save_peaks(spec)

import __main__
s = __main__.main_session
tga(s)
print('Done.')
