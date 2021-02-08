# transfer peaks to all spectra and save peak lists
# Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# 01/20/2021
#
# Usage:
#    In NMRFAM-SPARKY,
#         Click "Load_Module" button and select this file
#         Make sure the source spectrum is selected,
#         and type:
#         tga.tga(s)

import os.path
from sputil import parse_assignment, split_group_atom

def save_peaks(spec):
  a, b = os.path.splitext(spec.save_path)
  path = a + '.list'
  f = open(path, 'w')
  for peak in spec.peak_list():
    line = '%-13s %7.3f %7.3f %10.3f %10.3f %11d %5d\n' % \
        (peak.assignment, peak.frequency[0], peak.frequency[1],
        peak.frequency[0] * spec.hz_per_ppm[0], 
        peak.frequency[1] * spec.hz_per_ppm[1],
        peak.data_height, int(peak.data_height / spec.noise))
    f.write(line)
  f.close()

def tpa(session):
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