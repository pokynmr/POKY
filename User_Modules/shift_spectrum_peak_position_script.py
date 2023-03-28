#
# This is an example script for shifting peak positions on a spectrum.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

# POKY libraries
import __main__
s = __main__.session

from sputil import name_to_spectrum

specname = s.show_spectrumselectiondialog('Select a spectrum', 0)
sp = name_to_spectrum(specname, s)

if sp == None:
  raise SystemError

ndim = sp.dimension
sht = ', '.join(['0',]*ndim)

ans = s.show_inputdialog('Peak Shift',
  'Please specify shift amounts for peaks (PPMs for all dimensions)',
  sht)

ans_list = ans.split(',')
if len(ans_list) != ndim:
  raise SystemError

shifts = list(map(lambda x: float(x.strip()), ans_list))

for p in sp.peak_list():
  new_freqs = list(map(lambda i: p.frequency[i] + shifts[i], range(ndim)))
  p.frequency = tuple(new_freqs)

s.show_message('Finished', 
               'Peak shifting is completed. ')
