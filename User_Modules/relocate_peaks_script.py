#
# This is an example script to relocate peaks to the proper position.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

spec_names = s.show_spectrumselectiondialog('Select spectra', 1)
if spec_names == '':
  raise SystemError
specname_list = spec_names.split('\t')

tols = s.show_inputdialog('Tolerances',
          'Tolerances (C, N, H order in ppms) to navigate:', '3, 3, 0.3')
tol_list = tols.split(',')

try:
  Ctol = float(tol_list[0].strip())
  Ntol = float(tol_list[1].strip())
  Htol = float(tol_list[2].strip())
except:
  print('Tolerance format is incorrect.')
  raise SystemError

from sputil import relocate_peaks, name_to_spectrum
spec_list = list(map(lambda sp: name_to_spectrum(sp, s), specname_list))
for spec in spec_list:
  relocate_peaks(spec, peaks=None, Ctol=Ctol, Ntol=Ntol, Htol=Htol)