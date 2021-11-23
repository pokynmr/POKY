#
# This is an example script to convert ucsf to nmrpipe file.
# Adapted from NMRGlue exmple.
# https://nmrglue.readthedocs.io/en/latest/jbnmr_examples/s1_sparky_to_nmrpipe.html

#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# This script can take long if your data size is large. Be patient!
# It is a good idea to save your project before you run this script.

print('\n\n\n------------------------------------------------------')
print('POKY UCSF2PIPE conversion')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

from sputil import name_to_spectrum
import __main__
s = __main__.main_session

# User parameter
# Choose either the ucsf path or spectrum name
#ucsfpath = s.open_filedialog('Select a spectrum', 'Any (*);; UCSF (*.ucsf)', '')
specname = s.show_spectrumselectiondialog('Select a spectrum', 0)

# In case, if you wish to convert multiple files
# pipepath will be always considered "None"
#specs = s.show_spectrumselectiondialog('Select spectra', 1)
#specnames = specs.split('\t')

# Choose either the pipe path to write or None for writing in the same
# directory with .pipe extension
#pipepath = '/path/to/your/spectrum.pipe'
#pipepath = None
pipepath = s.save_filedialog('New pipe path', 'Any (*);; PIPE (*.pipe)', '')

# Processing start
try:
  getattr("specnames")
  pipepath = None
except:
  specnames = [specname,]

import os
import nmrglue as ng
for specname in specnames:
  try:
    spec = name_to_spectrum(specname, s)
    ucsfpath = spec.data_path
  except:
    pass

  if pipepath == None:
    from pathlib import Path
    pre, ext = os.path.splitext(ucsfpath)
    pipepath = pre + ".pipe"

  # read in the Sparky file
  print('Convert %s to %s' % (ucsfpath, pipepath))
  sdic, sdata = ng.sparky.read(ucsfpath)

  # convert to NMRPipe format
  C = ng.convert.converter()
  C.from_sparky(sdic, sdata)
  pdic, pdata = C.to_pipe()

  # write results to NMRPipe file
  ng.pipe.write(pipepath, pdic, pdata, overwrite=True)
