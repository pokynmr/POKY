#
# This is an example script to run Poky-Potenci analysis.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

# Poky-Potenci extensively uses POTENCI from Mulder group.
# Reference:
#   POTENCI: prediction of temperature, neighbor and pH-corrected
#           chemical shifts for intrinsically disordered proteins
#   Jakob Toudahl Niesen and Frans A A Mulder
#   Journal of Biomolecular NMR. 2018 Mar; 70(3):141-165
#   doi.: 10.1007/s10858-018-0166-5
#
# POTENCI license: MIT
# Poky-Potenci license: BSD-2 for non-commercial academic work
#                       Poky Commercial License for commercial work


#
# Parameters
pH =              6.0
temperature =     292
ionic_strength =  0.5

sequence = '''
MQIFVKTLTGKTITLEVEPSDTIENVKAK
IQDKEGIPPDQQRLIFAGKQLEDGRTLSD
YNIQKESTLHLVLRLRGG
'''

# output
printout = True
writeout = False  # Set to True if you wish to write a file
#

if writeout:
  import __main__
  s = __main__.main_session
  outName = s.save_filedialog('Name a potenci output name.',
                              'Text file (*.txt);; Any (*)', '')
  if outName == '':
    raise SystemExit

# processing & run poky_potenci
sequence = ' '.join(sequence.split())
sequence = sequence.upper()

import poky_potenci
shiftdct, pkadict = poky_potenci.run(sequence, pH, temperature, ionic_strength)

if printout:
  poky_potenci.printOutput(shiftdct)
if writeout:
  from os.path import exists
  poky_potenci.writeOutput(shiftdct, outName)
  if exists(outName):
    import __main__
    s = __main__.main_session
    s.open_notepad(outName)

