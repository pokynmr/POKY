#
# This is an example script for adding/deleting random peaks on a spectrum.
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

opt_list = ('-100%', '-90%', '-80%', '-70%', '-60%',
            '-50%', '-40%', '-30%', '-20%', '-10%', 
            '+10%', '+20%', '+30%', '+40%', '+50%', 
            '+60%', '+70%', '+80%', '+90%', '+100%', 'Cancel')
ans = s.show_selectionexdialog('Add or delete random peaks', 'Select an option', 
                         opt_list)

if ans in [-1, len(opt_list)-1]:
  raise SystemError
ans = opt_list[ans]

specnames = s.show_spectrumselectiondialog('Select a spectrum', 1)
specnames_list = specnames.split('\t')
if len(specnames_list) == 0:
  raise SystemError

specs = list(map(lambda x: name_to_spectrum(x, s), specnames_list))

for p in s.selected_peaks():
  p.selected = 0

msg = ''
import random
for spec in specs:
  peaks = spec.peak_list()
  npeaks = len(peaks)
  
  if npeaks == 0:
    msg += f'failed: {spec.name} has no peaks.\n'
    continue
  ndiff = int(npeaks * float(ans[1:-1]) / 100. + 0.5)
  if ans[0] == '+':
    for i in range(ndiff):
      freq = []
      for dim in range(spec.dimension):
        f = random.uniform(spec.region[0][dim], spec.region[1][dim])
        freq.append(f)
      spec.place_peak(freq)
    msg += f'added {ndiff} peaks on {spec.name}: {npeaks} -> {npeaks+ndiff}\n'
  elif ans[0] == '-':
    for i in range(ndiff):
      while True:
        idx = random.randint(0, npeaks-1)
        if peaks[idx].selected == 0:
          peaks[idx].selected = 1
          break
    s.command_characters("")
    msg += f'deleted {ndiff} peaks on {spec.name}: {npeaks} -> {npeaks-ndiff}\n'

s.show_message('Finished.', msg)
