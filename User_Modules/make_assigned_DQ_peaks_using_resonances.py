#
# This is an example script to make assigned DQ peaks using resonances.
# DQ options adapted from Goldbourt's JACS Au paper work.
# https://pubs.acs.org/doi/10.1021/jacsau.4c00549
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

from sputil import name_to_condition, name_to_spectrum
import __main__
s = __main__.main_session

DQ_options = {}
DQ_options['G'] = [["C","CA"]]
DQ_options['A'] = DQ_options['C'] = DQ_options['S'] = [["C","CA"], ["CA","CB"]]
DQ_options['F'] = DQ_options['Y'] = [["C","CA"], ["CA","CB"], ["CB","CG"], 
    ["CG","CD1"], ["CG","CD2"], ["CD1","CE1"], ["CD2","CE2"], ["CE1","CZ"],
    ["CE2","CZ"]]
DQ_options['W'] = [["C","CA"], ["CA","CB"], ["CB","CG"], ["CG","CD1"],
    ["CG","CD2"], ["CD2","CE2"], ["CD2","CE3"], ["CE2","CZ2"], 
    ["CE3","CZ3"], ["CZ2","CH2"], ["CZ3","CH2"]]
DQ_options['H'] = [["C","CA"], ["CA","CB"], ["CB","CG"], ["CG","CD2"]]
DQ_options['V'] = [["C","CA"], ["CA","CB"], ["CB","CG1"], ["CB","CG2"]]
DQ_options['L'] = [["C","CA"], ["CA","CB"], ["CB","CG"], ["CG","CD1"],
                    ["CG","CD2"]]
DQ_options['I'] = [["C","CA"], ["CA","CB"], ["CB","CG1"], ["CB","CG2"],
                    ["CG1","CD1"]]
DQ_options['K'] = [["C","CA"],["CA","CB"],["CB","CG"],["CG","CD"],["CD","CE"]]                
DQ_options['R'] = DQ_options['E'] = DQ_options['Q'] = DQ_options['P'] = \
    [["C","CA"], ["CA","CB"], ["CB","CG"], ["CG","CD"]]
DQ_options['D'] = DQ_options['N'] = DQ_options['M'] = \
    [["C","CA"], ["CA","CB"], ["CB","CG"]]
DQ_options['T'] = [["C","CA"], ["CA","CB"],["CB","CG2"]]


specname = s.show_spectrumselectiondialog('Select a target spectrum', 0)
spec = name_to_spectrum(specname, s)
if spec == None:
    raise SystemError

if spec.dimension != 2 or '13C' not in spec.nuclei:
    s.show_message('Error', 'Selected spectrum not supported.')
    raise SystemError

dimDQ = s.show_selectionexdialog('Dimension', 'DQ dimension: ', 
                            ['w1', 'w2', 'Cancel'])
if dimDQ == 2:
    raise SystemError

dim = 1 - dimDQ

clist = s.project.condition_list()
if len(clist) > 1:
    cname = s.show_conditionselectiondialog('Select a condition.', 0)
    c = name_to_condition(cname, s)
else:
    c = clist[0]

from myseq import BuildSequence
seq_list = BuildSequence(s, 3)

freqs = [0, 0]
atms = ['', '']
for nseq, a in seq_list:
    try:
        option_list = DQ_options[a]
    except:
        continue
    for atm1, atm2 in option_list:
        r1 = c.find_resonance(f'{a}{nseq}', atm1)
        r2 = c.find_resonance(f'{a}{nseq}', atm2)

        if None in [r1, r2]:
            continue
        
        atmDQ = f'{atm1}&{atm2}'
        csDQ = r1.frequency + r2.frequency
        
        freqs[dimDQ] = csDQ
        atms[dimDQ] = atmDQ
        
        for r in [r1, r2]:
            freqs[dim] = r.frequency
            atms[dim] = r.atom.name
            peak = spec.place_peak(freqs)
            peak.assign(0, f'{a}{nseq}', atms[0])
            peak.assign(1, f'{a}{nseq}', atms[1])
            peak.show_assignment_label()

s.show_message('Done', 'Finished.')