#
# This is an example script to calculates the surface area
#   in square Angstroms of the selected residue(s).
# by Mehdi Rahimi, Ph.D. (mehdi.rahimi@ucdenver.edu)
# based on "ponderosa_rmsd" script by
#  Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import os.path
import re
import numpy as np
import __main__
from pymol import cmd

s = __main__.main_session

PDB_PO = s.show_selectiondialog('Do you want to use the PONDEROSA result?',
                                'Or, do you have a PDB file you want to use?',
                                ('Use a PDB file', 'Use PONDEROSA result'))

if PDB_PO == 1:
    pdir = s.open_directorydialog('Select the PONDEROSA result directory', '')
    if pdir == '':
        raise SystemExit

    pbdir = os.path.join(pdir, 'BestEvaluated')
    pdbpath = os.path.join(pbdir, 'final_water_refined.pdb')
    if not os.path.exists(pdbpath):
        pdbpath = os.path.join(pbdir, 'final.pdb')

    if not os.path.exists(pdbpath):
        s.show_message('Error', 'This is not a PONDEROSA result directory')
        raise SystemExit

else:
    pdbpath = s.open_filesdialog('Select the PDB file',
                                 'PDB File (*.pdb);; Any (*);;', '.')
    if pdbpath == '':
        raise SystemExit

in_msg = 'Select a residue number or a "range" using a dash (ex: 5-8)'
selection = s.show_inputdialog('Select the residue numbers', in_msg, '0')

correct_selection = False
if re.match(r"^\d+$|^\d+-\d+$", selection):
    correct_selection = True

if selection == '0' or not correct_selection:
    s.show_message('Error', 'The residue numbers selection format is invalid.')
    print("The residue numbers selection format is invalid.")
    raise SystemExit

in_msg = 'Select the sampling density from 1-4. \n\
Higher density means higher accuracy but slower performance.'
density = int(s.show_inputdialog('Select the Sampling density', in_msg, '2'))

if density < 1 or density > 4:
    s.show_message('Error', 'The sampling density is invalid.')
    raise SystemExit

in_msg = 'Select the model number. \n\
Use 1 if there\'s only one model. Use 0 for all models.'
model = int(s.show_inputdialog('Select the model', in_msg, '1'))

cmd.delete('all')
cmd.load(pdbpath, 'for_area')
cmd.split_states('for_area')
cmd.set('dot_solvent', 1)
cmd.set('dot_density', density)

if model != 0:  # for a specific model
    try:
        area = cmd.get_area('resi ' + selection +
                            ' and model for_area_' + '{:04d}'.format(model))
    except:
        s.show_message('Error', 'The model number was wrong.')
        raise SystemExit

    print(area)

    message = 'The surface area for residue ' + selection + ' is ' + \
               '{:.3f}'.format(area) + ' Angstroms^2'
    if area == 0:
        message += "\nYou'd also get a zero if the residue number is invalid."

    s.show_message('Finished', message)


else: # for all

    areas = []
    for i in range(1, 1000):
        try:
            area = cmd.get_area('resi ' + selection +
                                ' and model for_area_' + '{:04d}'.format(i))
        except:
            break
        areas.append(area)
        
    areas_mean = np.mean(areas)
    areas_std = np.std(areas)
    
    print('Area Mean: ' + str(areas_mean))
    print('Area StDev: ' + str(areas_std))

    message = 'The surface area for residue ' + selection + \
              ' for all models, has a mean of ' + \
              '{:.3f}'.format(areas_mean) + \
              ' with the standard deviation of ' + \
              '{:.3f}'.format(areas_std) + ' Angstroms^2'
    s.show_message('Finished', message)

