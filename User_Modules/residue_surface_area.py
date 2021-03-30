#
# This is an example script to calculates the surface area in square Angstroms of the selected residue(s).
# by Mehdi Rahimi, Ph.D. (mehdi.rahimi@ucdenver.edu)
# based on "ponderosa_rmsd" script by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
import os.path
import re 
from pymol import cmd

s = __main__.main_session

PDB_selection = s.show_selectiondialog('Do you want to use the PONDEROSA result?', 'Or, do you have a PDB file you want to use?', ('Use a PDB file', 'Use PONDEROSA result'))

if PDB_selection == 1:
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
    pdbpath = s.open_filesdialog('Select the PDB file', 'PDB File (*.pdb);; Any (*);;', '.')
    if pdbpath == '':
        raise SystemExit
  
selection = s.show_inputdialog('Select the residue numbers','Select a residue number or a "range" using a dash (ex: 5-8)', '0')

correct_selection = False
if re.match(r"^\d+$|^\d+-\d+$", selection):
    correct_selection = True
	
if selection == '0' or not correct_selection:
    s.show_message('Error', 'The residue numbers selection format is invalid.')
    print("The residue numbers selection format is invalid.")
    raise SystemExit

cmd.load(pdbpath, 'for_area')
cmd.set('dot_solvent', 1)
area = cmd.get_area('resi ' + selection)

if area == 0:
    s.show_message('Error', 'The residue number is invalid.')
    print("The residue number is invalid.")
    raise SystemExit
    
print(area)
s.show_message('Finished', 'The surface area for residue ' + selection + ' is ' + '{:.3f}'.format(area) + ' Angstroms^2')

