#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
This program shows the protein contact map for a PDB file

By Mehdi Rahimi, Ph.D. (mehdi.rahimi@ucdenver.edu)

To run this script in Poky:
    In Poky Notepad,
    File -> Run Python Module
"""

from math import sqrt
import matplotlib.pyplot as plt
import __main__
s = __main__.main_session


def distance3D(atom1, atom2):
    """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                   (26.913, 26.639, -3.51))
        returns the distance
    """
    return sqrt((atom1[0] - atom2[0]) ** 2 +
                (atom1[1] - atom2[1]) ** 2 +
                (atom1[2] - atom2[2]) ** 2)



def readPDB(pdb, modelnumber):
    """
    Parameters
    ----------
    pdb : address to a .pdb file
    modelnumber : the model number in the pdb (starting from 1)

    Returns
    -------
    pdb_list: a list of each [atom's coordinates, aa, atom]

    """
    pdbLines, modelList = [], []
    tempLines = open(pdb, 'r').readlines()

    # clean
    for line in tempLines:
        if line[0:4] in ['MODE', 'ATOM', 'ENDM']:
            pdbLines.append(line)

    # fill modelList
    for line in pdbLines:
        if line[0:5] == 'MODEL':
            modelList.append([])
        if line[0:4] != 'ATOM':
            continue
        if line[12:16].strip()[0] not in ['C', 'N', 'H', 'O']:
            continue
        aaa = line[17:20].strip()
        atm = line[12:16].strip()
        nSeq = int(line[23:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        modelList[-1].append( [nSeq, x, y, z, aaa, atm] )

    return modelList[modelnumber - 1]



def createContactMap(pdb, contact, dist_cutoff):
    """
    Generates protein contact map matrix

    Parameters
    ----------
    pdb : a list of each [atom's coordinates, aa, atom]
    contact : contact type ('CA' or 'CB')
    dist_cutoff : the contact distance cutoff in Angstrom (6-12)

    Returns
    -------
    cm : the 2D matrix of the contact map
    res_range : x,y axis values for the map

    """
    if contact[-1] == '*':
        contacts = [row for row in pdb if row[5].startswith(contact[:-1])]
    else:
        contacts = [row for row in pdb if row[5]==contact]

    res_range = []
    for n in contacts:
        res_range.append(n[0])
    res_range = list(set(res_range))

    cm = [([0] * len(res_range)) for i in range(len(res_range))]

    for row_i in contacts:
        i_pos = (row_i[1], row_i[2], row_i[3])

        for row_j in contacts:
            j_pos = (row_j[1], row_j[2], row_j[3])

            distance = distance3D(i_pos, j_pos)

            if distance < dist_cutoff:
                cm[res_range.index(row_i[0])][res_range.index(row_j[0])] = 1
                # the pcolor matrix cannot have empty rows/columns
                # and should match the res_range

    return cm, res_range


pdb_file = s.open_filedialog('Select a PDB file.',
                             'PDB file (*.pdb)', '')

if pdb_file == '' or not pdb_file.endswith('.pdb'):
    print("You didn't select a PDB file")
    raise SystemExit


model_number = s.show_inputdialog('Model number',
                                  'Select the model number in PDB',
                                  '1')

msg = 'Select the contact type for the contact map'
atom_list = ('CA', 'CB', 'H', 'HA*', 'HB*', 'H*')
contact_selection = s.show_selectionexdialog(msg, '', atom_list)
try:
    contact_type = atom_list[contact_selection]
except:
    print("You didn't select the contact type")
    raise SystemExit


msg = 'Select the contact distance cutoff in Angstrom (6-12)'
dist = s.show_inputdialog('Distance cutoff', msg, '8')


pdb_list = readPDB(pdb_file, int(model_number))
contact_map, residue_range = createContactMap(pdb_list,
                                              contact_type,
                                              float(dist))


cmap = plt.get_cmap("binary")
residue_range.append(residue_range[-1] + 1)  # pcolor needs one more value
plt.pcolor(residue_range, residue_range, contact_map, cmap=cmap)
plt.show(block=False)

