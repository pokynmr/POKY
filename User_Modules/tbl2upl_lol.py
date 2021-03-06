#!/usr/bin/env python3
# pylint: disable=invalid-name
"""
This program converts XPLOR TBL (simple or ambi) to DIANA UPL and LOL
The TBL file must be single-line single-restraint.

By Mehdi Rahimi, Ph.D. (mehdi.rahimi@ucdenver.edu)
Based on a script by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)

To run this script in Poky:
    In Poky Notepad,
    File -> Run Python Module
To run in the Terminal:
    python tbl2upl_lol.py [TBL FILE] [PDB FILE] [Sequence FILE]
"""

import sys
import os.path
import __main__
import re
import math

#from atomnames import s12s32
#from myseq import A_dict, AAA_dict

s12s32 = (('G', 'HA1', 'HA3'), ('C', 'HB1', 'HB3'), ('D', 'HB1', 'HB3'),
          ('E', 'HB1', 'HB3'), ('E', 'HG1', 'HG3'), ('F', 'HB1', 'HB3'),
          ('H', 'HB1', 'HB3'), ('I', 'HG11', 'HG13'), ('K', 'HB1', 'HB3'),
          ('K', 'HD1', 'HD3'), ('K', 'HG1', 'HG3'), ('K', 'HE1', 'HE3'),
          ('L', 'HB1', 'HB3'), ('M', 'HB1', 'HB3'), ('M', 'HG1', 'HG3'),
          ('N', 'HB1', 'HB3'), ('P', 'HB1', 'HB3'), ('P', 'HD1', 'HD3'),
          ('P', 'HG1', 'HG3'), ('Q', 'HB1', 'HB3'), ('Q', 'HG1', 'HG3'),
          ('R', 'HB1', 'HB3'), ('R', 'HD1', 'HD3'), ('R', 'HG1', 'HG3'),
          ('S', 'HB1', 'HB3'), ('W', 'HB1', 'HB3'), ('Y', 'HB1', 'HB3'))
A_dict={'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
        'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
        'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
        'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
AAA_dict={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def one_or_three(seq_file):
    tmp = open(seq_file, 'r').read()
    tmp2 = tmp.upper()
    regex = re.compile('[^A-Z]')
    seq = regex.sub('', tmp2)
    if len(seq) % 3 != 0:
        return 3, seq
    for i in range(int(len(seq) / 3)):
        subseq = seq[i*3:i*3+3]
        try:
            AAA_dict[subseq]
        except IndexError:
            return 1, seq
    return 3, seq


def extract_seq_file(seq_file):
    seq_mode, seq = one_or_three(seq_file)
    a_list = list(map(lambda x: '', range(1000)))
    if seq_mode == 1:
        for i, row in enumerate(seq):
            a_list[i] = row
        return a_list
    regex = re.compile('[^A-Z+0-9+ ]')
    lines = open(seq_file, 'r').readlines()
    iCur = 1
    for line in lines:
        if line.strip() == '':
            continue
        if line.strip()[0] == '#':
            continue
        tmp = regex.sub('', line.upper()).split()
        if len(tmp) == 2:
            try:
                iCur = int(tmp[0])
                a_list[iCur] = AAA_dict[tmp[1]]
                iCur = iCur + 1
            except:
                iCur = int(tmp[1])
                a_list[iCur] = AAA_dict[tmp[0]]
                iCur = iCur + 1
                continue
        if len(tmp) == 1:
            a_list[iCur] = AAA_dict[tmp[0]]
            iCur = iCur + 1
            continue

    return a_list


def translate12to23(aaa, atm):
    if atm == 'HN':
        return 'H'
    if atm == 'CO':
        return 'C'
    if len(aaa) == 3:
        a = AAA_dict[aaa]
    else:
        a = aaa
    for grp in s12s32:
        if (grp[0] == a) and (grp[1] == atm):
            return grp[2]
    return atm


def tbl2upl_lol(tbl_file, seq_file):
    """
    Reads the TBL file and create the UPL and LOL as text

    """
    lines = open(tbl_file, 'r').readlines()
    upl_text = ''
    lol_text = ''
    a_list = extract_seq_file(seq_file)

    tbl_list = []
    for line in lines:
        line = line.replace('(','').replace(')','').replace('"','')
        if line == '' or line[0] in '!#/':
            continue
        line = line.replace('#','').replace('*','')
        tbl_list.append(line.split())

    last_assign = 0
    for i, row in enumerate(tbl_list):
        n1 = int(row[5])
        n2 = int(row[13])
        a1 = translate12to23(a_list[n1], row[8])
        a2 = translate12to23(a_list[n2], row[16])
        s1 = A_dict[a_list[n1]]
        s2 = A_dict[a_list[n2]]

        if row[0] == 'assign':
            last_assign = i
        else:   # for a 'or' line in ambi TBL
            row = tbl_list[last_assign]

        upl_dist = float(row[17]) + float(row[19])
        lol_dist = float(row[17]) - float(row[18])

        # check the PDB to see if the constraint makes sense
        for peak in pdb_list:
            if peak[0] == n1 and peak[2] == n2 and \
               peak[1] == a1 and peak[3] == a2:
                if upl_dist > peak[4] and lol_dist < peak[4]:
                # to make sure the distance is within the bounds

                    upl_text += '%4d %4s %4s %4d %4s %4s %8.2f\n' % \
                                (n1, s1, a1, n2, s2, a2, upl_dist)

                    lol_text += '%4d %4s %4s %4d %4s %4s %8.2f\n' % \
                                (n1, s1, a1, n2, s2, a2, lol_dist)

                break

    return upl_text.strip(), lol_text.strip()


def read_pdb(pdb_file):
    """ open a PDB file and get the atom positions """

    modelList = []
    tempLines = open(pdb_file, 'r').readlines()
    pdbLines = []

    # clean
    for line in tempLines:
        if len(line) < 4:
            continue
        try:
            if not line[0:4] in ['MODE', 'ATOM', 'ENDM']:
                continue
        except IndexError:
            continue
        pdbLines.append(line.strip())

    # fill modelList
    for line in pdbLines:
        if line[0:5] == 'MODEL':
            modelList.append([])
        if line[0:4] != 'ATOM':
            continue
        if not line[12:16].strip()[0] in ['C', 'N', 'H', 'O']:
            continue
        aaa = line[17:20].strip()
        atm = line[12:16].strip()
        nSeq = int(line[23:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        try:
            modelList[-1].append( [nSeq, x, y, z, aaa, atm] )
        except IndexError:
            modelList.append([])
            modelList[-1].append( [nSeq, x, y, z, aaa, atm] )
    return modelList[0]


def distance3D(atoms):
    """ takes a tuple of two coordinates: ([26.266, 25.413, 2.842],
                                            26.913, 26.639, 3.531])
        returns the distance
    """
    return math.sqrt((atoms[0][0] - atoms[1][0]) ** 2 +
                     (atoms[0][1] - atoms[1][1]) ** 2 +
                     (atoms[0][2] - atoms[1][2]) ** 2)


def get_distances(pdb):
    """
    Takes the pdb list (proccessed from PDB file)

    Returns
    -------
    distances_list : List
        [Carbon_Atom1,
         Carbon_Atom1_Residue#,
         Carbon_Atom2,
         Carbon_Atom2_Residue#,
         3D Distance between them,
         Atom1_postion,
         Atom2_postion,
         AminoAcid1,
         AminoAcid2]

    Sorted by the distance
    """
    c_list = ['C', 'CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2',
              'CE', 'CE1', 'CE2', 'CE3', 'CH', 'CH1', 'CH2', 'CZ']

    distances_list = []

    for atom1 in pdb:
        if atom1[5] in c_list:
            for atom2 in pdb:
                if atom2[5] in c_list:
                    distance = distance3D([[atom1[1],atom1[2],atom1[3]],
                                           [atom2[1],atom2[2],atom2[3]]])
                    distances_list.append([atom1[0],
                                           atom1[5],
                                           atom2[0],
                                           atom2[5],
                                           distance,
                                           [atom1[1],atom1[2],atom1[3]],
                                           [atom2[1],atom2[2],atom2[3]],
                                           AAA_dict[atom1[4]],
                                           AAA_dict[atom2[4]]
                                           ])

    return sorted(distances_list, key=lambda x: x[4])


def list_cutoff(distances_list):
    """
    Takes the distance list and remove row above the cutoff
    """
    DISTANCE_CUTOFF = 7.8

    final_row = len(distances_list)
    for i, atom in enumerate(distances_list):
        if atom[4] > DISTANCE_CUTOFF:
            final_row = i
            break

    return distances_list[:final_row]


def sort_and_write(text, output_file):
    """
    Gets the text for UPL or LOL file, sorts it and writes to disk

    """
    # sorting
    text_list = text.split('\n')
    for i, row in enumerate(text_list):
        text_list[i] = row.split()
    text_list.sort(key=lambda x: (int(x[0]), int(x[3]), x[1], x[4], x[2], x[5], x[6]))

    output = ''
    for r in text_list:
        line = '%4d %4s %4s %4d %4s %4s %8.2f' % \
            (int(r[0]), r[1], r[2], int(r[3]), r[4], r[5], float(r[6]))
        output += line
        output += '\n'


    f = open(output_file, 'w')
    f.write(output)
    f.close()


#################################################################

try:  
    # if  we're in Poky
    s = __main__.session
    in_poky = True
except AttributeError: 
    # if we're running this code in the Terminal
    in_poky = False


if in_poky:
    dist_file = s.open_filedialog('Select a distance restraint file.',
                                  'XPLOR TBL (*.tbl)', '')
    if dist_file == '' or not dist_file.endswith('.tbl'):
        print("You didn't select a TBL file")
        raise SystemExit

    p = os.path.dirname(dist_file)
    sequence_file = s.open_filedialog('Select a sequence file.',
                                      'Sequence file (*.seq *.fasta *.txt);; Any (*)', p)
    pdb_file = s.open_filedialog('Select the PDB file.',
                                 'PDB File (*.pdb)', '')
    if pdb_file == '':
        print("You didn't select a PDB file")
        raise SystemExit

    upl_output_file = s.save_filedialog('Save the UPL file as...',
                                        'DIANA UPL (*.upl);; Any (*)', p)
    lol_output_file = s.save_filedialog('Save the LOL file as...',
                                        'DIANA LOL (*.lol);; Any (*)', p)
    
    if not upl_output_file.endswith('.upl'):
           upl_output_file += '.upl'
           
    if not lol_output_file.endswith('.lol'):
           lol_output_file += '.lol'
                 
    s.show_message('The program is running', 'Click OK and wait for the program to finish')                                    
    pdb_list = read_pdb(pdb_file)
         
else:   
    if len(sys.argv) != 4:
        print("python tbl2upl_lol.py [TBL FILE] [PDB FILE] [Sequence FILE]")
        print("e.g. python tbl2upl_lol.py ambi.tbl 5mwv.pdb 34088.seq")
        sys.exit()

    dist_file = sys.argv[1]
    pdb_list = read_pdb(sys.argv[2])
    sequence_file = sys.argv[3]

    upl_output_file = 'dist.upl'
    lol_output_file = 'dist.lol'


# dist_file = 'ambi.tbl'
# pdb_list = read_pdb('5mwv.pdb')
# sequence_file = '34088.seq'

dist_list = get_distances(pdb_list)
pdb_list = list_cutoff(dist_list)

upl, lol = tbl2upl_lol(dist_file, sequence_file)

# UPL file
sort_and_write(upl, upl_output_file)

# LOL file
sort_and_write(lol, lol_output_file)

print("The UPL and LOL files are saved.")

