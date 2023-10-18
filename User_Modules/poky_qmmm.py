#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 13:25:26 2022

@author: abbychiu
"""

# This is a POKY Notepad script to assist the QMMM input preparation.
# by  Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#     Abigail Chiu, B.S. (abigail.chiu@ucdenver.edu)  
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

poky_qmmm='''
import os
from pymol import cmd
cmd.delete('pdb_place_holder')
cmd.load('pdb_place_holder')

def find_serial_number(object_name, cap=False):
  atom_list = cmd.get_model(object_name).atom
  if not cap:
    return list(map(lambda atm: '%9d' % atm.index, atom_list))
  return list(map(lambda atm: '%10d %4s    0.731' % (atm.index, atm.name), 
                              atom_list)) 

def parse_serial_number(lines, begin_index, end_index):
  result_list = []
  for line in lines[begin_index+1:end_index]:
    sp_list = line.split()
    if len(sp_list) == 0:
      continue
    result_list.append(sp_list[0].strip())
  return result_list

def find_qmmm_section(lines, section, start_index=0):
  for i in range(start_index, len(lines)):
    line = lines[i].strip()
    if line == section:
      index = i
      return index
  return -1

def insert_qmmm_lines(lines, dest_lines, top_bound, bottom_bound):
  begin_index = find_qmmm_section(lines, top_bound)
  end_index = find_qmmm_section(lines, bottom_bound, begin_index+1)
  if -1 not in [begin_index, end_index]:
    return lines[:begin_index+1] + dest_lines + lines[end_index:]
  return None

def export_qmmm(qm_object_name='QMATOM', cap_object_name='CAPATOM', 
                dat_file_name='dat_place_holder'):
  
  qmatoms = find_serial_number(qm_object_name)
  capatoms = find_serial_number(cap_object_name, True)

  if os.path.exists(dat_file_name):
    f = open(dat_file_name, 'r')
    dat_lines = f.readlines()
    f.close()
    
    
    idx_qmatm = find_qmmm_section(dat_lines, 'QMATOM')
    idx_capatm = find_qmmm_section(dat_lines, 'CAPATOM')
    idx_embed = find_qmmm_section(dat_lines, 'EMBED mechanical')
    idx_qmmm = find_qmmm_section(dat_lines, '*QM/MM')
      
    # CAPATOM section exists 
    if idx_capatm != -1:
        dat_lines = insert_qmmm_lines(dat_lines, capatoms, 'CAPATOM', 'END')
    # EMBED mechanical exists
    elif idx_embed != -1:
        dat_lines.insert(idx_embed+1, '')
        dat_lines.insert(idx_embed+2, '    CAPATOM')
        dat_lines.insert(idx_embed+3, '    END')
        dat_lines = insert_qmmm_lines(dat_lines, capatoms, 'CAPATOM', 'END')
    # *QM/MM section exists
    elif idx_qmmm != -1:
        dat_lines.insert(idx_qmmm+1, '')
        dat_lines.insert(idx_qmmm+2, '    CAPATOM')
        dat_lines.insert(idx_qmmm+3, '    END')
        dat_lines = insert_qmmm_lines(dat_lines, capatoms, 'CAPATOM', 'END')
    # Not a good dat file. Complaint and quit.
    else:
        print('ERROR! .DAT file could not be created.')
        
    # QMATOM section exists                
    if idx_qmatm != -1:
        dat_lines = insert_qmmm_lines(dat_lines, qmatoms, 'QMATOM', 'END')
    # EMBED mechanical exists
    elif idx_embed != -1:
        dat_lines.insert(idx_embed+1, '')
        dat_lines.insert(idx_embed+2, '    QMATOM')
        dat_lines.insert(idx_embed+3, '    END')
        dat_lines = insert_qmmm_lines(dat_lines, qmatoms, 'QMATOM', 'END')
    # *QM/MM section exists
    elif idx_qmmm != -1:
        dat_lines.insert(idx_qmmm+1, '')
        dat_lines.insert(idx_qmmm+2, '    QMATOM')
        dat_lines.insert(idx_qmmm+3, '    END')
        dat_lines = insert_qmmm_lines(dat_lines, qmatoms, 'QMATOM', 'END')
    # Not a good dat file. Complaint and quit.
    else:
        print('ERROR! .DAT file could not be created.')
    
  else:
    dat_lines = ['    QMATOM',] + qmatoms + ['    END','']
    dat_lines = dat_lines + ['    CAPATOM',] + capatoms + ['    END']

  # Remove right null characters
  dat_lines = list(map(lambda x: x.rstrip(), dat_lines))
  
  f = open(dat_file_name, 'w')
  f.write('\\n'.join(dat_lines))
  f.close()

  print(dat_file_name + ' saved.')


def view_qmmm(qm_object_name='QMATOM', cap_object_name='CAPATOM', 
                dat_file_name='dat_place_holder', color='red'):
  
  if not os.path.exists(dat_file_name):
    print(f'ERROR: {dat_file_name} does not exist.')
    return
  f = open(dat_file_name, 'r')
  dat_lines = f.readlines()
  f.close()
  
  qm_begin_index = find_qmmm_section(dat_lines, 'QMATOM')
  qm_end_index = find_qmmm_section(dat_lines, 'END', qm_begin_index)
  qm_list = parse_serial_number(dat_lines, qm_begin_index, qm_end_index)

  cap_begin_index = find_qmmm_section(dat_lines, 'CAPATOM')
  cap_end_index = find_qmmm_section(dat_lines, 'END', cap_begin_index)
  cap_list = parse_serial_number(dat_lines, cap_begin_index, cap_end_index)

  all_atom_list = qm_list + cap_list
  if len(all_atom_list) == 0:
    print('ERROR: NO QMATOM AND CAPATOM FOUND.')
    return
  
  qm_atoms = '+'.join(qm_list)
  cmd.select('QMATOM', 'idx. ' + qm_atoms)
  cmd.color(color, 'QMATOM')
  cmd.show('sphere', 'QMATOM')
  if len(cap_list) != 0:
    cap_atoms = '+'.join(cap_list)
    cmd.select('CAPATOM', 'idx. ' + cap_atoms)
    cmd.color(color, 'CAPATOM')
    cmd.show('sphere', 'CAPATOM')
  cmd.zoom('QMATOM')  
'''

import __main__
s = __main__.main_session

import os
temp_dir = os.path.dirname(__file__)

# Select a PDB File
pdb_file = s.open_filedialog('Select the PDB file', 
                  'PDB (*.pdb);; mmCIF (*.cif);; Any (*)', '')
if pdb_file == '':
  raise SystemError
poky_qmmm = poky_qmmm.replace('pdb_place_holder', pdb_file)

# Select a DAT File
dat_file = s.open_filedialog('Select the DAT file', 
                  'DAT (*.dat);; Any (*)', '')
if dat_file == '':
  dat_file = s.save_filedialog('Save the DAT file as ...', 
                  'DAT (*.dat);; Any (*)', '')
if dat_file == '':
  raise SystemError
poky_qmmm = poky_qmmm.replace('dat_place_holder', dat_file)

f = open(os.path.join(temp_dir, 'poky_qmmm.py'), 'w')
f.write(poky_qmmm)
f.close()

clipboard_text = f'''
sys.path.append('{temp_dir}'); from poky_qmmm import export_qmmm as qmmm; from poky_qmmm import view_qmmm
'''

s.set_clipboard(clipboard_text)
s.open_pymol('')
s.show_message('PyMOL', 'CTRL+V in the PyMOL commandline. \n'
    'Select QM atoms and rename the object name to QMATOM. \n'
    'Select CAP atoms and rename the object name to CAPATOM. \n'
    'Then, type qmmm() to save the DAT file. \n'
    'You can type view_qmmm() to view QM/CAP atoms in PyMOL.\n'
    'Don\'t close this window until you run qmmm(). ')

s.open_notepad(dat_file)
