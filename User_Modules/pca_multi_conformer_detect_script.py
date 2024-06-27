#
# This is an example script to detect multi conformers in a PDB ensemble.
# This script runs on POKY BUILD 11/27/23q or newer
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# You can demonstrate this with PDB 7SA5.
#       https://www.rcsb.org/structure/7SA5
#

import __main__
s = __main__.main_session

# test with 7sa5
import numpy as np

try:
    from pokypdb import modelPDB, readPDB, createContactMap
    from pokydihe import get_all_phi_psi_by_file
    from wlutil import simple_pca, simple_cluster, simple_pca_plot
except:
    s.show_message('Error', 'Please update POKY to use this script.')
    s.open_url('https://poky.clas.ucdenver.edu')
    raise SystemError


ans = s.show_selectionexdialog('Approach', 'Choose a method', 
            ('Distance', 'Angle with geoPCA',  'Angle with dPCA', 'Cancel'))

if ans == 3:
    raise SystemError

pdb_file = s.open_filedialog('Select a PDB file.',
                             'PDB file (*.pdb);; CIF file (*.cif)', '')

if pdb_file == '' or \
        not (pdb_file.endswith('.pdb') or pdb_file.endswith('.cif')):
    print("You didn't select a PDB file")
    raise SystemExit

if pdb_file.endswith('.cif'):
    from wlutil import get_temp_file_name
    from cif2pdb import cif2pdb
    pdb_file2 = get_temp_file_name('.pdb')
    cif2pdb(pdb_file, pdb_file2)
else:
    pdb_file2 = pdb_file


# interdistance approach
nmodel = modelPDB(pdb_file2)
contact_type = 'CA'

ncluster = int(s.show_inputdialog('Clusters', 
                              '# of conformations: ', '2'))

if ans == 0:
    cm_list = []
    for i in range(1, nmodel+1):
        pdb_list = readPDB(pdb_file2, i)
        contact_map, residue_range = createContactMap(pdb_list,
                                                contact_type)
        cm_list.append(np.array(contact_map))
        try:
            data = np.vstack((data, cm_list[-1].flatten()))
        except:
            data = cm_list[-1].flatten()

# phi, psi approach- 
# 1: polar to cartesian and then pca
# 2: phi -> cos(phi), sin(phi).... then pca
elif ans in [1, 2]:
    phi_angles, psi_angles = get_all_phi_psi_by_file(pdb_file2)
    for i in range(1, nmodel+1):
        if ans == 1:
            x = np.multiply(np.cos(phi_angles), np.cos(psi_angles))
            y = np.multiply(np.cos(phi_angles), np.sin(psi_angles))
            z = np.sin(phi_angles)
            xyz = np.hstack((x, y, z))
        else:
            xyz = np.hstack((np.cos(phi_angles), np.sin(phi_angles),
                             np.cos(psi_angles), np.sin(psi_angles)))
        try:
            data = np.vstack((data, xyz))
        except:
            data = xyz

converted_data, loadings, explained = simple_pca(data)
print('---------------------------------')
print('  POKY MULTICONFORMER DETECTION  ')
print('---------------------------------')
print(' MODEL PC1      PC2      CLUSTER ')
print('---------------------------------')
cluster_assignments = simple_cluster(converted_data, ncluster=ncluster)
colors = {1: 'blue', 2: 'green', 3: 'purple', 4: 'magenta', 5: 'orange', 
        6: 'pink', 7: 'cyan', 8: 'yellow', 9: 'gray', 10: 'lightblue'} 
plt = simple_pca_plot(converted_data, loadings, explained, cluster_assignments,
                      ncluster=ncluster, colors=colors)

if s.show_message_yes_no('PyMOL', 'Do you want to open PyMOL to view ' + \
    'clustered conformers? You can simply CTRL+V in the PyMOL commandline.'):
    # make pymol text
    cmd = f'delete pokymc; load {pdb_file2}, pokymc; split_states pokymc\n'
    for i in range(ncluster):
        filtered = list(filter(lambda x: cluster_assignments[x-1] == i+1, 
                    range(1, len(cluster_assignments)+1)))
        models = list(map(lambda x: f'pokymc_{x:04}', filtered))
        m = '|'.join(models)
        cmd += f'join_states pokymc{i+1}, {m}, 0; '
        cmd += f'color {colors[i+1]}, pokymc{i+1};'
        cmd += f'intra_fit pokymc{i+1}///CA, 1\n'
    cmd += 'color gray, pokymc; delete pokymc_*\n'
    cmd += 'delete pokymc; set all_state, on; bg white\n'
    s.set_clipboard(cmd)
    s.open_pymol('')