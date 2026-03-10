#
# This is an example script to simulate 3D NOESY using a PDB and resonances.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#


from pymol import cmd
import numpy as np
import sputil
import time
import os

AA_MAP = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M', 'PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def get_poky_atom_name(heavy_atom_name, pdb_h_name, res_poky_data):
    """
    Attempts to map a generic H name (H01) to a Poky name based on its heavy parent.
    Example: Attached to CA -> HA. Attached to N -> HN.
    """
    # 1. Determine the 'Type' based on the heavy atom
    # Removing numbers from heavy atom (CA1 -> CA, CB -> CB)
    h_type = ''.join([i for i in heavy_atom_name if not i.isdigit()])
    
    # Map Heavy type to Proton Prefix
    prefix_map = {'N': 'H', 'CA': 'HA', 'CB': 'HB', 'CG': 'HG', 'CD': 'HD', 'CE': 'HE', 'CZ': 'HZ'}
    prefix = prefix_map.get(h_type, 'H')
    
    # 2. Try to find a match in Poky data for this residue
    # Check for direct prefix (HN, HA)
    if prefix == 'H' and 'HN' in res_poky_data: return 'HN'
    if prefix in res_poky_data: return prefix
    
    # Check for numbered versions (HA1, HA2, etc.)
    potential_matches = [a for a in res_poky_data.keys() if a.startswith(prefix)]
    if potential_matches:
        # If there are multiple (like HB1, HB2), we'd ideally match by PDB index
        # but for now, we return the first available match
        return potential_matches[0]
        
    return None

def calculate_r6(source_coords_list, target_coords_list):
    all_dist_r6 = []
    for s_coords, t_coords in zip(source_coords_list, target_coords_list):
        for sc in s_coords:
            for tc in t_coords:
                dist = max(np.linalg.norm(sc - tc), 0.1)
                all_dist_r6.append(dist**-6)
    return np.mean(all_dist_r6)**(-1/6) if all_dist_r6 else 999.0

def run_noesy_generic_fix(s):
    dist_input = s.show_inputdialog('NOE Prediction', 'Distance Cutoff (A):', '5.0')
    if not dist_input: return
    dist_limit = float(dist_input)

    file_path = s.open_filedialog('Select Structure', 'Structure (*.pdb *.cif);;All (*)', '')
    if not file_path: return

    obj_name = f"struct_{int(time.time())}"
    cmd.delete(obj_name)
    cmd.load(file_path, obj_name)
    
    # Standardize hydrogens
    # 'h_add' with 'infer=1' tries to name them better, but if it fails, we map manually
    cmd.h_add(obj_name)

    # 1. Index Poky Resonances
    condition = s.project.condition_list()[0]
    poky_data = {str(r.group.name).strip(): {} for r in condition.resonance_list()}
    for r in condition.resonance_list():
        poky_data[str(r.group.name).strip()][str(r.atom.name).strip()] = r.frequency

    # 2. Extract and Geometric Mapping
    num_states = cmd.count_states(obj_name)
    states_to_process = list(range(1, num_states + 1)) if num_states > 1 else [1]
    
    ensemble_coords = {} 
    site_meta = {}

    for i, state in enumerate(states_to_process):
        model = cmd.get_model(obj_name, state=state)
        # Separate heavy atoms and protons
        heavy_atoms = [a for a in model.atom if a.symbol.upper() in ['N', 'C']]
        protons = [a for a in model.atom if a.symbol.upper() == 'H']
        
        for h in protons:
            res_num = str(int(h.resi)).strip()
            group_id = f"{AA_MAP.get(h.resn.upper().strip(), '')}{res_num}"
            if group_id not in poky_data: continue

            # Find heavy atom parent by distance (< 1.2 A)
            h_coord = np.array(h.coord)
            parent = None
            for a in heavy_atoms:
                if a.resi == h.resi and np.linalg.norm(h_coord - np.array(a.coord)) < 1.3:
                    parent = a; break
            
            if not parent: continue
            
            # Map Generic Name (H01) to Poky Name (HA, HN, etc.)
            target_h = get_poky_atom_name(parent.name.strip(), h.name.strip(), poky_data[group_id])
            
            if target_h:
                key = (group_id, target_h)
                if key not in ensemble_coords: ensemble_coords[key] = [[] for _ in states_to_process]
                ensemble_coords[key][i].append(h_coord)
                
                if key not in site_meta:
                    site_meta[key] = {'heavy': parent.name.strip(), 'symbol': parent.symbol.upper()}

    # 3. Process Peaks
    final_sites = []
    for (group_id, h_name), coords_list in ensemble_coords.items():
        if (group_id, h_name) not in site_meta: continue
        meta = site_meta[(group_id, h_name)]
        w2 = poky_data[group_id][h_name]
        w1 = poky_data[group_id].get(meta['heavy']) or poky_data[group_id].get(meta['heavy'] + "_s")
        
        if w2 != 0.0 and w1 is not None:
            final_sites.append({'group': group_id, 'h_name': h_name, 'coords': coords_list, 'w1': w1, 'w2': w2, 'symbol': meta['symbol'], 'heavy': meta['heavy']})

    lists = {'15N': [], '13C': []}
    for i, s1 in enumerate(final_sites):
        for j, s2 in enumerate(final_sites):
            if i == j: continue
            dist = calculate_r6(s1['coords'], s2['coords'])
            if dist <= dist_limit:
                label = f"{s1['group']}{s1['heavy']}-{s1['h_name']}-{s2['group']}{s2['h_name']}"
                line = f"{label:28} {s1['w1']:8.3f} {s1['w2']:8.3f} {s2['w2']:8.3f} {dist:6.2f}"
                lists['15N' if s1['symbol'] == 'N' else '13C'].append(line)

    # 4. Save
    out_dir = os.path.expanduser("~")
    for k, content in lists.items():
        if content:
            path = s.save_filedialog(f'Save {k}', 'Sparky List (*.list)', os.path.join(out_dir, f'mapped_{k}.list'))
            if path:
                with open(path, 'w') as f:
                    f.write("Assignment                  w1       w2       w3       Dist\n\n" + "\n".join(content))

    cmd.delete(obj_name)
    s.show_message('Success', f"Mapped generic H-names to {len(final_sites)} Poky resonances.")

if __name__.startswith('pyrun_'):
    import __main__
    run_noesy_generic_fix(__main__.main_session)
