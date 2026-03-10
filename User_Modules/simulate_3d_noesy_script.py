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

def calculate_r6(source_coords_list, target_coords_list):
    all_dist_r6 = []
    for s_coords, t_coords in zip(source_coords_list, target_coords_list):
        for sc in s_coords:
            for tc in t_coords:
                dist = max(np.linalg.norm(sc - tc), 0.1)
                all_dist_r6.append(dist**-6)
    return np.mean(all_dist_r6)**(-1/6) if all_dist_r6 else 999.0

def run_noesy_full(s):
    # 1. Setup
    dist_input = s.show_inputdialog('Full NOE Prediction', 'Distance Cutoff (A):', '5.0')
    if not dist_input: return
    dist_limit = float(dist_input)

    file_path = s.open_filedialog('Select Structure', 'Structure (*.pdb *.cif);;All (*)', '')
    if not file_path: return

    obj_name = f"struct_{int(time.time())}"
    cmd.delete(obj_name)
    cmd.load(file_path, obj_name)
    
    # 2. Handle Models
    num_states = cmd.count_states(obj_name)
    states_to_process = []
    if num_states > 1:
        is_ensemble = s.show_message_yes_no('Mode Selection', f'Found {num_states} models. Use Ensemble Mode?')
        states_to_process = list(range(1, num_states + 1)) if is_ensemble else [1]
    else:
        states_to_process = [1]

    if cmd.count_atoms(f"{obj_name} and symbol H") == 0:
        cmd.h_add(obj_name)

    # 3. Resonance Data
    condition = s.project.condition_list()[0]
    poky_data = {str(r.group.name).strip(): {} for r in condition.resonance_list()}
    for r in condition.resonance_list():
        poky_data[str(r.group.name).strip()][str(r.atom.name).strip()] = r.frequency

    # 4. Global Proton Collection
    # We collect every proton that exists in Poky, regardless of residue
    ensemble_coords = {} 
    site_meta = {}

    for i, state in enumerate(states_to_process):
        model = cmd.get_model(obj_name, state=state)
        for h in [a for a in model.atom if a.symbol.upper() == 'H']:
            res_num = str(int(h.resi)).strip()
            group_id = f"{AA_MAP.get(h.resn.upper().strip(), '')}{res_num}"
            
            if group_id not in poky_data: continue

            # Map PDB H-name to Poky H-name
            h_pdb = h.name.strip()
            target_h = None
            if h_pdb in poky_data[group_id]: target_h = h_pdb
            elif h_pdb == "H" and "HN" in poky_data[group_id]: target_h = "HN"
            else:
                prefix = ''.join([c for c in h_pdb if not c.isdigit()])
                for p_atom in poky_data[group_id]:
                    if p_atom.startswith(prefix) and "H" in p_atom:
                        target_h = p_atom; break
            
            if target_h:
                key = (group_id, target_h)
                if key not in ensemble_coords: ensemble_coords[key] = [[] for _ in states_to_process]
                ensemble_coords[key][i].append(np.array(h.coord))
                
                # Heavy atom mapping
                if key not in site_meta:
                    for a in model.atom:
                        if a.symbol.upper() in ['N', 'C'] and a.resi == h.resi:
                            if np.linalg.norm(np.array(h.coord) - np.array(a.coord)) < 2.0:
                                site_meta[key] = {'heavy': a.name.strip(), 'symbol': a.symbol.upper(), 'res': group_id}
                                break

    # 5. Build Final Sites
    final_sites = []
    for key, coords_list in ensemble_coords.items():
        if key not in site_meta: continue
        group_id, h_name = key
        meta = site_meta[key]
        w2 = poky_data[group_id][h_name]
        w1 = poky_data[group_id].get(meta['heavy']) or poky_data[group_id].get(meta['heavy'] + "_s")
        
        if w2 != 0.0 and w1 is not None:
            final_sites.append({'group': group_id, 'h_name': h_name, 'coords': coords_list, 'w1': w1, 'w2': w2, 'symbol': meta['symbol'], 'heavy_name': meta['heavy']})

    # 6. Global Comparison (Intra + Inter)
    lists = {'15N': [], '13C': []}
    for i, s1 in enumerate(final_sites):
        for j, s2 in enumerate(final_sites):
            if i == j: continue # Allow i != j to get full cross-peaks
            
            dist = calculate_r6(s1['coords'], s2['coords'])
            if dist <= dist_limit:
                # Label: SourceRes-SourceHeavy-SourceH-TargetH (e.g., Q2N-H-V17HA)
                # To match standard Poky format: Assignment w1 w2 w3 Dist
                label = f"{s1['group']}{s1['heavy_name']}-{s1['h_name']}-{s2['group']}{s2['h_name']}"
                line = f"{label:26} {s1['w1']:8.3f} {s1['w2']:8.3f} {s2['w2']:8.3f} {dist:6.2f}"
                
                if s1['symbol'] == 'N': lists['15N'].append(line)
                else: lists['13C'].append(line)

    # 7. Save
    out_dir = os.path.expanduser("~")
    for k, content in lists.items():
        if content:
            path = s.save_filedialog(f'Save {k} Full List', 'Poky List (*.list)', os.path.join(out_dir, f'full_noe_{k}.list'))
            if path:
                with open(path, 'w') as f:
                    f.write("Assignment                  w1       w2       w3       Dist\n\n" + "\n".join(content))

    cmd.delete(obj_name)
    s.show_message('Success', f"Full prediction done! Generated {sum(len(v) for v in lists.values())} peaks.")

if __name__.startswith('pyrun_'):
    import __main__
    run_noesy_full(__main__.main_session)
