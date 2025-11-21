#
# This is an example script for extracting peaks from nef.
#
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Nov 21, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

# POKY libraries
import __main__
s = __main__.main_session

import os
import io
import pandas as pd

# Dictionary to map 3-letter amino acid codes to 1-letter codes
AMINO_ACID_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
def convert_three_to_one_letter(three_letter_code):
    """Converts a three-letter amino acid code to its one-letter code."""
    return AMINO_ACID_CODES.get(three_letter_code.upper(), three_letter_code)

def extract_nef_peak_block(filename):
    if not os.path.exists(filename):
        s.show_message('Error', f"Error: File not found at path: {filename}")
        return None

    f = open(filename, 'r')
    lines = f.readlines()
    lines = list(map(lambda x: x.strip(), lines))
    f.close()

    # check number of spectrum
    spectrum_list = list(filter(lambda x: x.startswith('save_nef_nmr_spectrum_'), lines))
    specname_list = list(map(lambda x: x.replace('save_nef_nmr_spectrum_', ''), spectrum_list))
    if len(spectrum_list) == 0:
        return None
    
    if len(specname_list) > 1:
        ans = s.show_customselectiondialog('Select a spectrum', specname_list, 0)
        if ans == '':
            return None
        chosen_spec = 'save_nef_nmr_spectrum_' + ans
    else:
        chosen_spec = spectrum_list[0]
    
    lines = lines[lines.index(chosen_spec):]

    start, end = -1, -1
    for i, line in enumerate(lines):
        if line.strip().startswith('_nef_peak.index'):
            start = i-1
        if start != -1 and line.strip() == 'stop_':
            end = i
            break
    if -1 in [start, end]:
        return None
  
    return '\n'.join(lines[start+1:end])


def generate_poky_list_nd(nef_data):
    """
    Converts N-dimensional NEF formatted peak data into a poky peak list format.
    (This is the robust function from the previous step, unchanged).
    """
    # ... [The complete logic for parsing and formatting from the previous answer goes here] ...
    # 1. Prepare and parse the data
    lines = [line.strip() for line in nef_data.split('\n') if line.strip()]
    header_lines = [line for line in lines if line.startswith('_nef_peak.')]
    data_lines = [line for line in lines if not line.startswith('_nef_peak.')]

    # Map NEF column names to simpler names
    columns = [col.replace('_nef_peak.', '') for col in header_lines]
    data_str = '\n'.join(['\t'.join(columns)] + data_lines)

    try:
        df = pd.read_csv(io.StringIO(data_str), sep=r'\s+', engine='python')
    except Exception as e:
        return f"Error parsing data with pandas: {e}"

    # 2. Determine Peak Dimensionality and Extract Data
    max_dim = 2
    data_columns = []
    assignment_parts = []
    
    # Check for dimensions up to 4
    for dim in range(1, 5):
        pos_col = f'position_{dim}'
        res_col = f'residue_name_{dim}'
        seq_col = f'sequence_code_{dim}'
        atom_col = f'atom_name_{dim}'
        
        if pos_col in df.columns:
            max_dim = dim
            
            # Extract chemical shift (wN) and round it
            try:
                w_col = df[pos_col].astype(float).round(3)
                data_columns.append(w_col)
            except ValueError:
                 print(f"Warning: Position column {dim} contains non-numeric data. Skipping this dimension.")
                 continue

            # Extract assignment parts for this dimension
            if all(col in df.columns for col in [res_col, seq_col, atom_col]):
                res_code = df[res_col].apply(convert_three_to_one_letter)
                seq_code = df[seq_col].astype(str)
                atom_name = df[atom_col].astype(str)
                dim_assignment = res_code + seq_code + atom_name
                assignment_parts.append(dim_assignment)
            else:
                print(f"Warning: Assignment columns for dimension {dim} are missing. Using Peak ID instead.")
                assignment_parts.append(df['peak_id'].astype(str))
        else:
            break
            
    if max_dim < 2:
        return "Error: Could not find required 'position_1' and 'position_2' columns."

    # 3. Construct the poky Output
    final_assignment = assignment_parts[0]
    for part in assignment_parts[1:]:
        final_assignment += '-' + part
        
    poky_data = {'Assignment': final_assignment}
    for i, col in enumerate(data_columns):
        poky_data[f'w{i+1}'] = col

    poky_df = pd.DataFrame(poky_data)
    poky_output = poky_df.to_string(header=True, index=False)
    lines = poky_output.split('\n')
    lines.insert(1, '')
    poky_output = '\n'.join(lines)
    return poky_output

nef_file = s.open_filedialog('Select an NEF peak list file.', 'NEF file (*.nef)', '')
if nef_file == '':
    raise SystemError

nef_data_block = extract_nef_peak_block(nef_file)

if nef_data_block:
    poky_output = generate_poky_list_nd(nef_data_block)
    list_file = s.save_filedialog('Save as', 'POKY peak list file (*.list)', nef_file[:-4] + '.list')
    if list_file != '':
        with open(list_file, 'w') as f:
            f.write(poky_output)
    s.show_message('Done.', f'{list_file} created. Read peak list using "rp".')        
else:
    s.show_message('Failed.', 'No peaks.')
    
