#
# This is an example script for making the data height absolute 
# bruker spectrum.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: April 21, 2022
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# POKY 01/14/22e or higher is required to run this script.
#

import numpy as np
# POKY libraries
import __main__
s = __main__.main_session

if s.show_message_yes_no('Warning', 
    'This script will overwrite original data. Okay?') == 0:
  raise SystemError

org_file = s.open_filedialog('Select a spectrum', 
                '2rr (2rr);; 1r (1r);; 3rrr (3rrr);; Any (*)', '')

if org_file == '':
  raise SystemError

import os.path
org_path = os.path.dirname(org_file)

import nmrglue as ng
dic, data = ng.bruker.read_pdata(org_path)
new_data = np.absolute(data)

ng.bruker.write_pdata(org_path, dic, new_data, overwrite=True)

if s.show_message_yes_no('Load data', 'Do you want to load the processed spectrum?'):
  s.open_spectrum(org_file)
