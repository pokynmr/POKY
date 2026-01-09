#
# This is an example script for removing specific lines with a keyword in a file.
#
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: January 9, 2026
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

# POKY libraries
import __main__
s = __main__.main_session

in_file = s.open_filesdialog('Select a file', 'Any (*)', '')
if in_file == '':
  raise SystemError

keyword = s.show_inputdialog('Keyword', 'Type a keyword to detect', '999')
if keyword == '':
  raise SystemError

out_file = s.save_filesdialog('Save as', 'Any (*)', '')
if out_file == '':
  raise SystemError

with open(in_file, 'r') as f_in:
  lines = f_in.readlines()

with open(out_file, 'w') as f_out:
  for line in lines:
    if keyword not in line:
        f_out.write(line)

s.show_message('Done.', f'File saved: {out_file}')