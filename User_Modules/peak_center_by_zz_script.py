#
# This is an example script for adding two-letter-code 'zz' for peak centering.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
# 

import __main__
s = __main__.session

def center_selected():
  for p in s.selected_peaks():
    p.center()

s.add_command('zz', 'Center selected peaks', center_selected)