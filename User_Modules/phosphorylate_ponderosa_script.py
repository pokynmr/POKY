#
# This is an example script for phosphorylating PONDEROSA results.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: Jan. 14, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

######################
# USER PARAMETER START
######################

# Ponderosa Job ID from the received URL.
jobid = '250101_010101_001' 

# Residues to phosphorylate. S, T, Y are supported. 
residues = '3,5,6'

######################
# USER PARAMETER END
######################

# POKY libraries
import __main__
import serverside
s = __main__.main_session

residues = residues.replace(' ', '')
fields = {'ponderosaid': jobid, 'pres': residues}
serverside_result = serverside.send('PHOSPHOPOND', files=[], parameters=fields)

s.show_message('Submitted', 'Job submitted. Check your email inbox.')