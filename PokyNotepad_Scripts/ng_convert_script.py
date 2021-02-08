#
# This is an example script to convert NMRPIPE file to POKY.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Using NMRGlue.
#
# If the data is large, it will freeze for a few minutes!! Save before running this script.
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module


# You can drag and drop into this notepad text area to get file path easy.
# You must delete "file://" in the beginning of the dropped text though.
input_files = '/path/to/test%03d.ft3'   # 3D data
#input_files = '/path/to/test.ft2'      # 2D data

output_file = '/path/to/save/file.ucsf'

import nmrglue as ng
dic,data = ng.pipe.read_lowmem(input_files)

# Set the parameters
u = ng.pipe.guess_udic(dic,data)

# create the converter object and initilize with varian data
C = ng.convert.converter()
C.from_pipe(dic,data,u)

# create pipe data and then write it out
ng.sparky.write(output_file,*C.to_sparky(),overwrite=True)

print('Done.')
