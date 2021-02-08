#
# This is an example script to process Bruker 2D data and convert to POKY.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# Using NMRGlue. Read this reference guide:
#   https://nmrglue.readthedocs.io/en/latest/reference/proc_base.html
#
#
# If the data is large, it will freeze for a few minutes!! Save before running this script.
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module


magnet = 'varian' # bruker or varian
LoadData = True  # Load data in POKY


######################################################
# You will end up turning knobs for phasing down there
import nmrglue as ng

# POKY libraries
import __main__
s = __main__.main_session

input_path = s.open_directorydialog('Select your FID directory', '')
output_file = s.save_filedialog('Save your new .ucsf file as...',
                                'UCSF NMR (*.ucsf);; Any (*)',
                                input_path)

C = ng.convert.converter()
if magnet == 'bruker':
  dic,data = ng.bruker.read(input_path)
elif magnet == 'varian':
  dic,data = ng.varian.read(input_path)

# direct dimension
data = ng.proc_base.sp(data, off=0.45, end=0.98, pow=2)
data = ng.proc_base.zf(data)
data = ng.proc_base.fft(data)
data = ng.proc_base.ps(data, p0=-145.4, p1=159.0)
data = ng.proc_base.di(data)

# indirect dimension
data = ng.proc_base.decode_States(data)   # acq2D
data = ng.proc_base.tp(data)
data = ng.proc_base.sp(data, off=0.45, end=0.9, pow=1)
data = ng.proc_base.zf(data)
data = ng.proc_base.ifft(data)
data = ng.proc_base.ps(data, p0=0.0, p1=0.0)
data = ng.proc_base.di(data)

# polynomial baseline correction has not been implemented in NMRGlue yet
# median baseline correction is used instead
data = ng.proc_bl.med(data)
data = ng.proc_base.tp(data)
data = ng.proc_bl.med(data)

# Converter
if magnet == 'bruker':
  C.from_bruker(dic, data)
elif magnet == 'varian':
  C.from_varian(dic, data)
# create pipe data and then write it out
ng.sparky.write(output_file,*C.to_sparky(),overwrite=True)

if LoadData:
  s.open_spectrum(output_file)

print('Done.')

