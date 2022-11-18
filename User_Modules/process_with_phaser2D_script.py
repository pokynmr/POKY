#
# This is an example script to process Varian/Bruker data and convert to POKY.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
#
# If the data is large, it will freeze for a few minutes!! 
# Save before running this script.
# Basically, this script works on 2D data. 
# However, PPM won't be shown.  
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

##################################################################
# USER PARAMETER START
##################################################################
LoadData = True       # Load data in POKY
DigitalFilter = True  # Remove digital filter (Bruker)
# Encodings: None, 'undefined', 'magnitude', 
# 'tppi', 'states', 'states-tppi', 'echo-antiecho'
w1_encoding, w2_encoding, w3_encoding = None, None, None
# Baseline correction: (T/F, T/F, T/F)
median_baseline_correction = (False, False, False)
polyfit_baseline_correction = (True, True, True)
# Zero-filling: (True/False, 'auto' or size)
w1_zf, w2_zf, w3_zf = (True, 1024), (True, 1024), (True, 1024)
# Linear prediction
w1_lp, w2_lp, w3_lp = False, False, False
# Reverse data
w1_reverse, w2_reverse, w3_reverse = False, False, False
# Extract: None, 'left', 'right'
w1_ext, w2_ext, w3_ext = None, None, None
# Apodization SP
w1_off, w1_end, w1_pow, w1_c = 0.45, 0.98, 1, 1.0
w2_off, w2_end, w2_pow, w2_c = 0.45, 0.98, 2, 0.5
w3_off, w3_end, w3_pow, w3_c = 0.45, 0.98, 2, 0.5
# Change below once found. Use None to activate Poky Phaser 
# Also, 'auto' can be used to run acme auto-phasing instead.
w1_p0, w1_p1 = None, None
w2_p0, w2_p1 = None, None
w3_p0, w3_p1 = None, None

##################################################################
# USER PARAMETER END
##################################################################


encoding_list = [w1_encoding, w2_encoding, w3_encoding]
zf_list = [w1_zf, w2_zf, w3_zf]
lp_list = [w1_lp, w2_lp, w3_lp]
reverse_list = [w1_reverse, w2_reverse, w3_reverse]
ext_list = [w1_ext, w2_ext, w3_ext]
sp_list = [ {'off': w1_off, 'end': w1_end, 'pow': w1_pow, 'c': w1_c}, 
            {'off': w2_off, 'end': w2_end, 'pow': w2_pow, 'c': w2_c},
            {'off': w3_off, 'end': w3_end, 'pow': w3_pow, 'c': w3_c}]
p_list = [[w1_p0, w1_p1], [w2_p0, w2_p1], [w3_p0, w3_p1]]

import os
import nmrglue as ng
import numpy as np
    
# POKY libraries
import __main__
s = __main__.main_session
try:
  import pokyphaser
  from pokyphaser import show_phaser2D
except:
  s.show_message('Error', 'Please update POKY to use this module.\n' + \
              'Visit https://poky.clas.ucdenver.edu')
  raise SystemError

input_path = s.open_directorydialog('Select your FID directory', '')
ser_path = os.path.join(input_path, 'ser')
fid_path = os.path.join(input_path, 'fid')
procpar_path = os.path.join(input_path, 'procpar')

C = ng.convert.converter()
if os.path.exists(ser_path) or os.path.exists(fid_path):
  magnet = 'bruker'
  dic,data = ng.bruker.read(input_path)
  if DigitalFilter:
    data = ng.bruker.remove_digital_filter(dic, data)
  u = ng.bruker.guess_udic(dic, data)
  C.from_bruker(dic, data, u)
elif os.path.exists(procpar_path):
  magnet = 'varian'
  dic,data = ng.varian.read(input_path)
  u = ng.varian.guess_udic(dic, data)
  C.from_varian(dic, data, u)
else:
  print('Only Bruker/Varian is supported.')
  raise SystemError
ndim = len(data.shape)

for dim in range(ndim):
  if encoding_list[dim] != None:
    u[dim]['encoding'] = encoding_list[dim]
dic, data = C.to_pipe()
# Processing Start
# Checking intermediate data:
#   pokyphaser.plot_data(data, show_as=1)
# Don't use plot_data with pokyphaser.show_phaser(s)
# You can stop running this code before that by adding
#   raise SystemError
for dim in range(ndim):
  # Time Domain Correction
  if dim == 0:
    try:
      data = pokyphaser.td_bl(data)
    except:
      import numpy as np
      data -= np.mean(np.atleast_2d(data)[..., int(data.shape[-1] / -4.):])
  # Apodization
  dic, data = ng.pipe_proc.sp(dic, data, 
    off=sp_list[dim]['off'], end=sp_list[dim]['end'], 
    pow=sp_list[dim]['pow'], c=sp_list[dim]['c'])
  # Zero-Filling
  if zf_list[dim][0]:
    if zf_list[dim][1] == 'auto':
      dic, data = ng.pipe_proc.zf(dic, data, auto=True)
    else:
      dic, data = ng.pipe_proc.zf(dic, data, size=zf_list[dim][1])
  # Linear Prediction
  if lp_list[dim]:
    dic, data = ng.pipe_proc.lp(dic, data)
  # Fourier Transform
  dic, data = ng.pipe_proc.ft(dic, data, auto=True)
  # Phasing
  if 'auto' in p_list[dim]:
    data, opt = ng.proc_autophase.autops(data, 'acme', return_phases=True)
    print(f'w{dim} acme: p0->{opt[0]}')
    print(f'w{dim} acme: p1->{opt[1]}')    
  elif None not in p_list[dim]:
    data = ng.proc_base.ps(data, p0=p_list[dim][0], p1=p_list[dim][1])
  else:
    noise = pokyphaser.get_noise(data.real)
    s.show_message(f'w{dim+1} Phasing', 
        f'Replace w{dim+1}_p0/w{dim+1}_p1 in the user parameter and rerun.')
    uc = ng.pipe_proc.make_uc(dic, data)
    tmp_xdata = uc.ppm_scale()
    xdata = np.linspace(tmp_xdata[0], tmp_xdata[-1], 
                    num=data.shape[len(data.shape)-1])
    uc2 = ng.pipe_proc.make_uc(dic, data, dim=len(data.shape)-2)
    tmp_ydata = uc2.ppm_scale()
    ydata = np.linspace(tmp_ydata[0], tmp_ydata[-1], 
                    num=data.shape[len(data.shape)-2])
    d = show_phaser2D(s)
    d.set_data(data, xdata, ydata, noise, index=0)
    print(f'!!! Update w{dim+1}_p0/w{dim+1}_p1 and rerun.')
    raise SystemExit  
  # Delete Imaginary
  dic, data = ng.pipe_proc.di(dic, data)
  # Extraction
  if ext_list[dim] == 'left':
    dic, data = ng.pipe_proc.ext(dic, data, left=True)
  elif ext_list[dim] == 'right':
    dic, data = ng.pipe_proc.ext(dic, data, right=True)
  # Baseline Correction
  if median_baseline_correction[dim]:
    data = ng.proc_bl.med(data)
    data = ng.proc_base.tp(data)
    data = ng.proc_bl.med(data)
  if polyfit_baseline_correction[dim]:
    try:
      data -= pokyphaser.fd_poly_bl(data, deg=2)
    except:
      print('Warning. Update your POKY to use fd_poly_bl.')
      continue
  # Reverse Data
  if reverse_list[dim]:
    dic, data = ng.pipe_proc.rev(dic, data, sw=True)
  # Transpose
  dic, data = ng.pipe_proc.tp(dic, data)

# Convert and Save
C.from_pipe(dic, data)

# create ucsf data and then write it out
output_file = s.save_filedialog('Save your new .ucsf file as...',
                                'UCSF NMR (*.ucsf);; Any (*)',
                                input_path)

ng.sparky.write(output_file,*C.to_sparky(),overwrite=True)

# Open in POKY
if LoadData:
  s.open_spectrum(output_file)

print('Done.')