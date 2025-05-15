#
# This is an example script for aligning two spectra.
#
# This uses image phase_cross_correlation. 
# Aligns current viewing X and Y dimensions. 
# So that viewing axis order must match.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: May 15, 2025
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
#

import numpy as np
import nmrglue as ng
from skimage.registration import phase_cross_correlation
from skimage.transform import resize

from sputil import name_to_view

# POKY libraries
import __main__
s = __main__.main_session

viewname = s.show_viewselectiondialog('Select a reference view', 0)
ref_v = name_to_view(viewname, s)

if ref_v == None:
  raise SystemError
ref_sp = ref_v.spectrum  

viewname = s.show_viewselectiondialog('Select a target view to align', 
                                          0)
aln_v = name_to_view(viewname, s)
aln_sp = aln_v.spectrum

if aln_v == None:
  raise SystemError

if ref_v == aln_v:
  s.show_message('Error', 'Identical spectrum cannot be aligned.')
  raise SystemError

ref_xy = ref_sp.nuclei[ref_v.axis_order[0]] + ref_sp.nuclei[ref_v.axis_order[1]]
aln_xy = aln_sp.nuclei[aln_v.axis_order[0]] + aln_sp.nuclei[aln_v.axis_order[1]]

if ref_xy != aln_xy:
  s.show_message('Error', 
    f'Nuclei do not match. {ref_xy} <> {aln_xy}. \nUse "xx" and "xy" to match first.')
  raise SystemError

def project_to_2d_max_abs_simple(array, projection_axes):
  ndim = array.ndim
  if ndim == 2:
    return np.transpose(array, axes=projection_axes)
  else:
    transposed_array = array.transpose(projection_axes)
    max_abs_array = np.max(np.abs(transposed_array), axis=ndim-1)
    if ndim == 3:
      return max_abs_array
    return np.max(max_abs_array, axis=2)

def preprocessing(data):
  data = np.abs(data)
  std = np.std(data)
  data = np.divide(data, np.sqrt(std))
  min_value = np.amin(data)
  max_value = np.amax(data)
  if max_value - min_value != 0:
    data = (data - min_value) / (max_value - min_value)
  return data

def shift_and_align_images_with_preprocessing(data, data2):
  data = preprocessing(data)
  data2 = preprocessing(data2)
  
  h1, w1 = h1_, w1_ = data.shape
  
  if h1 < 2048:
    h1 = 2048

  if w1 < 2048:
    w1 = 2048


  if (h1, w1) != (h1_, w1_):
    data = resize(data, (h1, w1), anti_aliasing=True)

  data2 = resize(data2, (h1, w1), anti_aliasing=True)
        
  # Perform subpixel registration
  shift, error, _ = phase_cross_correlation(data, data2, upsample_factor=1)
  shift = [shift[0] * h1_ / h1, shift[1] * w1_ / w1] 

  return shift, error

x_ref_ppm_per_pt = ref_sp.spectrum_width[ref_v.axis_order[0]] / \
  ref_sp.data_size[ref_v.axis_order[0]]

y_ref_ppm_per_pt = ref_sp.spectrum_width[ref_v.axis_order[1]] / \
  ref_sp.data_size[ref_v.axis_order[1]]

ref_xy = y_ref_ppm_per_pt / x_ref_ppm_per_pt

x_aln_ppm_per_pt = aln_sp.spectrum_width[aln_v.axis_order[0]] / \
  aln_sp.data_size[aln_v.axis_order[0]]

y_aln_ppm_per_pt = aln_sp.spectrum_width[aln_v.axis_order[1]] / \
  aln_sp.data_size[aln_v.axis_order[1]]

aln_xy = y_aln_ppm_per_pt / x_aln_ppm_per_pt


ref_dic, ref_data = ng.sparky.read(ref_sp.data_path)
ref_proj = project_to_2d_max_abs_simple(ref_data, ref_v.axis_order)

aln_dic, aln_data = ng.sparky.read(aln_sp.data_path)
aln_proj = project_to_2d_max_abs_simple(aln_data, aln_v.axis_order)

ratio = aln_xy / ref_xy

if ratio < 1.0:
  aln_proj = resize(aln_proj, [int(aln_proj.shape[0] / ratio + 0.5), 
                             aln_proj.shape[1]], anti_aliasing=True)
else:
  aln_proj = resize(aln_proj, [aln_proj.shape[0], 
                             int(aln_proj.shape[1] * ratio + 0.5)], 
                             anti_aliasing=True)
                     
shift, error = shift_and_align_images_with_preprocessing(ref_proj, aln_proj)

x_offset = shift[0] * x_ref_ppm_per_pt
y_offset = shift[1] * y_ref_ppm_per_pt

offset = list(aln_sp.scale_offset)
offset[aln_v.axis_order[0]] = ref_sp.scale_offset[ref_v.axis_order[0]] + x_offset
offset[aln_v.axis_order[1]] = ref_sp.scale_offset[ref_v.axis_order[1]] + y_offset

print(offset)
aln_sp.scale_offset = tuple(offset)
s.show_message('Finished.', f'Offset: {offset} applied.')