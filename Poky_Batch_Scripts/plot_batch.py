#
# This is an example script to open a spectrum and 
# save a PostScript file from Poky by a batch script.
# It will save in the same place where the data exists.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this batch script:
#   In the terminal,
#     $ poky -b plot_batch.py


# User parameters
# spectrum to open
#user_data_list = ['/path/to/my/spectrum.ucsf',
#                  '/path/to/my/spectrum2.ucsf',
#                  '/path/to/my/spectrum3.ucsf']
user_data_list = ['/Users/whlee/work/test/Sparky_2A/Spectra/wt2A_HSQC.ucsf', ]
quit_after_batch = True


import __main__
s = __main__.main_session

for user_data in user_data_list:
  if s.open_file(user_data) != 1:
    s.quit_poky() # close the program if opening a file failed.

proj = s.project

import nmrglue as ng
import numpy as np

from matplotlib import use as matplotlib_use
matplotlib_use('Agg')

import os
import matplotlib.pyplot as plt
from sputil import name_to_view

# process color
def process_color(tki, user_color, safe_color='red'):
  # color process
  try:
    rgb = tki.winfo_rgb(user_color)
    color = '#%02x%02x%02x' % (int(rgb[0]/256), 
                    int(rgb[1]/256), int(rgb[2]/256))
  except:
    color = safe_color
  return color

# Use NMRGlue to draw a spectrum in matplotlib
def postscript_spectrum(s, user_view_name, out_name):
  # , draw_peaks=False, draw_labels=False, 
  #                draw_lines=False, draw_grids=False
  v = name_to_view(user_view_name, s)
  if v == None:
    return False
  sp = v.spectrum
  se = v.session
  fig = plt.figure()
  offset = sp.scale_offset
  data_path = sp.data_path
  pos_levels = v.positive_levels
  neg_levels = v.negative_levels
  center = v.center
  # User can play with region to draw
  # default is full extent
  region = sp.region
  #tregion = list(v.region)
  #region = list(list(l1) for l1 in tregion)
  #sregion = sp.region
  #for i in range(len(center)):
  #  region[0][i] = max(region[0][i], sregion[0][i]) 
  #  region[1][i] = min(region[1][i], sregion[1][i])
  pos_cl = pos_levels.lowest * pos_levels.factor \
            ** np.arange(pos_levels.levels)
  pos_color = pos_levels.color
  temp_neg_cl = neg_levels.lowest * neg_levels.factor \
            ** np.arange(neg_levels.levels)
  neg_cl = np.flipud(temp_neg_cl)
  neg_color = neg_levels.color
  dic, data = ng.sparky.read_lowmem(data_path)
  
  # color process
  pos_color = process_color(se.tk, pos_color, 'red')
  neg_color = process_color(se.tk, neg_color, 'green')
  
  # for unit conversion
  #
  x_idx, y_idx = v.axis_order[0:2] 
  x_center = center[x_idx]
  x_width = region[1][x_idx] - region[0][x_idx]
  x_min, x_max = region[1][x_idx], region[0][x_idx]
  y_min, y_max = region[0][y_idx], region[1][y_idx]

  uc_x = ng.sparky.make_uc(dic, data, x_idx)
  uc_y = ng.sparky.make_uc(dic, data, y_idx)
  if len(v.axis_order) > 2:
    z_idx = v.axis_order[2]
    uc_z = ng.sparky.make_uc(dic, data, z_idx)
  if len(v.axis_order) > 3:
    a_idx = v.axis_order[3]
    uc_a = ng.sparky.make_uc(dic, data, a_idx)
    
  xmin_idx = uc_x(x_min - offset[x_idx], "ppm")
  xmax_idx = uc_x(x_max - offset[x_idx], "ppm")
  ymin_idx = uc_y(y_min - offset[y_idx], "ppm")
  ymax_idx = uc_y(y_max - offset[y_idx], "ppm")

  if len(v.axis_order) > 2:
    z_plane = center[z_idx]
    zc_idx = uc_z(z_plane - offset[z_idx], "ppm")
    
  if len(v.axis_order) > 3:
    a_plane = center[a_idx]
    ac_idx = uc_a(a_plane - offset[a_idx], "ppm")

  if ymin_idx > ymax_idx:
    ymin_idx, ymax_idx = ymax_idx, ymin_idx

  # extract data
  if len(v.axis_order) == 2:
    if y_idx == 1:
      strip = data[xmin_idx:xmax_idx+1, ymin_idx:ymax_idx+1]
    else:
      strip = data[ymin_idx:ymax_idx+1, xmin_idx:xmax_idx+1]
  elif len(v.axis_order) == 3:  
    if [x_idx, y_idx, z_idx] == [0, 1, 2]:
      strip = data[xmin_idx:xmax_idx+1, ymin_idx:ymax_idx+1, zc_idx]
    elif [x_idx, z_idx, y_idx] == [0, 1, 2]:
      strip = data[xmin_idx:xmax_idx+1, zc_idx, ymin_idx:ymax_idx+1]
    elif [y_idx, x_idx, z_idx] == [0, 1, 2]:
      strip = data[ymin_idx:ymax_idx+1, xmin_idx:xmax_idx+1, zc_idx]
    elif [y_idx, z_idx, x_idx] == [0, 1, 2]:
      strip = data[ymin_idx:ymax_idx+1, zc_idx, xmin_idx:xmax_idx+1]
    elif [z_idx, x_idx, y_idx] == [0, 1, 2]:
      strip = data[zc_idx, xmin_idx:xmax_idx+1, ymin_idx:ymax_idx+1]
    elif [z_idx, y_idx, x_idx] == [0, 1, 2]:
      strip = data[zc_idx, ymin_idx:ymax_idx+1, xmin_idx:xmax_idx+1]
  elif len(v.axis_order == 4):
    from itertools import permutations
    dict1 = {0: list(range(xmin_idx, xmax_idx+1)), 
             1: list(range(ymin_idx, ymax_idx+1)),
             2: zc_idx,
             3: ac_idx}
    perm = permutations([0, 1, 2, 3])
    for i in list(perm):
      rl = list[i]
      if reflist == [x_idx, y_idx, z_idx, a_idx]:
        strip = data[dict1[rl[0]], dict1[rl[1]], dict1[rl[2]], dict1[rl[3]]]
        break

  # determine ppm limits of contour plot
  uc_xs = uc_x.ppm_scale()
  strip_ppm_x = uc_xs[xmin_idx:xmax_idx+1]
  
  uc_ys = uc_y.ppm_scale()
  strip_ppm_y = uc_ys[ymin_idx:ymax_idx+1]
  strip_x, strip_y = np.meshgrid(strip_ppm_x + offset[x_idx], 
                                strip_ppm_y + offset[y_idx])
  
  # add contour plot of strip to figure
  # positive
  ax = fig.add_subplot(1, 1, 1)
  if y_idx > x_idx:
    strip = strip.transpose()

  try:  
    ax.contour(strip_x, strip_y, strip, pos_cl, 
          colors=pos_color, linewidths=0.5, linestyles='solid')
  except:
    pass

  try:
    ax.contour(strip_x, strip_y, strip, neg_cl, 
          colors=neg_color, linewidths=0.5, linestyles='solid')
  except:
    pass
  ax.set_title(v.name)
  ax.invert_yaxis()  # flip axes since ppm indexed
  ax.invert_xaxis()
  # label and put ticks on first strip plot
  ax.tick_params(axis='both', labelbottom=True, bottom=True, top=False, 
      labelleft=True, left=True, right=False, direction='out')
  #ax.set_xlabel('%.2f'%(x_center)) #, size=6)

  nuc = sp.nuclei[y_idx]
  nuc2 = sp.nuclei[x_idx]
  ax.set_ylabel("%s (ppm)" % (nuc))
  ax.set_xlabel("%s (ppm)" % (nuc2))

  if v.show_ornaments:
    # draw peaks
    if v.show_peaks:
      for peak in sp.peak_list():
        # convert to indices
        # z, a matches
        # plot
        #ax.plot([x0, x1], [y0, y1], 'k')
        pass
  
    # draw labels
    if v.show_labels:
      for label in sp.label_list():
        # convert to indices
        # z, a matches
        # plot
        #ax.text(x1 + 1, y0, level.text, 'k')
        pass

    # draw lines
    if v.show_lines:
      for line in sp.line_list():
        pass

    # draw grids
    if v.show_grids:
      for grid in sp.grid_list():
        pass
  plt.savefig(out_name)
  return True    


for user_data in user_data_list:
  view_name = os.path.splitext(os.path.basename(user_data))[0]
  out_name = os.path.splitext(user_data)[0] + '.ps'
  postscript_spectrum(s, view_name, out_name)

if quit_after_batch:
  s.quit_poky()

# If you wish to let Poky opened and use after, uncomment below.
# s.set_verbose(True)