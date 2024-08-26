#
# This is an example script to save postscripts of selected peaks
# view region will be from selected view.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#

import __main__
s = __main__.main_session

import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np
import os

fileext = 'pdf'

# we will use region width/height from selected view.
# we only support 2D for now.
v = s.selected_view()
if v == None:
  s.show_message('Error', 'Please select a view.')
  raise SystemError
if v.spectrum.dimension != 2:
  s.show_message('Error', 'Please select a 2D spectrum.')
  raise SystemError
if len(s.selected_peaks()) == 0:
  s.show_message('Error', 'Please select peaks.')
  raise SystemError

for peak in s.selected_peaks():
  if peak.assignment == '' or peak.assignment == '?-?-?':
    s.show_message('Error', 'Unassigned peak is selected.')
    raise SystemError
  
save_dir = s.open_directorydialog('Select a directory to save', '')

sw = v.region[1][0] - v.region[0][0]
sh = v.region[1][1] - v.region[0][1]

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

plot_view_list = [v,]
for ov in v.spectrum.session.project.overlay_list():
  if ov.overlay_onto == v:
    plot_view_list.append(ov.overlay_from)

x_idx, y_idx = v.axis_order[0:2]

spec_attr = {}
for i in range(len(plot_view_list)):
  pv = plot_view_list[i]
  sp = pv.spectrum
  offset = sp.scale_offset
  data_path = sp.data_path
  pos_levels = pv.positive_levels
  neg_levels = pv.negative_levels
  pos_cl = pos_levels.lowest * pos_levels.factor \
          ** np.arange(pos_levels.levels)
  pos_color = process_color(s.tk, pos_levels.color, 'red')
  temp_neg_cl = neg_levels.lowest * neg_levels.factor \
            ** np.arange(neg_levels.levels)
  neg_cl = np.flipud(temp_neg_cl)
  neg_color = process_color(s.tk, neg_levels.color, 'green')
  dic, data = ng.sparky.read(data_path)
  uc_x = ng.sparky.make_uc(dic, data, x_idx)
  uc_y = ng.sparky.make_uc(dic, data, y_idx)
  uc_xs = uc_x.ppm_scale()
  uc_ys = uc_y.ppm_scale()
  
  spec_attr[i] = pv, sp, offset, data_path, pos_levels, neg_levels, \
      pos_cl, pos_color, neg_color, dic, data, uc_x, uc_y, uc_xs, uc_ys

for i in range(len(s.selected_peaks())):
  peak = s.selected_peaks()[i]
  fig = plt.figure(peak.assignment)
  fig.clf()
  ax = fig.add_subplot(1, 1, 1)
  ax.clear()
  center = peak.frequency
  region = ((center[0] - sw/2, center[1] - sh/2), 
            (center[0] + sw/2, center[1] + sh/2))
  
  for j in range(len(plot_view_list)):
    pv, sp, offset, data_path, pos_levels, neg_levels, pos_cl, pos_color, \
      neg_color, dic, data, uc_x, uc_y, uc_xs, uc_ys = spec_attr[j]
    
    xmin_idx = uc_x(region[1][x_idx] - offset[x_idx], "ppm")
    xmax_idx = uc_x(region[0][x_idx] - offset[x_idx], "ppm")
    ymin_idx = uc_y(region[1][y_idx] - offset[y_idx], "ppm")
    ymax_idx = uc_y(region[0][y_idx] - offset[y_idx], "ppm")

    if y_idx == 1:
      strip = data[xmin_idx:xmax_idx+1, ymin_idx:ymax_idx+1]
    else:
      strip = data[ymin_idx:ymax_idx+1, xmin_idx:xmax_idx+1]

    # determine ppm limits of contour plot
    strip_ppm_x = uc_xs[xmin_idx:xmax_idx+1]
    strip_ppm_y = uc_ys[ymin_idx:ymax_idx+1]
    strip_x, strip_y = np.meshgrid(strip_ppm_x + offset[x_idx], 
                                  strip_ppm_y + offset[y_idx])
  
    # add contour plot of strip to figure
    # positive
    if y_idx > x_idx:
      strip = strip.transpose()
      strip_x, strip_y = strip_y, strip_x
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
  
  ax.set_title(peak.assignment)
  ax.invert_yaxis()  # flip axes since ppm indexed
  ax.invert_xaxis()

  # label and put ticks on first strip plot
  ax.tick_params(axis='both', labelbottom=True, bottom=True, top=False, 
      labelleft=True, left=True, right=False, direction='out')

  nuc = sp.nuclei[y_idx]
  nuc2 = sp.nuclei[x_idx]
  ax.set_ylabel("%s (ppm)" % (nuc))
  ax.set_xlabel("%s (ppm)" % (nuc2))

  from tkutil import get_view_geometry, get_dpi
  dpi = get_dpi(s, pv)
  w, h, x, y = get_view_geometry(s, pv)
  wi = w / dpi
  hi = h / dpi
  fig.set_figwidth(wi)
  fig.set_figheight(hi)
  filename = os.path.join(save_dir, peak.assignment + f'.{fileext}')
  plt.savefig(filename)

s.show_message('Done', 'Saving completed.')