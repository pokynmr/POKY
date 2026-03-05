#
# This is an example script for generating HSQC0 of the 
# selected spectra in the project.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: March 5, 2026
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# POKY 01/14/22e or higher is required to run this script.
#

import numpy as np
import nmrglue as ng
from skimage.restoration import (denoise_wavelet, estimate_sigma)
from sputil import name_to_spectrum
import os.path
import __main__

noise_cut = 1.0
denoise = False

# POKY setup
s = __main__.main_session
specname = s.show_spectrumselectiondialog('Select all HSQC0 spectra', 1)
specname_list = specname.split('\t')
num_spec = len(specname_list)

if num_spec < 2:
    raise SystemError("At least 2 spectra are required for extrapolation.")

sp_list = list(map(lambda x: name_to_spectrum(x, s), specname_list))
new_path = s.save_filedialog('Save HSQC0 Spectrum', 'UCSF (*.ucsf)', 
                            os.path.dirname(sp_list[0].data_path))

# --- Data Loading & Sorting ---
data_list = []
for sp in sp_list:
    dic, data = ng.sparky.read(sp.data_path)
    # Optional Denoising could be added here loop-wise
    data_list.append((np.average(data), data, dic))

# Sort by intensity (Descending: y1 is highest/shortest delay, yn is lowest)
data_list.sort(key=lambda x: x[0], reverse=True)

# Stack data: Shape (n, rows, cols)
full_stack = np.array([item[1] for item in data_list])
base_dic = data_list[0][2]
y1_data = full_stack[0]

# --- Vectorized Linear Regression (Any N) ---
# x represents the increment index (1, 2, ..., n)
x = np.arange(1, num_spec + 1)
n = float(num_spec)

# Estimate noise from the first (strongest) spectrum
noise = np.std(y1_data[:50, :50])

# 1. Log-transform (Stabilizes math and prevents overflow)
log_stack = np.log(np.maximum(full_stack, 1e-5))

# 2. Linear Least Squares for Intercept (at x=0)
# Formula: intercept = (sum(y)*sum(x^2) - sum(x)*sum(xy)) / (n*sum(x^2) - sum(x)^2)
sum_x = np.sum(x)
sum_xx = np.sum(x**2)
sum_y = np.sum(log_stack, axis=0) # Sum across the spectrum stack
sum_xy = np.sum(x[:, None, None] * log_stack, axis=0)

denom = (n * sum_xx - sum_x**2)
intercept = (sum_y * sum_xx - sum_x * sum_xy) / denom

# 3. Re-exponentiate to find intensity at t=0
y_zero = np.exp(intercept)

# --- Logic Filtering ---
# Check if the trend is generally decaying (y1 > yn)
decay_mask = (y1_data > full_stack[-1]) & (y1_data > noise * noise_cut)

# Apply results: if decaying, use y_zero; else, keep original y1
result = np.where(decay_mask, y_zero, y1_data)

# Final Sanity Check: Prevent massive noise spikes (cap at 5x original)
result = np.minimum(result, y1_data * 5)

# --- Save ---
ng.sparky.write(new_path, base_dic, result.astype('float32'), overwrite=True)

if s.show_message_yes_no('Complete', f'Processed {num_spec} spectra. Load now?'):
    s.open_spectrum(new_path)
