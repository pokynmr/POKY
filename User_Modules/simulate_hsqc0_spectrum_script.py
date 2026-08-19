#
# This is an example script for generating HSQC0 of the
# selected spectra in the project.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: August 18, 2026
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# POKY 01/14/22e or higher is required to run this script.
#
# ---------------------------------------------------------------------------
# Changes in this revision (2026-08-18)
#
# 1. The final cap no longer amplifies negative points.
#
#    It was:   result = np.minimum(result, y1_data * 5)
#
#    That assumes every intensity is positive. A real 2D plane is roughly
#    half negative, and for a measured value of -1.5e8 the "cap" evaluates
#    to -7.5e8 -- so np.minimum does not bound that point, it multiplies it
#    by five. Because the line ran *after* np.where(), it also modified
#    points the decay mask had already declined to extrapolate, so the
#    "keep the measured value" branch did not actually keep it.
#
#    Now the *magnitude* is bounded, and only where the fit was used.
#
# 2. The logarithm floor scales with the data.
#
#    It was:   log_stack = np.log(np.maximum(full_stack, 1e-5))
#
#    1e-5 is an absolute constant. Real HSQC intensities here are ~1e12, so
#    it never engages; on normalized or scaled data it swallows the signal
#    entirely. It is now a small fraction of the measured noise, which keeps
#    the intent (keep log() finite) at any intensity scale.
#
# 3. scikit-image is imported only if denoising is switched on.
#
#    The import was unconditional while `denoise` was never used, so the
#    script failed at import for anyone without scikit-image installed.
# ---------------------------------------------------------------------------

import numpy as np
import nmrglue as ng
from sputil import name_to_spectrum
import os.path
import __main__

noise_cut = 1.0
denoise = False

# Ceiling on the extrapolation, as a multiple of the measured intensity.
# This is load-bearing, not cosmetic: extrapolating an exponential through
# three noisy points is violently unstable wherever there is no real signal.
# On a real three-spectrum series the uncapped fit reaches 5.4e22 against a
# true maximum of 1.2e12.
max_gain = 5.0

if denoise:
    from skimage.restoration import (denoise_wavelet, estimate_sigma)

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
full_stack = np.array([item[1] for item in data_list], dtype=float)
base_dic = data_list[0][2]
y1_data = full_stack[0]

# --- Vectorized Linear Regression (Any N) ---
# x represents the increment index (1, 2, ..., n)
x = np.arange(1, num_spec + 1)
n = float(num_spec)

# Estimate noise from the first (strongest) spectrum
noise = np.std(y1_data[:50, :50])

# 1. Log-transform (Stabilizes math and prevents overflow)
# Floor scaled to the noise rather than fixed at 1e-5 -- see note 2 above.
floor = max(abs(noise) * 1e-6, np.finfo(float).tiny)
log_stack = np.log(np.maximum(full_stack, floor))

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

# Final Sanity Check: bound the extrapolation by magnitude, and apply it
# only where the fit was actually used, so points that failed the decay or
# noise test come back exactly as measured -- see note 1 above.
limit = np.abs(y1_data) * max_gain
y_zero_capped = np.clip(y_zero, -limit, limit)

# Apply results: if decaying, use the capped extrapolation; else keep y1
result = np.where(decay_mask, y_zero_capped, y1_data)

# --- Save ---
ng.sparky.write(new_path, base_dic, result.astype('float32'), overwrite=True)

if s.show_message_yes_no('Complete', f'Processed {num_spec} spectra. Load now?'):
    s.open_spectrum(new_path)
