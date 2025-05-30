#
# This is an example script to plot CSPs.
#
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
import __main__
s = __main__.main_session

print('\n\n\n------------------------------------------------------')
print('POKY Plot CSPs')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

weight_dict = {'1H': 1.0, '15N': 6.0, '13C': 8.0}

spec_names = s.show_spectrumselectiondialog('Select twp spectra', 1)

try:
  spec_name1, spec_name2 = spec_names.split('\t')
except:
  s.show_message('Error', 'Select two spectra.')
  raise SystemError

from numpy import std
from sputil import name_to_spectrum, sort_peaks_by_assignment
import matplotlib.pyplot as plt
from matplotlib import use as matplotlib_use
matplotlib_use('TkAgg')

spec1 = name_to_spectrum(spec_name1, s)
spec2 = name_to_spectrum(spec_name2, s)

if spec1.dimension != 2 or spec2.dimension != 2:
  s.show_message('Error', 'Select 2D NMR spectra.')
  raise SystemError
  

peaks1 = spec1.peak_list()
peaks2 = spec2.peak_list()

peaks1 = list(filter(lambda x: x.is_assigned, peaks1))
peaks2 = list(filter(lambda x: x.is_assigned, peaks2))

peaks1 = sort_peaks_by_assignment(peaks1, None)
peaks2 = sort_peaks_by_assignment(peaks2, None)

B1 = weight_dict[spec1.nuclei[0]], weight_dict[spec1.nuclei[1]]
B2 = weight_dict[spec2.nuclei[0]], weight_dict[spec2.nuclei[1]]

n1 = list(map(lambda x: x.resonances()[0].group.number, peaks1))
n2 = list(map(lambda x: x.resonances()[0].group.number, peaks2))
cs1 = list(map(lambda x: (x.frequency[0]/B1[0], x.frequency[1]/B1[1]), 
               peaks1))
cs2 = list(map(lambda x: (x.frequency[0]/B2[0], x.frequency[1]/B2[1]), 
               peaks2))
cs2_dict = dict(zip(n2, cs2))

xdata, ydata = [], []

print('\nResidue\tCSP(ppm)\n')
for i, nseq in enumerate(n1):
  d11, d12 = cs1[i]
  try:
    d21, d22 = cs2_dict[nseq]
  except:
    continue
  y = min(((d11-d21)**2 + (d12-d22)**2)**.5, ((d11-d22)**2 + (d12-d21)**2)**.5)
  xdata.append(nseq)
  ydata.append(y)
  print(f'{nseq}\t{y}')


def plot_poky_csp(residue_numbers, csps):
    if len(residue_numbers) != len(csps):
        print("Error: X and Y data lists must have the same length.")
        return

    # Calculate standard deviation
    std_dev = std(csps)

    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(residue_numbers, csps, color='blue', alpha=0.7, edgecolors='w', s=50)

    # Add horizontal lines at 1 and 2 standard deviations
    plt.axhline(y=std_dev, color='red', linestyle='--', label='1σ')
    plt.axhline(y=2 * std_dev, color='green', linestyle='--', label='2σ')

    # Annotate the lines
    plt.text(residue_numbers[-1], std_dev, '1σ', color='red', verticalalignment='bottom', horizontalalignment='right')
    plt.text(residue_numbers[-1], 2 * std_dev, '2σ', color='green', verticalalignment='bottom', horizontalalignment='right')

    # Set plot title and axis labels
    plt.title("POKY CSP PLOT", fontsize=16)
    plt.xlabel("Residue Number", fontsize=12)
    plt.ylabel("CSP (ppm)", fontsize=12)

    # Add a grid
    plt.grid(True, linestyle='--', alpha=0.6)

    # Customize ticks
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # Add legend
    plt.legend()

    # Add tight layout
    plt.tight_layout()

    # Display the plot
    plt.show(block=False)

plot_poky_csp(xdata, ydata)