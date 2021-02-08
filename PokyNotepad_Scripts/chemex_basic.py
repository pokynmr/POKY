#
# This is an example script to conduct ChemEx analysis.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

#
# ChemEx Reference:
#   Hsp70 biases the folding pathways of client proteins
#   Ashok Sekhar, Rina Rosenzweig, Guillaume Bouvignies, and Lewis E. Kay
#  PNAS May 17, 2016 113 (20) E2794-E2801; https://doi.org/10.1073/pnas.1601846113
#
#
# TOML keys and their meanings
# [experiment]
#   name          = experiment name
#   carrier       = RF carrier of studied nuclei in ppm
# [conditions]
#   h_larmor_freq	= magnetic field strength in MHz
#   label         =	labeling scheme of the sample
# [data]
#   path          = directory containing data files
#   error         = directory containing data files
#

# example toml_input
toml_input = """
[experiment]
name         = "dcest_15n"
time_t1      = 0.2
carrier      = 118.663
b1_frq       = 20.0
b1_inh_res   = 5
sw           = 800.0
pw90         = 36.2e-6

[conditions]
h_larmor_frq = 1000.3
# sample = "G48A Fyn SH3"
# temperature = 25.0

[data]
path = "../Data/800Hz/"
error = "file"
filter_offsets = [[0.0, 20.0]]
profiles = [
    ["F4N", "F4N-HN.out"],
    ["E5N", "E5N-HN.out"],
    ["A6N", "A6N-HN.out"],
    ["L7N", "L7N-HN.out"],
    ["Y8N", "Y8N-HN.out"],
]
"""


# Data files
# The location of data files is specified in experiment files.
# Data files typically contain three columns with the following information:
#
# Experiment            First column    Second column   Third column
# ------------
# CPMG                  ncyc            Intensity       Uncertainty
# CEST/DCEST/COSCEST    Offset (Hz)     Intensity       Uncertainty
# Relaxation            Time (s)        Intensity       Uncertainty
#

data_file = """
#Offset (Hz)        Intensity    Uncertainty
  -2.000e+04    1.6713620e+07  7.5500000e+04
  -4.000e+02    9.1779790e+06  7.5500000e+04
  -3.500e+02    5.3376120e+06  7.5500000e+04
  -3.000e+02    8.3675620e+05  7.5500000e+04
  -2.500e+02    8.5308660e+06  7.5500000e+04
  -2.000e+02    1.0703890e+07  7.5500000e+04
  -1.500e+02    1.1374080e+07  7.5500000e+04
  -1.000e+02    1.1372830e+07  7.5500000e+04
  -5.000e+01    1.1546840e+07  7.5500000e+04
   0.000e+00    1.1448260e+07  7.5500000e+04
   5.000e+01    1.1345880e+07  7.5500000e+04
   1.000e+02    1.1018610e+07  7.5500000e+04
   1.500e+02    8.8288670e+06  7.5500000e+04
   2.000e+02    5.9639350e+06  7.5500000e+04
   2.500e+02    9.3583360e+06  7.5500000e+04
   3.000e+02    1.0234890e+07  7.5500000e+04
   3.500e+02    1.0114990e+07  7.5500000e+04
   4.000e+02    9.1075310e+06  7.5500000e+04
"""

# Parameter files
# The parameter files (indicated with -p) contain initial estimate of parameters
# to be used during the fitting process, which typically looks like this:
#
parameter_file = """
[GLOBAL]
PB     =    0.6
KEX_AB =  130.0
TAUC_A =    4.0
"""

# If certain parameter is required but not included in the parameter files,
# a default value will be used to initialize, the initial value depends on
# each specific module. Due to the multidimensional feature of the
# minimization process, it is essential to set suitable initial parameters
# to avoid being trapped in a local minimum.


# Kinetic models
# The kinetic model (indicated with -d) indicates the type of exchange model
# to be used for the data analysis, which can be one of the following:

#  2st          2-state exchange model (default)
#  3st          3-state exchange model
#  4st          4-state exchange model
#  2st_rs       2-state residue-specific exchange model
#  2st_hd       2-state exchange model for H/D solvent exchange study
#  2st_eyring   2-state exchange model for temperature-dependent study
#  3st_eyring   3-state exchange model for temperature-dependent study
#  2st_binding	2-state exchange model for binding study
#  4st_hd       4-state exchange model for simutaneous normal and H/D solvent exchange study

kinetic_model = "2st"

import chemex.chemex
print(chemex.LOGO)

chemex.main()
