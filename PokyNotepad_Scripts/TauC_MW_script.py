#
# This is an example script to approximate M.W. from TauC.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module


print('\n\n\n------------------------------------------------------')
print('POKY M.W. approximation using TauC')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# User parameters
TauC = 6.0 # unit: ns
plotMW = True # Plot on top of NESG data

# http://www.nmr2.buffalo.edu/nesg.wiki/NMR_determined_Rotational_correlation_time
# NESG's empirical correlation
# between TauC (ns) and M.W. (kDa) for rigid proteins.
#
# TauC = 0.58 * MW + 0.472

# Tc (ns) = 1 / (4*pi*Nresonance_freq) * (6*T1/T2 - 7)**.5
# This means, when TauC is determined, MW can be predicted.

MW = (TauC - 0.472) / 0.58
print('TauC (ns): %f' % (TauC))
print('Approx. M.W. (kDa): %f' % (MW))


# Plot MW on NESG data
if plotMW:
  NESG_data = '''PsR76A  (NC5)	7.2	478.0	128.0	5.10
  VfR117 (NC)	11.2	605.0	119.0	6.30
  SyR11 (NC5)	12.4	630.0	104.0	7.10
  ER541-37-162 (NC5)	15.8	729.0	66.5	10.0
  ER540 (NC5)	18.8	909.0	66.5	11.3
  SoR190 (NC)	13.8	697.5	100.9	7.70
  TR80 (NC5)	10.5	612.8	102.9	7.00
  Ubiquitin (NC)	9.0	441.8	144.6	4.40
  HR2873B (NC)	10.7	492.0	115.0	5.70
  B-domain (NC)	7.2	423.5	153.3	4.05
  BcR97A (NC)	13.1	705.8	80.6	8.80
  PfR193A (NC)	13.6	733.9	80.9	9.00
  MvR76 (NC)	20.2	1015.0	64.5	12.2
  DvR115G (NC)	10.9	608.7	115.6	6.50
  MrR110B (NC5)	11.8	707.0	99.2	7.80
  VpR247 (NC5)	12.5	661.2	88.3	8.05
  BcR147A (NC)	11.9	645.0	104.0	7.20
  WR73 (NC5)	21.9	1261.0 41.3 13.0
  NsR431C (NC5)	16.8	855.5	71.2	10.6
  StR82 (NC)	9.2	537.3	100.4	6.6'''

  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.optimize import curve_fit
  from matplotlib import use as matplotlib_use
  matplotlib_use('TkAgg')

  ln = NESG_data.split('\n')
  MW_list = list(map(lambda x: float(x.strip().split()[2]), ln))
  Tc_list = list(map(lambda x: float(x.strip().split()[5]), ln))
  xdata = np.array(sorted(MW_list))
  ydata = np.array(list(map(lambda i: Tc_list[MW_list.index(xdata[i])],
                      range(len(xdata)))))

  def func(x, a, b):
    return a * x + b
  popt, pcov = curve_fit(func, xdata, ydata)
  #point_sd = list(np.sqrt(np.diag(pcov)))

  plt.figure()
  plt.plot(xdata, ydata, 'bo', label='NESG data')
  plt.plot(xdata, func(xdata, *popt), 'k-',
            label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
  plt.plot([MW,], [TauC], 'ro',
            label='%.3f kDa' % (MW), markersize=7)
  plt.xlabel('MW (kDa)')
  plt.ylabel('TauC (ns)')
  plt.legend()
  plt.pause(0.1)
  plt.show(block=False)
