#
# This is an example script for performing TRACT minimization analysis.
#   by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
#   Original code by Scott A. Robson, Ph.D.
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#     and select a directory that has fid file in.
#
# Reference:
#   Scott A. Robson, Çağdaş Dağ, Hongwei Wu & Joshua J. Ziarek
#     TRACT revisited: an algebraic solution for determining overall 
#     rotational correlation times from cross-correlated relaxation rates. 
#     J Biomol NMR 75, 293–302 (2021). 
#     https://doi.org/10.1007/s10858-021-00379-5
#
#   https://github.com/nomadiq/TRACT_analysis
#

import numpy as np
from scipy.optimize import minimize_scalar, minimize
import __main__
s = __main__.session

# field in MHz
try:
    field = float(s.show_inputdialog('Field in MHz', 
                                     'Type magnetic strength (MHz)', '750'))    
except:
    s.show_message('Error', 'Value must be numberic.')
    raise SystemError

# measured relaxation rates
try:
    Rb = float(s.show_inputdialog('Relaxation Rates', 
                                     'Measured relaxation rate Rb (Hz) ?', '64'))    
except:
    s.show_message('Error', 'Value must be numberic.')
    raise SystemError

try:
    Ra = float(s.show_inputdialog('Relaxation Rates', 
                                     'Measured relaxation rate Ra (Hz) ?', '13'))    
except:
    s.show_message('Error', 'Value must be numberic.')
    raise SystemError

# constants
h = 6.62607004 * (1/np.power(10,34, dtype=np.longdouble))    # Plank's
mu_0 = 1.25663706 * (1/np.power(10,6, dtype=np.longdouble))  # vacuum permeability
gamma_H = 267.52218744 * np.power(10,6, dtype=np.longdouble) # proton gyromagnetic ratio
gamma_N = -27.116 * np.power(10,6, dtype=np.longdouble)      # 15N gyromagnetic ratio
r = 1.02 * (1/np.power(10,10, dtype=np.longdouble))          # internuclear distance
delta_dN = 160 * (1/np.power(10, 6))                         # diff in axially symetric 15N CS tensor
theta = 17*np.pi/180                                         # angle between CSA axis and N-H bond

# derived field value in Tesla
B_0 = field * np.power(10, 6, dtype=np.longdouble) * 2 * np.pi / gamma_H # in Tesla

# equation (5)
p = mu_0*gamma_H*gamma_N*h/(16*np.pi*np.pi*np.sqrt(2)*np.power(r,3))     # DD 1H-15N bond

# equation (6)
dN = gamma_N*B_0*delta_dN/(3*np.sqrt(2)) # 15N CSA

# equation (7)
w_N = B_0 * gamma_N                      # 15N frequency (radians/s)

args = (w_N, Rb, Ra, p, dN, theta)

# Equation 12 function
def objective_function(tau_c, *args):
    
    #spectral density function
    def J(w_N, tau_c):
        return 0.4*tau_c/(1+(w_N**2*tau_c**2))
    #print(*args)
    (w_N, Rb, Ra, p, dN, theta) = args # unpact these constants inside the function
    return np.abs((4*J(0, tau_c) + 3*J(w_N, tau_c)) - ((Rb - Ra)/(2*p*dN*(3*np.cos(theta)**2-1))))

# guess tau_c 
t = 10 * (1/np.power(10,8, dtype=np.longdouble)) # guess 10 ns

msg = '-- Brent minimization --\n'
res = minimize_scalar(objective_function, args=args, method='Brent')
msg += f'Given: field = {field} MHz, Rb = {Rb} Hz and Ra = {Ra} Hz\n'
msg += f'tau_c: {res.x} seconds -- Using Brent Minimization\n'

msg += '\n-- BFGS minimization --\n'
res = minimize(objective_function, t, args=args, method='BFGS')
msg += f'Given: field = {field} MHz, Rb = {Rb} Hz and Ra = {Ra} Hz\n'
msg += f'tau_c: {res.x} seconds -- Using BFGS Minimization\n'

msg += '\n-- Powell minimization --\n'
res = minimize(objective_function, t, args=args, method='Powell')
msg += f'Given: field = {field} MHz, Rb = {Rb} Hz and Ra = {Ra} Hz\n'
msg += f'tau_c: {res.x} seconds -- Using Powell Minimization\n'

msg += '\n-- TNC minimization --\n'
res = minimize(objective_function, t, args=args, method='TNC')
msg += f'Given: field = {field} MHz, Rb = {Rb} Hz and Ra = {Ra} Hz\n'
msg += f'tau_c: {res.x} seconds -- Using TNC Minimization\n'

print(msg)
s.show_message('Results', msg)