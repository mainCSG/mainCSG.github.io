'''
Exchange calculation for a potential
'''

import numpy as np
from scipy import special
from scipy.constants import hbar, e

import matplotlib.pyplot as plt

# adding the ...QuDiPy to path
from pathlib import Path
import sys

p = Path(__file__).parents[2]
sys.path.append(str(p))

from qudipy.utils import Constants
import qudipy.potential as pot


def hl_quartic(gparams, #  y_slice=0, 
                    material='vacuum'):
    '''
    Function calculates exchange values in Heitler-London approximation
    for an array of dots. Each pair is fitted to a quartic function for the
    calculation.
    y-slice is determined based on the positions of the dot minima.
    Formula borrowed from Burkard & Loss, PRB 89 (3), 1999, pgs. 2073-2074  

    Parameters:
    -------------
    
    gparams: GridParameters Object
        Must contain 2D potential, and x and y coordinates as a meshgrid object


    Keyword Arguments:
    -------------------
    material: String, optional
        The string defining the material that the dots are
        in. Default is 'vacuum'.

    Returns:
    --------------
    J: float or tuple of floats
        Exchange value(s) for the array of dots
    '''
    
    #material constants
    consts = Constants(material)
    m = consts.me
    eps = consts.eps    # epsilon * epsilon_0

    # obtaining parameters from fitting to a quartic potential
    fit_params = pot.manipulate.fit_quartic(gparams, material=material, 
                                                            return_params=True)
    a = fit_params['dot_sep'] / 2
    e_field_x = fit_params['e_field_x']
    omega_0 = fit_params['omega_0']

    # Bohr radius
    a_B = np.sqrt(hbar/ (m * omega_0))
    # dimensionless distance
    d = a/ a_B
    # interaction to confinement energy ratio
    c = np.sqrt(np.pi/ 2) * ( e**2 / (4 * np.pi * eps * a_B)) / (hbar * omega_0)
      
    J_hl = ((hbar * omega_0 / np.sinh(2 * d**2)
                * (c * (np.exp(-d**2) * special.iv(0, d**2) -1 ) 
                    + 3 / 4 * (1 + d**2))) 
                        + (1 / np.sinh(2 * d**2)) * (3 / (2 * d**2)) *
                            (e * e_field_x * a )**2 / (hbar * omega_0)
                )
    return J_hl
    

def hm_quartic(gparams, material='vacuum'):
    '''
    Function calculates exchange values in Hund-Mulliken approximation
    for an array of dots. Each pair is fitted to a quartic function for the
    calculation.
    y-slice is determined based on the positions of the dot minima.
    Formulas borrowed from Burkard & Loss, PRB 89 (3), 1999, pgs. 2073-2074  

    Parameters:
    -------------
    
    gparams: GridParameters Object
        Must contain 2D potential, and x and y coordinates as a meshgrid object

    Keyword Arguments:
    -------------------

    material: String, optional
        The string defining the material that the dots are
        in. Default is 'vacuum'.

    Returns:
    --------------
    J: float or tuple of floats
        Exchange value(s) for the array of dots
    '''
        
    # material constants
    consts = Constants(material)
    m = consts.me
    eps = consts.eps    # epsilon * epsilon_0

    # obtaining parameters from fitting to a quartic potential
    fit_params = pot.manipulate.fit_quartic(gparams, material=material, 
                                                        return_params=True)
           
    a = fit_params['dot_sep'] / 2
    e_field_x = fit_params['e_field_x']
    omega_0 = fit_params['omega_0']

    # Bohr radius
    a_B = np.sqrt(hbar/ (m * omega_0))
    # dimensionless distance
    d = a / a_B
    # interaction to confinement energy ratio
    c = np.sqrt(np.pi/ 2) * ( e**2 / (4 * np.pi * eps * a_B)) / (hbar * omega_0)

    # No out-of plane field
    b = 1

    # Hund-Mulliken related parameters
    S = np.exp( d**2 * (1 / b - 2 * b))
    g = (1 - np.sqrt(1 - S ** 2)) / S
    N = 1/ np.sqrt(1 - 2 * S * g + g **2)
    t = 3 / 8 * S / (1 - S**2) * (1 / b + d**2)
    
    #-------------------#

    F_1 = c * np.sqrt(b)
    F_2 = c * np.sqrt(b) * np.exp(- b * d**2) * special.iv(0, b * d**2)
    F_3 = c * np.sqrt(b) * np.exp(d**2 * (b - 1 / b)) * special.iv(0, d**2 * (b - 1 / b))
    F_4 = (c * np.sqrt(b) * np.exp(- d**2 / (4 * b)) *
                           sum([ (-1)**k * special.iv(2 * k, (d**2 * ( 2 * b - 1 / b) / 4)) * 
                                special.iv(2*k, 1j * d**2 / 2 * np.sqrt(b**2 - 1))
                                for k in (-1,0,1)]))
    
    V_plus = N**4 * (4 * g**2 * (1 + S**2) * F_1 + (1 + g**2) **2 * F_2 
                     + 4 * g**2 * F_3 - 16 * g**2 * F_4)

    V_minus = N**4 * (1 - g**2)**2 * (F_2 - S**2 * F_3)

    U = N**4 * ((1 + g**4 + 2 * g**2 * S**2) * F_1 + 2 * g**2 * F_2 + 
                2 * g**2 * S**2 * F_3 - 8 * g **2 * F_4)

    X = N**4 * (((1 + g**4) * S**2+ 2 * g**2 ) * F_1 + 2 * g**2 * F_2 + 
                2 * g**2 * S**2 * F_3 - 8 * g **2 * F_4)
     
    w = N**4 * (-g * (1 + g**2) *(1 + S**2) * F_1 - g * (1 + g**2) * F_2 
               -g * (1 + g**2) * S**2 * F_3 + (1 + 6 * g**2 + g**4) * S * F_4)

    t_H = t - w
    V = V_minus - V_plus
    U_H = U - V_plus + X
    
    # Final expression
    J_hm = hbar * omega_0 *(V - U_H / 2 + 1 / 2 * np.sqrt(U_H**2 + 16 * t_H**2)   +    
                ((1 / np.sinh(2 * d**2)) * (3 / (2 * d**2)) *
                            (e * e_field_x * a )**2 / (hbar * omega_0)**2))

    return np.real(J_hm)