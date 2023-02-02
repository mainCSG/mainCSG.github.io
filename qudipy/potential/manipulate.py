'''
Functions for the manipulation of potential data

@authors: Madi Schuetze, Bohdan Khromets
'''

from scipy.signal import find_peaks
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import numpy as np
from qudipy.utils.constants import Constants

def fit_quartic(gparams, material='vacuum',
                                        return_params=True):
    '''
    Returns the potential data of best fit to a quartic potential. 
    The data will be a two dot system  fitted 
    to a quartic well based on the positions of local maxima/minima and the
    potential energy values there.

    Note: RMS fitting of the potential landscape to a quartic function potential
    was used previously but was deprecated as it introduces very large fitting 
    error. Code can be found in qudipy/archive/rms_fit_quartic.py

    #TODO #64 update to the case of many dots

    Parameters:
    -------------

    gparams: GridParameters Object
        Must contain 2D potential, and x and y coordinates as a meshgrid object

    Keyword Arguments:
    -------------------
    
    material: String, optional
        The string defining the material that the dots are
        in. Default is 'vacuum'.

    return_params: Bool
        The option to return all of the fitting parameters
        as a dictionary with keys:
        * U_fit
        * dot_sep
        * e_field_x
        * omega_0
        * x_centre
        * U_0
        If False, only the potential landscape U_fit of best fit is returned

    Returns:
    --------------
    U_fit: 2D Arrayof(Float) 
    *or* result_dict: dictionary
        Contains Float values that map out the 1-D potential energy quartic fit
        of a two dot system.
    '''
    # Get x and y parameters
    x = gparams.x
    x_mesh = gparams.x_mesh

    y = gparams.y
    y_mesh = gparams.y_mesh
    
    U_data = gparams.potential

    # ** Find the maximum and minima **
    # find absolute minimum

    min_idx = np.unravel_index(U_data.argmin(), U_data.shape)
    #draw a y_0 through the minimum
    y_0 = y[min_idx[0]]     # initial guess

    # Find U_data_values of y_0 used
    try:
        abs_val_array = np.abs(y - y_0)
    except TypeError:
        print('y_0 must be the constant y-value at which the\
         wells are lined up on. Currently, y_0 = {}'.format(y_0))
    index = abs_val_array.argmin()
    U_data_slice =  U_data[index]

    # Maxima
    maxima = np.sort(find_peaks(U_data_slice)[0])
    # Minima
    minima = np.sort(find_peaks(-1 * U_data_slice)[0])


    # ** Determine Inital Guesses of Values **
    # Constants
    consts = Constants(material)
    e = consts.e # electron charge
    m = consts.me # mass of free electron
    
    # If the number of minima is less than 2, and/or 
    # number of maxima is less than 1, exchange cannot be calculated

    if (len(minima) != 2 or len(maxima) != 1  or x[minima[0]] >= x[maxima[0]]
                                               or x[maxima[0]]>= x[minima[-1]]):
        x_max = x_1 = x_2 = np.nan
        E_max = E1 = E2 = np.nan
    else:
        # Indicies
        i0 = maxima[0]
        i1 = minima[0]
        i2 = minima[-1]
        
        # discrete positions of extrema
        x_max_discrete = x[i0] # x-value of local maximum (initial guess)
        x_1_discrete = x[i1] # x-value of first minimum
        x_2_discrete = x[i2] # x-value of second minimum

        # find precise value of positions by interpolation + minimization
        # the potential 
        x_max, E_max = _precise_min(-1* U_data_slice, x, x_max_discrete, 
                        [x_max_discrete - 0.3* (x_max_discrete -x_1_discrete), 
                        x_max_discrete + 0.3* (x_2_discrete- x_max_discrete)])
        E_max *= -1


        x_1, E1 = _precise_min(U_data_slice, x, x_1_discrete, 
                [x_1_discrete - 0.3* (x_max_discrete - x_1_discrete), 
                    x_1_discrete + 0.3* (x_max_discrete - x_1_discrete)])
        
        x_2, E2 = _precise_min(U_data_slice, x, x_2_discrete, 
                [x_2_discrete - 0.3 * (x_2_discrete - x_max_discrete),
                     x_2_discrete + 0.3 * (x_2_discrete - x_max_discrete)])


    d = (x_2-x_1)/2 # half dot seperation (initial guess)
    
    e_field_x = ((E2 - E1) / (2 * e * d) ) # Electric field (initial guess)

    x_0 = (x_max + x_1 + x_2)/3             # Center of symmetry shifted 
                                                # due to electric field 

    k = 4 / d ** 2 * (2 * E_max - (E1 + E2) )
    omega_0 = np.sqrt((k + np.sqrt(k**2 - 6 * (e_field_x *e)**2 ))/(2 * m))
                                            # well width parameter 
                                                # (initial guess)

    U_0 = (E1 + E2 + E_max) /3 - (m * omega_0**2 * d**2)/24
                                            # potential at x_0 
                                                # (initial guess)

    # grouping all initial guesses
    min_params = np.array([e_field_x, omega_0, d, x_0, y_0, U_0])
    
    
    # Create function for best fit
    def quartic_function(parameters, x_grid, y_grid):
        '''
        Creates a 2D Array of fit data

        Parameters:
        ------------
        parameters: 1D Arrayof(Float)
            Array containing values that define the well (in proper order)
            normalized by initial guesses:

            np.array([e_field_x, omega_0, d, x_0, U_0]) / guess 

        Returns:
        ------------
        U_fit: 2D Arrayof(Float)
            Contains the quartic potential of best fit
            depending on x and y 
        '''

        e_field_x, omega_0, d, x_0, y_0, U_0 = parameters 
        U_fit = (m * omega_0 ** 2) / 2 *\
                (np.square(np.square(x_grid - x_0) - d ** 2) / (4 * d**2)
                                + np.square(y_grid - y_0)) +\
                                        e_field_x * e * x_grid + U_0
        return U_fit
    
    #___________________
    # HERE, a minimization routine used to be (see archive)
    #____________________

    min_e_field_x, min_omega_0, min_d, min_x_0, min_y_0, min_U_0 = min_params

    U_best_fit = quartic_function(min_params, x_mesh, y_mesh)
    
    if return_params:
        result_dict = {'U_fit': U_best_fit, \
                        'dot_sep': 2* min_d, 'e_field_x': min_e_field_x,\
                            'omega_0': min_omega_0, 'x_centre': min_x_0, \
                                'y_centre': min_y_0, 'U_0': min_U_0
                        }
        return result_dict
    else:
        return U_best_fit



def _precise_min(fun_array, arg_array, guess, bounds):
    """
    Finds minimum of a function precisely by interpolating data array and minimizing it

    Args:
        fun_array (1d float array): array of function values
        arg_array (1d float array): [description]
        guess (float): guess of an argument, near which the minimum is located
        bounds (list): 2-element list/tuple with lower and upper bound of 
        minimization
    """
    bound_mask = np.logical_and(np.greater_equal(arg_array, bounds[0]), 
                                            np.less_equal(arg_array, bounds[1]))

    # converting potential landscape in eV, and distance in nanometer
    # for more precise minimization
    fun = interp1d(1e9 * arg_array[bound_mask], 
                        1.6e19* fun_array[bound_mask], kind='cubic')

    arg_min = minimize(fun, 1e9* guess, bounds=[1e9* np.array(bounds)], 
                                    method='L-BFGS-B', options={'ftol':1e-12})
                     
    result = (arg_min.x[0] / 1e9, fun(arg_min.x[0] )/1.6e19)
    return result