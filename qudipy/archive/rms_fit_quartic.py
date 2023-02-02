

rel_err = [0]
# Create function to minimize
def minimize_func(parameters):
    '''
    Defines function to minimize

    Parameters:
    ------------
    parameters: 1D Arrayof(Float)
        Array containing values that define the well
        [e_field_x, omega_0, d, x_0, U0] - order important

    Returns:
    ------------
    Float: ||U_fit - U_data|| within bounds
        The bounds will move through the interpolation. 

    '''
    e_field_x, omega_0, d, x_0, y_0, U0 = parameters * guess
    xi = np.sqrt(hbar/ (m * omega_0))

    # coordinate bounds for fitting
    #TODO verify shift
    x_min = x_1 - e* e_field_x / (m * omega_0**2) - 0.5* xi
    x_max = x_2 - e* e_field_x / (m * omega_0**2) +  0.5* xi

    y_min = y_0 - 0.5* xi
    y_max = y_0 + 0.5* xi

    # bounding mesh arrays for fitting
    x_bool = np.logical_and(
                                np.greater_equal(x, x_min), 
                                    np.less_equal(x, x_max))
    x_bound = x[x_bool]

    y_bool = np.logical_and(
                                np.greater_equal(y, y_min), 
                                    np.less_equal(y, y_max))

    y_bound = y[y_bool]

    xy_bool = y_bool[:, np.newaxis] * x_bool

    x_mesh_bound, y_mesh_bound = np.meshgrid(x_bound, y_bound)

    U_data_bound = U_data[xy_bool].reshape(len(y_bound), len(x_bound)) 

    U_fit_bound = quartic_function(parameters, x_mesh_bound, y_mesh_bound)

    
    # Update relative error of fitting
    rel_err_temp = (np.linalg.norm(U_fit_bound - U_data_bound) / 
                                            np.linalg.norm(U_data_bound))
    
    # print('Error: ', rel_err_temp, '\t xi: ', xi)
    rel_err[0] = rel_err_temp
        

    # return np.linalg.norm(U_fit_bound - U_data_bound)
    return rel_err_temp

# ** Find Optimal Fit **
# Create initial guess for minimization: 
guess_normalized = np.ones(6, dtype=np.float64)

# Find the minimum
res = minimize(minimize_func, x0=guess_normalized, tol=1e-8)

# Return the optimized potential and parameters
min_params = res.x