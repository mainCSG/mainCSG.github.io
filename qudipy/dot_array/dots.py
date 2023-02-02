
import copy
import itertools
import os
import pickle as pk
from itertools import compress

# python modules
import matplotlib.pyplot as plt
import numpy as np
from IPython.core.pylabtools import figsize
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import find_peaks
from tqdm import tqdm

import qudipy.exchange as ex
import qudipy.potential as pot
import qudipy.qutils as qt
import qudipy.starkshift as ss
# QuDiPy modules
import qudipy.utils.helpers as hp
from qudipy.utils.constants import Constants


class DotArray:
    '''
    
    This class is designed to process the potential landscape data set for a 
    given device dot geometry to exctract various meta data or effective 
    parameters.

    Example: find the the number of dots for the potential lanscape dataset and 
    then calculate the g-factor deviation for each dot and exchange for all
    pars of adjacent dots.

    Atributes
    ----------

    Methods
    ----------

    '''
    def singleton(self, ctrl_vals):
        '''
        Identify singleton dimensions and return boolean mask
        '''
        mask = np.zeros((len(ctrl_vals)))
        for idx, dim in enumerate(ctrl_vals):
            if len(dim) != 1:
                mask[idx] = 1

        mask = list(map(bool,mask))

        return mask

    def nosingleton(self, ctrl_vals):
        '''
        Remove single dimensional data
        '''
        self.mask = self.singleton(ctrl_vals)

        new_coords = tuple(compress(ctrl_vals, self.mask))
        new_names = tuple(compress(self.ctrl_names, self.mask))

        # return new_response, new_coords, self.mask
        return new_coords, list(new_names)

    def __init__(self, precalculated_dir, data_dir, preprocessed_dir, n_dots,
                    ctrl_names, ctrl_vals, file_prefix='dots_data', 
                        calc_eff=False, save=False, pts=[6,6,6]):
        '''
        Initializes the Dots class object.
        
        Parameters
        ----------
        precalculated_dir: string
            Data directory where pre-calculated data is stored/saved.
        data_dir: string
            Directory where nextnano raw data is save
        preprocessed_dir: string
            Directory where 2D slices are loaded from or saved too.
        n_dots: int
            Expected number of dots in the system
        ctrl_names: list
            list of all control names from nextnano input script.
        ctrl_vals: list of 1d arrays
            list of control values for all control names
        
        Keyword Arguments
        -----------------
        file_prefix: string
            Prefix added to saved/loaded files 
        calc_eff: bool
            Trigger the calculation of effective parameters.
        save: bool
            Save 2D potential slices to Pre-processed directory.
        pts: list
            Number of interpolation points per control name.
   
        Returns
        -------
        None.
        '''
        self.precalculated_dir = precalculated_dir
        self.data_dir = data_dir
        self.preprocessed_dir = preprocessed_dir
        
        self.n_dots = n_dots        
        # indices of accessible dots: will change when the object 
        # is split or sliced
        self.visible_dots = list(range(1, n_dots + 1))

        self.ctrl_names = ctrl_names
        self.ctrl_vals = ctrl_vals

        self.file_prefix = file_prefix
        self.calc_eff =  calc_eff
        self.save = save
        self.pts = pts

        self.csts = Constants('Si/SiO2')

        # getting potential information 
        # create file path for potential precalculated data
        self.eff_file = os.path.join(self.precalculated_dir,
                                        f"{self.file_prefix}_eff_{self.pts}.pkl")

        self.pot_file = os.path.join(self.precalculated_dir,
                                        f"{self.file_prefix}_pot_{self.pts}.pkl")
        self.pot_int_file = os.path.join(self.precalculated_dir,
                                                    f"pot_dict_{self.pts}.pkl")
        
        #create pickle file if it doesn't exist already
        # Importing potential data or pre-calculated potential variable
        if not os.path.exists(self.pot_file):
            potential = pot.process_nextnano.import_dir(self.data_dir, show_files=True)

            print('Done')
            f = open(self.pot_file,"wb")
            # write the python object (dict) to pickle file
            pk.dump(potential,f)
            f.close()
        else:
            file_to_read = open(self.pot_file, "rb")
            potential = pk.load(file_to_read)

        # Desired Slice
        z = -0.2
        _, nearest_slice = hp.find_nearest(potential[0]['coord']['z'], z)

        print(f'The nearest slice for z = {z} is: {nearest_slice}')
        # Specifying control values and names
        pot_dir = self.preprocessed_dir + \
                    '_for_nearest_slice{:.3e}'.format(nearest_slice) + '/'

        self.ctrl_vals_nn = pot.process_nextnano.get_ctrl_vals(potential) 
                                # Must make sure that ctrl_vals are sorted
        
        # Now we define the field types we wish to write 2D slices for: 
        # either potential or the electric field.
        if self.save:
            pot.process_nextnano.write_data(self.data_dir, 
                self.preprocessed_dir, slice=z, f_type=['potential','field'])

        loaded_data_pot = pot.load_potentials(self.ctrl_vals_nn, self.ctrl_names,
                                        f_type='pot', f_dir=pot_dir,
                                        f_dis_units='nm', f_pot_units='V')

        loaded_data_field = pot.load_potentials(self.ctrl_vals_nn, self.ctrl_names,
                                        f_type='electric', f_dir=pot_dir,
                                        f_dis_units='nm', f_pot_units='V/nm')

        self.potential_data = loaded_data_pot
        self.e_field_data = loaded_data_field

        self.potential = pot.build_interpolator(loaded_data_pot, 
                                                        constants=self.csts)
        self.e_field = pot.build_interpolator(loaded_data_field, 
                                                        constants=self.csts)

        # potential interpolator of the unmasked system: never gets modified
        self.potential_unmasked = copy.deepcopy(self.potential)
        # self.e_field_unmasked = copy.deepcopy(self.e_field)

        # coordinate values
        self.x = (self.potential).x_coords
        self.y = (self.potential).y_coords

        # remove any singleton dimensions
        self.ctrl_vals, self.ctrl_names = self.nosingleton(self.ctrl_vals)

    def split(self, group='single'):
        '''
        Method that defines new Dots objects where potentials of quantum dots 
        are separated either individually or as pairs.
    
        Each of the new Dots object inherits the potential and 
        electric field data from the parent object, 
        but the interpolating function changes: when called, the potential 
        landscape appears to have only one dot (or one pair)
        
        Keyword Arguments
        -----------------
        group: string
            Whether to get a list of individual dots ('single'/'singles'), 
            or a list of adjacent dot pairs ('pair'/'pairs')
           
        Returns
        -------
        split_dots: list 
            list of Dots objects corresponding to individual dots or their pairs
        '''

        # below, updating interpolator functionality for each subgroup of dots
        if group.lower() in ('single', 'singles'):
                
            # determine the number of dot subgroups
            n_subgroups = len(self.visible_dots)

            # i=i is needed because it's essentially a nested 
            # list comprehension: with i=i, lambda uses the variable from
            # the local scope
            dot_interps = [(lambda ctrl_vals, i=i: self._separator(ctrl_vals)[i])
                                for i in range(n_subgroups)]
            visible_dots = [[n] for n in range(1, n_subgroups + 1)]

        elif group.lower() in ('pair', 'pairs'):
            
            n_subgroups = len(self.visible_dots) - 1
            # individual dot potentials are overlayed by taking element-wise 
            # minima of the arrays
            dot_interps = [(lambda ctrl_vals, i=i: 
                    np.minimum(self._separator(ctrl_vals)[i], 
                                    self._separator(ctrl_vals)[i+1]))
                                for i in range(n_subgroups)]
            visible_dots = list(zip(range(1, n_subgroups + 1), 
                                    range(2, n_subgroups + 2)))
        
        else:
            raise ValueError('Please specify a correct key for splitting '
                                'type: single(s) or pair(s)')

        # creating new identical Dots objects corresponding to 
        # the subgroups of dots
        split_dots = [copy.deepcopy(self) for n in range(n_subgroups)]

        # updating self.potential.__call__() function, and self.visible_dots
        # for each subgroup
        for idx, dot in enumerate(split_dots):
            dot.potential.call_method = dot_interps[idx]
            dot.visible_dots = visible_dots[idx]
        
        return split_dots
 
    def _separator(self, ctrl_vals):
        '''
        Function that calculates potential landscape for a given 
        list of voltage configurations, and outputs the landscapes of each 
        individual dot
         
        Parameters
        -----------------
        ctrl_vals: list     #TODO add functionality for list of lists
            Voltage configuration to calculate dot subgroup potentials at
        Returns
        -------
        dot_pots: list of 2D float arrays
            List of potential landscapes that corresponds to each individual dot
        '''
        
        # Potential landscape of the entire dot array
        array_pot = self.potential(ctrl_vals)
        # finding maxima along each  horizontal (y=y_i) slice
        pot_max = np.max(array_pot, axis=1)

        # Finding local minima and maxima
        # find indices of absolute minimum
        min_idx = np.unravel_index(array_pot.argmin(), array_pot.shape)
        y_min_idx, x_min_idx = min_idx

        # 1D slice along y=y_min_idx, where minima/maxima are to be searched for
        pot_slice =  array_pot[y_min_idx]

        # minima and maxima
        maxima = find_peaks(pot_slice)[0]
        minima = find_peaks((-1) * pot_slice)[0]

        # remove edge maxima if there are any
        if maxima[0] < minima[0]:
            maxima = np.delete(maxima, 0 )
        if maxima[-1] > minima[-1]:
            maxima = np.delete(maxima, -1)

        # If #(local minima) != n_dots, or #(local maxima) != n_dots -1, 
        # then quantum dots are not well defined,
        # and arrays with np.nans are outputted
        if (len(maxima) != (len(self.visible_dots) - 1) 
                            or len(minima) != len(self.visible_dots)):
            print('Number of actual dots is', len(minima), 
            '-> different from user-specified value', self.n_dots)
            nan_array = np.full(array_pot.shape, np.nan)
            dot_pots = [nan_array] * len(self.visible_dots)
            return dot_pots
        else:
            # filling each single dot potential with the maximum values
            # along each y-slice
            single_dot_pot = np.tile(pot_max, (len(self.x),1)).T

            # *defining boolean masks that correspond to individual dots*
            # mesh of indices along horizontal (x) axis
            x_idx_mesh = np.vstack([np.arange(len(self.x))] * len(self.y))

            # index positions of dot boundaries
            dot_bounds = np.concatenate(([0], maxima, [len(self.x)]))

            # building masks of type bound_1 <= index_along_x < bound_2
            dot_masks = [(np.greater_equal(x_idx_mesh, dot_bounds[i]) & 
                            np.less_equal(x_idx_mesh, dot_bounds[i + 1])) 
                                        for i in range(len(dot_bounds) - 1)]
            
            # make potentials the same as self.potential where dots are localized, 
            # and keep the rest of the values at value of pot_max

            dot_pots = [copy.deepcopy(single_dot_pot) for dot in range(len(self.visible_dots))] 
            for dot_pot, mask in zip(dot_pots, dot_masks):
                dot_pot[mask] = array_pot[mask]

            return dot_pots
            
    def g_factors(self, ctrl_vecs):
        '''
        Determines g-factors for each dot in the Dots object.
        Parameters
        -------
        ctrl_vals: list of lists
            Configuration of control voltages for various voltage vectors in 
            the voltage state space
        Returns
        -------
        g_factor_path: 2D float array
            Array of the g-factors of each dot for various voltage vectors in 
            the voltage state space
        '''

        # check if ctrl_vecs is a single voltage state space vector or multiple
        dim = len(np.shape(ctrl_vecs)) 
        if dim > 1:
            # multiple voltage vectors
            num_vecs = np.shape(ctrl_vecs)[0]
            g_factor_path = np.zeros((num_vecs, len(self.visible_dots)), dtype=float)
        else:
            # one voltage vector
            ctrl_vecs = [ctrl_vecs]
            g_factor_path = np.zeros((len(self.visible_dots)), dtype=float)
            

        # loop through all list elements containing points in voltage state
        # space
        for vec_idx, ctrl_vals in enumerate(ctrl_vecs):

            # number of g-factors is the same as the number of visible dots
            g_factors = np.zeros((len(self.visible_dots),), dtype=float)
            gparams = pot.GridParameters(self.x, self.y)

            # create Dots objects for individual dots
            single_dots = self.split(group='single')

            for idx, dot in enumerate(single_dots):
                gparams.update_potential(dot.potential(ctrl_vals))
                # evaluating ground state wavefunction
                _, e_vecsr = qt.solvers.solve_schrodinger_eq(self.csts, gparams, n_sols=1)
                wfr = e_vecsr[:,:,0]

                # g-factor calculation
                stark_shiftr = ss.starkshift.StarkShift(gparams, self.csts)
                dg = stark_shiftr.delta_g(self.e_field, [ctrl_vals], self.ctrl_names,
                                                    wavefuncs=[wfr])['delta_g_1'][0]
                
                # assign g-factor deviation per dot to dictionary
                g_factors[idx] = dg

            if dim > 1:
                g_factor_path[vec_idx, :] = g_factors
            else:
                g_factor_path = g_factors


        return g_factor_path

    def exchanges(self, ctrl_vecs, method='HM'):
        '''
        Determines exchange couplings between adjacent pairs of dots.

        Parameters
        -------
        ctrl_vals: list of lists
            Configuration of control voltages for various voltage vectors in 
            the voltage state space

        Keyword Arguments
        -----------------
        method: string
            Choice of 'HL', 'HM', or 'LCHO' method for exchange calculations
            (strings are case insensitive)
        Returns
        -------
        ex_path: 2D float array
            Array of the exchange couplings between each adjacent pair of dots
            for various voltage vectors in the voltage state space
        '''

        # check if ctrl_vecs is a single voltage state space vector or multiple
        dim = len(np.shape(ctrl_vecs)) 
        if dim > 1:
            # multiple voltage vectors
            num_vecs = np.shape(ctrl_vecs)[0]
            ex_path = np.zeros((num_vecs, len(self.visible_dots)-1), dtype=float)
        else:
            # one voltage vector
            ctrl_vecs = [ctrl_vecs]
            ex_path = np.zeros((len(self.visible_dots)-1), dtype=float)

        # loop through all list elements containing points in voltage state
        # space
        for vec_idx, ctrl_vals in enumerate(ctrl_vecs):

            # number of exchange parameters is less by 1 than the number of visible dots
            exchanges = np.zeros((len(self.visible_dots) - 1,), dtype=float)
            
            if len(self.visible_dots) == 1:
                print('Only one dot detected, exchange is not evaluated')
                return np.zeros((1,), dtype=float)
            else:
                gparams = pot.GridParameters(self.x, self.y)

                # create Dots objects for adjacent pairs
                dot_pairs = self.split(group='pairs')

                for idx, pair in enumerate(dot_pairs):
                    gparams.update_potential(pair.potential(ctrl_vals))
                    if method.lower() == 'hl':
                        exchanges[idx] = ex.hl_quartic(gparams, material='Si/SiO2')
                    elif method.lower() == 'hm':
                        exchanges[idx] = ex.hm_quartic(gparams, material='Si/SiO2')
                    else:
                        print('Specify correct key for exchange calculation method.'
                        'LCHO is temporarily unavailable')
                        exchanges[idx] = 0

                if dim > 1:
                    ex_path[vec_idx, :] = exchanges
                else:
                    ex_path = exchanges

        return ex_path


    def plot(self, ctrl_vecs, param='gfactor', ex_units='J', yscale='linear'):

        # convert exchange into proper units
        csts = Constants()
        if ex_units.lower() == 'j':
            units = 1
            unit_label = 'J'
        elif ex_units.lower() == 'ev':
            units = csts.e
            unit_label = 'eV'
        elif ex_units.lower() == 'mev':
            units = csts.e * 1e-3
            unit_label = 'meV'
        elif ex_units.lower() == 'uev':
            units = csts.e * 1e-6
            unit_label = '$\mu$eV'
        elif ex_units.lower() == 'nev':
            units = csts.e * 1e-9
            unit_label = 'neV'
        
        # y-scale
        if yscale.lower() not in ['linear', 'log']:
            Warning(f'y-scale: {yscale} is not a valid scale for y-axis.'
                        'Linear is used by default')

        # collect all varying controls, keep a singe value for the static ones
        # (preserving their indices)

        all_ctrls = list(enumerate(zip(*ctrl_vecs)))
        var_ctrls = [(idx, ctrl) for idx, ctrl in all_ctrls
                                     if ctrl.count(ctrl[0]) != len(ctrl)]
        static_ctrls = [(idx, ctrl[0])  for idx, ctrl in all_ctrls
                                     if ctrl.count(ctrl[0]) == len(ctrl)]
        # grouping index positions and values separately
        var_ids, var_vals = zip(*var_ctrls) 
        static_ids, static_vals = zip(*static_ctrls)
        # names of variable and constant voltage controls
        var_names = [self.ctrl_names[idx] for idx in var_ids]
        static_names = [self.ctrl_names[idx] for idx in static_ids]

        # constructing annotation with the values of static parameters
        annot_text = (''.join(f'{name} = {val} V\n' for name, val
                                in zip(static_names, static_vals)))[:-1]
                                        # [:-2] removes the extra \n at the end

        # 1D plot in case *only one* control parameter is varied
        # otherwise, 2D (3D, ...) plot must be implemented (TBD) 
        if len(var_ctrls) == 1:
            fig, ax = plt.subplots(1,1, figsize=(6, 4), dpi=120)

            # g-factor
            if param.lower() in ('gfactor', 'g_factor', 'g-factor'):
                g_factor_path = self.g_factors(ctrl_vecs)
                for g in range(np.shape(g_factor_path)[1]):
                    # the only varied control on the horizontal axis
                    ax.plot(var_vals[0], g_factor_path[:, g], '-', 
                                                label=f'Dot {g+1}')
                ax.set_ylabel('$\delta g$ deviation')
                ax.set_title('$g-$factor deviation per dot')

            # exchange
            elif param.lower() in ('hl', 'hl-exchange', 'hl_exchange', 
                                                            'hl-ex','hl_ex'):
                ex_path = self.exchanges(ctrl_vecs, method='hl')

                for j in range(np.shape(ex_path)[1]):
                    ax.plot(var_vals[0], ex_path[:, j] / units, '-', 
                                                label=f'Dot pair {j+1} & {j+2}')
                ax.set_ylabel(f'J [{unit_label}]')
                ax.set_title('Heitler-London Exchange')

            elif param.lower() in ('hm', 'hm-exchange', 'hm_exchange', 
                                                                'hm-ex','hm_ex'):
                ex_path = self.exchanges(ctrl_vecs, method='hm')

                for j in range(np.shape(ex_path)[1]):
                    ax.plot(var_vals[0], ex_path[:,j] / units, '-', 
                                                label=f'Dot pair {j+1} & {j+2}')
                ax.set_ylabel(f'J [{unit_label}]')
                ax.set_title('Hund-Mulliken Exchange')

            ax.annotate(annot_text, (0.07, 0.6), xycoords='axes fraction',
                                 fontsize=10)
            ax.set_xlabel(f'{var_names[0]} [V]')
            ax.set_yscale(yscale)
            ax.legend(loc='upper left')


            


    def get_eff_params(self, exchange_method='HM'):
        '''
        Method that calculates the masked potentials and effective parameters
        for all simulation data. Then interpolation objects are computed.
        
        Keyword Arguments
        -----------------
        exchange_method: string
            Choice of HL, HM, or LCHO method for exchange calculations
           
        Returns
        -------
        None.
        '''
        if not os.path.exists(self.pot_int_file) or self.calc_eff:
            param_list = list(itertools.product(*self.ctrl_vals))

            # Determine data set size for array initialization
            tup = ()
            for i in self.ctrl_vals:
                tup += (len(i),)

            self.gl_data = np.empty(tup)
            self.gr_data = np.empty(tup)
            self.exchange = np.empty(tup)

            # temporary assignment to be assigned to each dot that is addressed
            self.dg = np.empty(tup)

            # Define gridparameters for dot array and subdot array
            # _2d_coords = loaded_data_pot['coords']

            x = (self.potential).x_coords
            y = (self.potential).y_coords

            gparams = pot.GridParameters(x,y)

            # initialize delta-g, exchange dictionaries and mask potential -----
            #  dictionary 
            self.eff_dict_vals = {}
            dg_keys, J_keys = [], []
            ctrl_coord = []

            # initialize dictionary to store 2d potentials per
            # control coordinate
            pot_dict = {
                            
                            'coords': {0: x, 1: y},
                            'potentials' : [],
                            'ctrl_vals' : [],
                            'ctrl_names' : self.ctrl_names
                        }
            self.pot_int_dict = {'unmasked': pot_dict}

            for idx, param in enumerate(tqdm(param_list)):

                # TODO: update the above code
                # make each tuple of values a list for later computations
                param = list(param)

                # Copy potential data for later use
                array_pot = self.potential(param).copy()
                pot_min = (-1)*self.potential(param).copy()

                # store 2d potential for given control coordinate
                self.pot_int_dict['unmasked']['potentials'].append(array_pot)
                self.pot_int_dict['unmasked']['ctrl_vals'].append(param)

                # save ctrl values
                ctrl_coord.append(param)
                
                # ** find the maximum **
                # find absolute minimum
                min_idx = np.unravel_index(array_pot.argmin(), array_pot.shape)
                # draw a y_0 through the minimum
                y_0 = y[min_idx[0]]     
                U_data_slice =  array_pot[min_idx[0]]
                # maxima
                maximas = find_peaks(U_data_slice)[0]

                # ** find the minima **            
                U_data_slice =  pot_min[min_idx[0]]
                # maxima
                minimas = find_peaks(U_data_slice)[0]

                # TODO: track largest number of dots observed over entire param set
                # number of dots
                n_dots = len(minimas)

                # ignore erroneous maximas at the ends of dot arrays
                maximas = [x for x in maximas if minimas[0] < x and minimas[-1] > x]

                # initialization of empty dictionaries for delta-g, exchange,
                # masked potentials
                if idx == 0:
                    for i in range(n_dots):
                        # delta-g value
                        key = f'delta_g{i+1}'
                        dg_keys.append(key)
                        self.eff_dict_vals[key] = []
                        # single dot masked potential
                        self.pot_int_dict[f'masked_{key}'] = {
                                        
                                        'coords': {0: x, 1: y},
                                        'potentials' : [],
                                        'ctrl_vals' : [],
                                        'ctrl_names' : self.ctrl_names
                                    }

                    for i in range(n_dots-1):                        
                        # exchange value
                        key = f'J_{i+1}{i+2}'
                        J_keys.append(key)
                        self.eff_dict_vals[key] = []
                        # double dot masked potentials
                        self.pot_int_dict[f'masked_{key}'] = {
                                        
                                        'coords': {0: x, 1: y},
                                        'potentials' : [],
                                        'ctrl_vals' : [],
                                        'ctrl_names' : self.ctrl_names
                                    }

                # calculate masked potentials and effective parameters ---------
                if (len(maximas) != 1):
                    print("NaNs outputted by g-factor and exchange")

                    # build 2D nan array to store for potential interpolators
                    nan_pot = np.empty((len(y),len(x)))
                    nan_pot[:] = np.nan

                    # TODO: make sure to assign nans for maximum number of 
                    # dots previously observed 
                    # (HANDLE IF FIRST PARAM CAUSES LESS THAN EXPECTED TOTAL DOTS)
                    for key in dg_keys:
                        # store nans for dot masked potentials
                        self.pot_int_dict[f'masked_{key}']['potentials'].append(nan_pot)
                        self.pot_int_dict[f'masked_{key}']['ctrl_vals'].append(param)
                        # store nans for g-factor deviation
                        self.eff_dict_vals[key].append(np.nan)

                    for key in J_keys:
                        # store nans for dot masked potentials
                        self.pot_int_dict[f'masked_{key}']['potentials'].append(nan_pot)
                        self.pot_int_dict[f'masked_{key}']['ctrl_vals'].append(param)
                        # store nans for exchange
                        self.eff_dict_vals[key].append(np.nan)

                else:
                    # calculated single dot masked potentials and g-factor
                    # deviation
                    for dg_idx in range(n_dots):
                            # copy unmasked potential
                            gparams = pot.GridParameters(x,y)
                            dot_pot = self.potential(param).copy()
                            # masked single dot potentials
                            base = np.amax(dot_pot)
                            if dg_idx == 0:
                                # define left dot ------------------
                                dot_pot[:,maximas[dg_idx]:] = base
                            elif dg_idx == n_dots-1:
                                # define right dot -----------------
                                dot_pot[:,:maximas[dg_idx-1]] = base
                            else:
                                # define middle dots ---------------
                                dot_pot[:,:maximas[dg_idx-1]] = base
                                dot_pot[:,maximas[dg_idx]:] = base

                            # store single dot masked potentials
                            gparams.update_potential(dot_pot)
                            self.pot_int_dict[f'masked_{dg_keys[dg_idx]}']['potentials'].append(dot_pot)
                            self.pot_int_dict[f'masked_{dg_keys[dg_idx]}']['ctrl_vals'].append(param)

                            if self.calc_eff:
                                
                                _, e_vecsr = qt.solvers.solve_schrodinger_eq(self.csts, gparams, n_sols=1)
                                wfr = e_vecsr[:,:,0]

                                # g-factor
                                stark_shiftr = ss.starkshift.StarkShift(gparams, self.csts)
                                dg = stark_shiftr.delta_g(self.e_field, [param], self.ctrl_names,
                                                                                wavefuncs=[wfr])['delta_g_1'][0]
                                
                                # assign g-factor deviation per dot to dictionary
                                self.eff_dict_vals[dg_keys[dg_idx]].append(dg)

                    # calculate masked potentials for adjacent pair of dots and
                    # exchange
                    for j_idx in range(n_dots-1):
                            # copy unmasked potential
                            gparams = pot.GridParameters(x,y)
                            dot_pot = self.potential(param).copy()

                            # must not be a single dot
                            if len(maximas) > 1:
                                # masked double dot potentials
                                base = np.amax(dot_pot)
                                if j_idx == 0:
                                    # define left dot ------------------
                                    dot_pot[:,maximas[j_idx+1]:] = base
                                elif j_idx == n_dots-2:
                                    # define right dot -----------------
                                    dot_pot[:,:maximas[j_idx-1]] = base
                                else:
                                    # define middle dots ---------------
                                    dot_pot[:,:maximas[j_idx-1]] = base
                                    dot_pot[:,maximas[j_idx+1]:] = base

                            # store double dot masked potentials
                            gparams.update_potential(dot_pot)
                            self.pot_int_dict[f'masked_{J_keys[j_idx]}']['potentials'].append(dot_pot)
                            self.pot_int_dict[f'masked_{J_keys[j_idx]}']['ctrl_vals'].append(param)

                            if self.calc_eff:
                                # assign exchange per double dot pair to dictionary
                                if exchange_method.lower() == 'hl':
                                    self.eff_dict_vals[J_keys[j_idx]].append(ex.hl_quartic(gparams, material='Si/SiO2'))
                                elif exchange_method.lower() == 'hm':
                                    self.eff_dict_vals[J_keys[j_idx]].append(ex.hm_quartic(gparams, material='Si/SiO2'))

            # Create dictionary of masked potential interpolation objects
            self.pot_interp_dict = {}
            for key in self.pot_int_dict.keys():
                self.pot_interp_dict[key] = pot.build_interpolator(self.pot_int_dict[key])

            f = open(self.pot_int_file,"wb")
            # write the python object (dict) to pickle file
            pk.dump(self.pot_interp_dict,f)
            f.close()


            if self.calc_eff:
                # Create dictionary of interpolation objects for exchange and 
                # delta-g
                self.eff_interp_dict = {}
                for key in self.eff_dict_vals.keys():
                    self.eff_interp_dict[key] = LinearNDInterpolator(ctrl_coord,self.eff_dict_vals[key])
                    
                self.eff_interp_dict['ctrl_vals'] = self.ctrl_vals

                f = open(self.eff_file,"wb")
                # write the python object (dict) to pickle file
                pk.dump(self.eff_interp_dict,f)
                f.close()

        else:

            file_to_read = open(self.pot_int_file, "rb")
            self.pot_interp_dict = pk.load(file_to_read)

            if os.path.exists(self.eff_file):
                file_to_read = open(self.eff_file, "rb")
                self.eff_interp_dict = pk.load(file_to_read)
