import h5py
import numpy as np 
import os
import sys 

from chimes_utils import set_initial_chemistry_abundances 
from phys_const import proton_mass_cgs, boltzmann_cgs, seconds_in_a_Myr 
from shielding_utils import compute_jeans_shield_length, compute_colibre_shield_length 

class SnapshotData: 
    def __init__(self, driver_pars, global_pars, gas_pars): 
        self.driver_pars = driver_pars 
        self.global_pars = global_pars 
        self.gas_pars = gas_pars 

        # Declare particle arrays that 
        # will need to be read in. 
        self.nH_arr = None 
        self.temperature_arr = None 
        self.metallicity_arr = None 
        self.shieldLength_arr = None 
        self.init_chem_arr = None 
        self.ChimesFluxIon_arr = None 
        self.ChimesFluxG0_arr = None 
        self.gas_coords_arr = None 
        self.star_coords_arr = None 
        self.star_mass_arr = None 
        self.star_age_Myr_arr = None 
        self.HIIregion_delay_time = None 
        
        return 

    def set_shielding_array(self): 
        if self.driver_pars['shield_mode'] == None: 
            self.shieldLength_arr = np.ones(len(self.nH_arr), dtype = np.float64) 
        elif self.driver_pars['shield_mode'] == "read-in": 
            # Reads in the shielding column density, from which 
            # it calculates the corresponding shielding lengths. 
            with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
                self.shieldLength_arr = np.array(h5file[self.driver_pars["snapshot_column_density_array"]])/ self.nH_arr 
        elif self.driver_pars['shield_mode'] == "Jeans": 
            # Compute Jeans length assuming 
            # hydrogen mass fraction XH = 0.7 
            # mean molecular weight mu = 1 
            self.shieldLength_arr = np.zeros(len(self.nH_arr), dtype = np.float64) 

            for i in range(len(self.nH_arr)): 
                self.shieldLength_arr[i] = compute_jeans_shield_length(self.temperature_arr[i], self.nH_arr[i], self.driver_pars["shield_length_factor"], self.driver_pars["max_shield_length"]) 
                
            print('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
            print('self.shieldLength_arr[0] = ', self.shieldLength_arr[0])
            print('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
                
        elif self.driver_pars['shield_mode'] == "Colibre": 
            self.shieldLength_arr = np.zeros(len(self.nH_arr), dtype = np.float64) 
            
            for i in range(len(self.nH_arr)): 
                XH = 1.0 - (self.metallicity_arr[i, 0] + self.metallicity_arr[i, 1]) 
                self.shieldLength_arr[i] = compute_colibre_shield_length(self.temperature_arr[i], self.nH_arr[i], XH, self.driver_pars["shield_length_factor"], self.driver_pars["max_shield_length"], self.driver_pars["colibre_log_T_min"], self.driver_pars["colibre_log_T_max"], self.driver_pars["colibre_use_turbulent_jeans_length"], self.driver_pars["colibre_sigma_turb"])
        else: 
            raise Exception("ERROR: shield_mode %d not recognised. Aborting." % (self.driver_pars['shield_mode'], ))

        # HII regions 
        if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
            ind_HII = (self.HIIregion_delay_time > 0.0) 
            self.shieldLength_arr[ind_HII] = 1.0 
        
        return



    def load_USER(self): 
        with h5py.File(self.driver_pars['input_file'], 'r') as h5file:

            unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"]
            
            print()
            print('ccccccccccccccccccccccccccccccccccccccccc')
            print('unit_length_in_cgs = ', unit_length_in_cgs)
            print('ccccccccccccccccccccccccccccccccccccccccc')
            print()
            

            print("Reading in particle data\n" )
            sys.stdout.flush()
            
            # Read in metallicity array. 
            # Given as mass fractions relative 
            # to total in the order: 
            # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
            self.metallicity_arr = np.array(h5file['PartType0/Metallicity']) # metallicity_arr shape is (N_part, 11)
            
            #N_part = len(GFM_Z) 
            #self.metallicity_arr = np.zeros((N_part, 11)) 
            #self.metallicity_arr[:, 0] = GFM_Z 
            #self.metallicity_arr[:, 1] = GFM_metals[:, 1] 
            #self.metallicity_arr[:, 2] = GFM_metals[:, 2] 
            #self.metallicity_arr[:, 3] = GFM_metals[:, 3] 
            #self.metallicity_arr[:, 4] = GFM_metals[:, 4] 
            #self.metallicity_arr[:, 5] = GFM_metals[:, 5] 
            #self.metallicity_arr[:, 6] = GFM_metals[:, 6] 
            #self.metallicity_arr[:, 7] = GFM_metals[:, 7] 
            #self.metallicity_arr[:, 10] = GFM_metals[:, 8] 

            # nH array 
            self.nH_arr = np.array(h5file['PartType0/nHG_hfv']) #!!!!!!!!!!!!!

            # Read in the initial CHIMES abundance array
            self.init_chem_arr = np.array(h5file['PartType0/InitIonState_hfv']) #!!!!!!!!!!!!!
            
            # Temp array
            self.temperature_arr = np.array(h5file['PartType0/TempG_hfv']) #!!!!!!!!!!!!! 
                    
            if self.driver_pars["UV_field"] == "S04": 
                self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) #* unit_length_in_cgs  # gas_coords_arr shape is (N_part, 3)

            # Set the shielding length array 
            self.set_shielding_array() # Use the same approach as before for the Shielding Length (No need to do any thing here!)!

        return 



        
# Routine to check whether output 
# arrays already exist in the 
# output file, for IO_mode == snapshot. 
def snapshot_check_output_arrays(driver_pars): 
    return_value = False 

    if os.path.exists(driver_pars["output_file"]): 
        with h5py.File(driver_pars["output_file"], 'r') as h5file: 
            if driver_pars["driver_mode"] == "eqm_state": 
                array_name_list = ["%s/EqmChemistryAbundances" % (driver_pars["hdf5_output_group"], )] 
            elif driver_pars["driver_mode"] == "cooling_rates": 
                array_name_list = ["%s/log_cooling_rate" % (driver_pars["hdf5_output_group"], ), "%s/log_heating_rate" % (driver_pars["hdf5_output_group"], )] 
            elif driver_pars["driver_mode"] == "noneq_evolution": 
                array_name_list = ["%s/AbundanceEvolution" % (driver_pars["hdf5_output_group"], ), "%s/TemperatureEvolution" % (driver_pars["hdf5_output_group"], ), "%s/TimeArray_seconds" % (driver_pars["hdf5_output_group"], )]
            else: 
                raise Exception("ERROR: IO_mode %s cannot be run with driver_mode %s. Aborting." % (driver_pars["IO_mode"], driver_pars["driver_mode"])) 

            if driver_pars["UV_field"] == "StellarFluxes" and driver_pars["compute_stellar_fluxes"] == 1: 
                if driver_pars["snapshot_flux_ion_array"] != None: 
                    array_name_list.append(driver_pars["snapshot_flux_ion_array"]) 

                if driver_pars["snapshot_flux_G0_array"] != None: 
                    array_name_list.append(driver_pars["snapshot_flux_G0_array"]) 

            for array_name in array_name_list: 
                node_test = array_name in h5file 
                
                if node_test == True: 
                    return_value = True
                    print("%s already present in file %s" % (array_name, driver_pars["output_file"]))
                    sys.stdout.flush()

    return return_value 
