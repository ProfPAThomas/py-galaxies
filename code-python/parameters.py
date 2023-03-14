import astropy.constants as c
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import yaml

class C_parameters:
    """Read in yml parameters and store them.
    
    Simple class to read in and store parameters from the yml file. Simple
    methods included to print out the parameters etc.

    Data attributes:
    ----------------
    b_* : boolean
        Flag for each of the model parameters
    baryon_fraction : float
        Cosmic baryon fraction.
    D_param : dictionary 
        Dictionary containing contents of yml file.
    galaxy_file : str
        The filepath to the halo output HDF5 file.
    graph_file : str
        The filepath to the input graph HDF5 file.
    halo_file : str
        The filepath to the galaxy output HDF5 file.
    n_HDF5_io_rec : int
        IO HDF5 buffer size.
 
    Methods:
    --------
        __init__
        __str__ 
    
    """
    
    def __init__(self,param_file,available_option_file):
        """ 
        Key parameters for the model
        
        Input parameters
        ----------------
        param_file : str
           Filepath to the yml file containing the model parameters.
        verbosity : int
           The level of detail needed for debugging messages.
        b_debug : bool
           To print out debugging messages, true or false. 
        """
        
        # Read in yaml file and save entries within this parameter instance
        self.param_file = param_file
        self.D_param = yaml.load(open(param_file),Loader=yaml.Loader)
        # I/O files
        for key, value in self.D_param['io_files'].items():
            for file_type, file_name in value.items():
                exec('self.'+str(file_type)+'=file_name')
                print('self.'+str(file_type)+' =',eval('self.'+str(file_type)))
        # Cosmological parameters (ideally would be in graph input files)
        for key, value in self.D_param['cosmology'].items():
            Value = value['Value']
            Units = value['Units']
            if value['Units']=='None':
                exec('self.'+str(key)+'=Value')
            else:
                exec('self.'+str(key)+'=Value*eval(Units)')
            print('self.'+str(key)+' =',eval('self.'+str(key)))
        cosmology = FlatLambdaCDM(H0=self.H0, Om0=self.omega_m, Tcmb0=2.725)
        # Halo model
        for key, value in self.D_param['halo_model'].items():
            Value = value['Value']
            exec('self.'+str(key)+'=Value')
            print('self.'+str(key)+' =',eval('self.'+str(key)))
        # Performance
        for key, value in self.D_param['performance'].items():
            Value = value['Value']
            exec('self.'+str(key)+'=Value')
            print('self.'+str(key)+' =',eval('self.'+str(key)))        #
        # Units (see also conversion factors below)
        for io_type, value in self.D_param['units'].items():
            for quantity, props in value.items():
                Value = props['Value']
                Units = props['Units']
                if Units=='None':
                    exec('self.units_'+quantity+'_'+io_type+'=Value')
                else:
                    exec('self.units_'+quantity+'_'+io_type+'=Value*eval(Units)')
                print('self.units_'+quantity+'_'+io_type+' =',eval('self.units_'+quantity+'_'+io_type))
        # Model switches
        for key in self.D_param['model_switches']:
            exec('self.'+key+'='+str(self.D_param['model_switches'][key]['Value']))
            print('self.'+str(key)+' =',eval('self.'+str(key)))
        # Model parameters
        # This loops over all astrophysical paremeters, extracting them.  We could do the same for other yaml blocks.
        for key, value in self.D_param['model_parameters'].items():
            Value = value['Value']
            Units = value['Units']
            if Units=='None':
                exec('self.'+str(key)+'=Value')
            else:
                exec('self.'+str(key)+'=Value*eval(Units)')
            print('self.'+str(key)+' =',eval('self.'+str(key)))
        # Astrophysical parameters
        for key, value in self.D_param['astrophysics'].items():
            Value = value['Value']
            Units = value['Units']
            if Units=='None':
                exec('self.'+str(key)+'=Value')
            else:
                exec('self.'+str(key)+'=Value*eval(Units)')
            print('self.'+str(key)+' =',eval('self.'+str(key)))

        # Loop over all available options, creating False flag for each missing one
        self.D_option = yaml.load(open(available_option_file),Loader=yaml.Loader)
        for key in self.D_option:
            try:
                exec('self.'+key)
            except:
                exec('self.'+key+'='+str(self.D_option[key]['Value']))

        if self.b_display_parameters: self.__str__()
            
        # Set up units used for I/O and within the code
        # Conversion of masses for I/O
        self.mass_input_to_internal=(self.units_mass_input/self.units_mass_internal).si.value
        self.mass_internal_to_output=(self.units_mass_internal/self.units_mass_output).si.value
        # Conversion of lengths for I/O
        self.length_input_to_internal=(self.units_length_input/self.units_length_internal).si.value
        self.length_internal_to_output=(self.units_length_internal/self.units_length_output).si.value
        # Conversion of times for I/O
        self.time_input_to_internal=(self.units_time_input/self.units_time_internal).si.value
        self.time_internal_to_output=(self.units_time_internal/self.units_time_output).si.value
        # Speed unit and conversion of speed for I/O
        self.units_speed_internal = self.units_length_internal / self.units_time_internal
        self.speed_input_to_internal=(self.units_speed_input/self.units_speed_internal).si.value
        self.speed_internal_to_output=(self.units_speed_internal/self.units_speed_output).si.value
        # Conversion of temperatures for I/O
        self.temperature_input_to_internal=(self.units_temperature_input/self.units_temperature_internal).si.value
        self.temperature_internal_to_output=(self.units_temperature_internal/self.units_temperature_output).si.value
        # Other unit combinations used in the code
        self.units_energy_internal = self.units_mass_internal*self.units_speed_internal**2
        self.units_density_internal = self.units_mass_internal/self.units_length_internal**3
        self.units_lambda_internal = self.units_energy_internal * self.units_length_internal**3 / self.units_time_internal
        
        # Dimensionless versions of constants used in the code.
        # This is the minimum allowable timestep for mini-steps
        self.timestep_halo = (self.timestep_halo / self.units_time_internal).si.value
        self.timestep_gal = (self.timestep_galaxies / self.units_time_internal).si.value
        self.mass_minimum = (self.mass_minimum/self.units_mass_internal).si.value

        # Gravitational constant: units L^3/T^2M
        self.c_G=(c.G/self.units_length_internal**3*self.units_time_internal**2*self.units_mass_internal).si.value
        # To convert 3-d rms_speed to virial temperature T=mumH*sigma^2/k_B=mumH*sigma_3D^2/3k_B
        self.c_rms_speed_to_temperature = (self.mumH * self.units_speed_internal**2 / 
                                           (3. * c.k_B * self.units_temperature_internal ) ).si.value
        self.c_half_mass_virial_speed_to_temperature = (self.mumH * self.units_speed_internal**2 / 
                                           (2. * c.k_B * self.units_temperature_internal ) ).si.value
        # Constant used in cooling model as described in cooling.py.  Should be dimensionless.
        # If halo model is MEGA use 80 pi mumH k_B L_unit**3 T_unit / (M_unit Lambda_ratio Lambda_unit t_unit)
        # In L-Galaxies mode, replace 80 with 3/5 * 80 = 48 (energy rather than enthalpy).
        if self.b_lgalaxies:
            const = 48.
        else:
            const = 80.
        self.c_cooling = const * np.pi * self.mumH * c.k_B * self.units_temperature_internal / \
                          (self.units_density_internal * self.lambda_ratio * self.units_lambda_internal * self.units_time_internal)
        print('c_cooling = {:.3g}'.format(self.c_cooling.si))
        self.c_cooling = self.c_cooling.si.value
        # Star formation rate model from from Hen15 (arXiv:1410.0365) S1.6.
        # This constant should be dimensionless.
        self.c_sfr_Mcrit = (self.sfr_Mcrit0 / self.units_mass_internal) * (self.units_speed_internal * u.s / (200. *u.km)) * \
                          (self.units_length_internal / (10. * u.kpc))
        print('c_sfr_Mcrit = {:.3g}'.format(self.c_sfr_Mcrit.si))
        self.c_sfr_Mcrit = self.c_sfr_Mcrit.si.value


    def __str__(self):
        """ 
        Simple method to print out parameters.  
        Returns 
        """
        print('Inbuilt options:')
        for item in self.D_option:
            print("{:20s}: {}\n".format(item,self.D_option[item]))
        print('\nRuntime options:')
        for item in self.D_param:
            print("{:20s}: {}\n".format(item,self.D_param[item]))
        print('\n')
        return ''
    
