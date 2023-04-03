import astropy.constants as c
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import yaml

class C_parameters:
    """Read in yml parameters and store them.
    
    Simple class to read in and store parameters from the yml file. 
    Simple methods included to print out the parameters etc.
    All global parameters are stored here.

    Attributes
    ----------
    b_* : boolean
        Flag for each of the model parameters
    baryon_fraction : float
        Cosmic baryon fraction.
    cooling_table: C_cooling
        Stores cooling table and associated metal and temperature scales
    c_* : float
        Various dimensionless physical constants for use in astrophysics routines
    D_param : dictionary 
        Dictionary containing contents of yml file.
    n_HDF5_io_rec : int
        IO HDF5 buffer size.
 
    Methods:
    --------
        __init__
        __str__ 
    
    """
    
    def __init__(self,param_file):
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
        # Diagnostics
        for key, value in self.D_param['diagnostics'].items():
            Value = value['Value']
            exec('self.'+str(key)+'=Value')
            print('self.'+str(key)+' =',eval('self.'+str(key)))
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
        # Should really label these as _internal, but that makes the names too long.
        # Also, parameters that are not read in do not have the internal label.
        # However, it feels wrong to redefine input variables, so will stick with adding _internal label to input parameters.
        # This should only arise for quantities in the parameter class instance.
        # These are the minimum allowable timesteps for mini-steps
        self.timestep_halo_internal = (self.timestep_halo / self.units_time_internal).si.value
        self.timestep_gal_internal = (self.timestep_galaxies / self.units_time_internal).si.value
        # The minimum mass (values below this are regarded as equivalent to zero mass)
        self.mass_minimum_internal = (self.mass_minimum/self.units_mass_internal).si.value
        # The minimum length (values below this are regarded as not set)
        self.length_minimum_internal = (self.length_minimum/self.units_length_internal).si.value
        # Quantities that need to be in internal code units for use in functions
        self.Hen15_v_eject_internal = (self.Hen15_v_eject/self.units_speed_internal).si.value
        self.Hen15_v_reheat_internal = (self.Hen15_v_reheat/self.units_speed_internal).si.value
        self.BH_v_q_internal = (self.BH_v_q/self.units_speed_internal).si.value
        # Reincoportation. The following gets divided by mass, hence the units
        self.c_Hen15_reinc = (self.Hen15_gamma_reinc*1e10*c.M_sun/(self.units_time_internal*self.units_mass_internal)).si.value

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
        self.c_sfr_Mcrit = (self.sfr_Mcrit0 / self.units_mass_internal) * (self.units_speed_internal * u.s / (200. *u.km)) * \
                          (self.units_length_internal / (10. * u.kpc))
        print('c_sfr_Mcrit = {:.3g}'.format(self.c_sfr_Mcrit.si))
        self.c_sfr_Mcrit = self.c_sfr_Mcrit.si.value
        # Heating energy from Hen15 (arXiv:1410.0365) equation S16.
        self.c_Hen15_S16 = self.units_mass_internal*self.Hen15_v_snr**2/(2.*self.units_energy_internal)
        print('c_Hen15_S16 = {:.3g}'.format(self.c_Hen15_S16.si))
        self.c_Hen15_S16 = self.c_Hen15_S16.si.value
        # BH radio accretion rate from Hen15 eq S24
        self.c_BH_r = self.BH_f_r * self.units_mass_internal*self.units_time_internal/(1e11*1e8*c.M_sun**2)
        print('c_BH_r = {:.3g}'.format(self.c_BH_r.si))
        self.c_BH_r = self.c_BH_r.si.value
        # BH cooling suppresion from Hen15 eq S25/S26: 0.2 is 0.1 / 0.5
        self.c_BH_mheat_r = 0.2 * (c.c/self.units_speed_internal)**2
        print('c_BH_mheat_r = {:.3g}'.format(self.c_BH_mheat_r.si))
        self.c_BH_mheat_r = self.c_BH_mheat_r.si.value
        # Determination of Eddington accretion rate (note: strictly should be m_H)
        self.c_BH_Edd = 4.*np.pi*c.G*c.m_p*self.units_time_internal/(c.sigma_T*c.c)
        print('c_BH_Edd = {:.3g}'.format(self.c_BH_Edd.si))
        self.c_BH_Edd = self.c_BH_Edd.si.value
        
        if self.b_display_parameters: self.__str__()
            
    def __str__(self):
        """ 
        Simple method to print out parameters.  
        Returns 
        """
        for item in self.D_param:
            print("{:20s}: {}\n".format(item,self.D_param[item]))
        print('\n')
        # Double down on listing units as super-useful
        print('Internal unit:')
        print('  mass  ',self.units_mass_internal.to(c.M_sun))
        print('  length',self.units_length_internal)
        print('  time  ',self.units_time_internal)
        print('  speed ',self.units_speed_internal.to(u.km/u.s))
        print('  temp. ',self.units_temperature_internal)
        return ''

    def F_update_parameters(self,graph_file):
        for key, value in graph_file['Header'].attrs.items():
            if self.b_display_parameters: print(key,value)
            exec('self.'+key+'=value')
        # Store the total number of graphs in the input file:
        self.n_graph=len(graph_file['graph_lengths'])
        # Put code in here to either copy table of snapshot redshifts/times from graph_file,
        # Or calculate them if that does not exist.
        # Currently read from disk:
        self.snap_table=np.loadtxt(self.snap_file,usecols=[0,2,4],
            dtype=[('snap_ID',np.int32),('redshift',np.float32),('time_in_years',np.float32)])


