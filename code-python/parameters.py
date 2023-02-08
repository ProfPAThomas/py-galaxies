import astropy.constants as c
import astropy.units as u
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
    galaxy_output_file : str
        The filepath to the halo output HDF5 file.
    graph_input_file : str
        The filepath to the input graph HDF5 file.
    halo_output_file : str
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
        
        self.param_file = param_file
        self.D_param = yaml.load(open(param_file),Loader=yaml.Loader)

        # extract key variables for ease of later use
        self.graph_input_file = self.D_param['input_files']['graph_file']
        self.snap_input_file = self.D_param['input_files']['snap_file']
        self.halo_output_file = self.D_param['output_files']['halo_file']
        self.subhalo_output_file = self.D_param['output_files']['subhalo_file']
        self.galaxy_output_file = self.D_param['output_files']['galaxy_file']
        self.omega_m = self.D_param['cosmology']['omega_m']['Value']
        self.baryon_fraction = self.D_param['cosmology']['baryon_fraction']['Value']
        self.n_HDF5_io_rec = self.D_param['performance']['n_HDF5_io_rec']['Value']
        #   self.sub_halo = self.D_param['model_switches']['sub_halo']['Value']

        # Loop over model switches, creating a boolean flag for each
        for key in self.D_param['model_switches']:
            exec('self.'+key+'='+str(self.D_param['model_switches'][key]['Value']))

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
        self.units_mass_input = self.D_param['units']['input']['mass']['Value'] * eval(self.D_param['units']['input']['mass']['Units'])
        self.units_mass_internal = self.D_param['units']['internal']['mass']['Value'] * eval(self.D_param['units']['internal']['mass']['Units'])
        self.units_mass_output = self.D_param['units']['output']['mass']['Value'] * eval(self.D_param['units']['output']['mass']['Units'])
        self.mass_input_to_internal=(self.units_mass_input/self.units_mass_internal).si.value
        self.mass_internal_to_output=(self.units_mass_internal/self.units_mass_output).si.value
        # Conversion of lengths for I/O
        self.units_length_input = self.D_param['units']['input']['length']['Value'] * eval(self.D_param['units']['input']['length']['Units'])
        self.units_length_internal = self.D_param['units']['internal']['length']['Value'] * eval(self.D_param['units']['internal']['length']['Units'])
        self.units_length_output = self.D_param['units']['output']['length']['Value'] * eval(self.D_param['units']['output']['length']['Units'])
        self.length_input_to_internal=(self.units_length_input/self.units_length_internal).si.value
        self.length_internal_to_output=(self.units_length_internal/self.units_length_output).si.value
        # Conversion of times for I/O
        self.units_time_input = self.D_param['units']['input']['time']['Value'] * eval(self.D_param['units']['input']['time']['Units'])
        self.units_time_internal = self.D_param['units']['internal']['time']['Value'] * eval(self.D_param['units']['internal']['time']['Units'])
        self.units_time_output = self.D_param['units']['output']['time']['Value'] * eval(self.D_param['units']['output']['time']['Units'])
        self.time_input_to_internal=(self.units_time_input/self.units_time_internal).si.value
        self.time_internal_to_output=(self.units_time_internal/self.units_time_output).si.value
        # Speed unit and conversion of speed for I/O
        self.units_speed_input = self.D_param['units']['input']['speed']['Value'] * eval(self.D_param['units']['input']['speed']['Units'])
        self.units_speed_internal = self.units_length_internal / self.units_time_internal
        self.units_speed_output = self.D_param['units']['output']['speed']['Value'] * eval(self.D_param['units']['output']['speed']['Units'])
        self.speed_input_to_internal=(self.units_speed_input/self.units_speed_internal).si.value
        self.speed_internal_to_output=(self.units_speed_internal/self.units_speed_output).si.value
        # Conversion of temperatures for I/O
        self.units_temperature_input = self.D_param['units']['input']['temperature']['Value'] * eval(self.D_param['units']['input']['temperature']['Units'])
        self.units_temperature_internal = self.D_param['units']['internal']['temperature']['Value'] * eval(self.D_param['units']['internal']['temperature']['Units'])
        self.units_temperature_output = self.D_param['units']['output']['temperature']['Value'] * eval(self.D_param['units']['output']['temperature']['Units'])
        self.temperature_input_to_internal=(self.units_temperature_input/self.units_temperature_internal).si.value
        self.temperature_internal_to_output=(self.units_temperature_internal/self.units_temperature_output).si.value
        
        # Dimensionless versions of constants used in the code
        # Gravitational constant: units L^3/T^2M
        self.c_G=(c.G/self.units_length_internal**3*self.units_time_internal**2*self.units_mass_internal).si.value
        # To convert 3-d rms_speed to virial temperature T=mumH*sigma^2/k_B=mumH*sigma_3D^2/3k_B
        self.c_rms_speed_to_temperature = (
            self.D_param['astrophysics']['mumH']['Value'] * eval(self.D_param['astrophysics']['mumH']['Units']) *
            self.units_speed_internal**2 / (3. * c.k_B * self.units_temperature_internal ) ).si.value

    def __str__(self):
        """ 
        Simple method to print out parameters.  
        Returns 
        """
        print('Inbuilt options:')
        for item in self.D_option:
            print("{:20s}: {}".format(item,self.D_option[item]))
        print('\nRuntime options:')
        for item in self.D_param:
            print("{:20s}: {}".format(item,self.D_param[item]))
        print('\n')
        return ''
    
