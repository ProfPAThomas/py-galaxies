import numpy as np

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

        import yaml
        
        self.param_file = param_file
        self.D_param = yaml.load(open(param_file),Loader=yaml.Loader)

        # extract key variables for ease of later use
        self.graph_input_file = self.D_param['input_files']['graph_file']
        self.snap_input_file = self.D_param['input_files']['snap_file']
        self.halo_output_file = self.D_param['output_files']['halo_file']
        self.galaxy_output_file = self.D_param['output_files']['galaxy_file']
        self.omega_m = self.D_param['cosmology']['omega_m']['Value']
        self.baryon_fraction = self.D_param['cosmology']['baryon_fraction']['Value']
        self.n_HDF5_io_rec = self.D_param['performance']['n_HDF5_io_rec']['Value']
        #   self.sub_halo = self.D_param['model_switches']['sub_halo']['Value']

        # Loop over model switches, creating a boolen flag for each
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

        # Create the dtype that we will need to store galaxy properties.
        # Not strictly a parameter but makes sense to do it here.
        D_dtype_gal=[
            ('halo_graph',np.int32),
            ('halo_snap',np.int32),
            ('sub_graph',np.int32),
            ('sub_snap',np.int32),
            ('first_prog_snap',np.int32),
            ('next_prog_snap',np.int32),
            ('stellar_mass',np.float32),
            ('cold_gas_mass',np.float32)
            ]
        # Properties that are optional depending upon run-tiem flags:
        if self.b_HOD:
            D_dtype_gal.append(('HOD_mass',np.float32))
        self.dtype_gal=np.dtype(D_dtype_gal)

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
    
