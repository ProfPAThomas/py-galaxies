class C_parameters:
    """Read in yml parameters and store them.
    
    Simple class to read in and store parameters from the yml file. Simple
    methods included to print out the parameters etc.
    
    Attributes
    ----------
    param_file : str
        Filepath to the yml file containing the model parameters.
    verbosity : int
        The level of detail needed for debugging messages.
    b_debug : bool
        To print out debugging messages, true or false. 
    D_param : dictionary 
        Dictionary containing contents of yml file.
    graph_input_file : str
        The filepath to the input graph HDF5 file.
    galaxy_output_file : str
        The filepath to the halo output HDF5 file.
    halo_output_file : str
        The filepath to the galaxy output HDF5 file.
    baryon_fraction : float
        Cosmic baryon fraction.
    n_HDF5_io_rec : int
        IO HDF5 buffer size.
    omega_m : float
        Density parameter.
    
    methods:
        __init__
        __str__ 
    
    """
    
    def __init__(self,param_file,verbosity,b_debug):
        """ 
        Key parameters for the model
        
        Parameters
        ----------
        param_file : str
           Filepath to the yml file containing the model parameters.
        verbosity : int
           The level of detail needed for debugging messages.
        b_debug : bool
           To print out debugging messages, true or false. 
        """

        import yaml
        
        self.param_file = param_file
        self.verbosity = verbosity
        self.b_debug = b_debug
        self.D_param = yaml.load(open(param_file),Loader=yaml.Loader)

        # extract key variables for ease of later use
        self.graph_input_file = self.D_param['input_files']['graph_file']
        self.halo_output_file = self.D_param['output_files']['halo_file']
        self.galaxy_output_file = self.D_param['output_files']['galaxy_file']
        self.omega_m = self.D_param['cosmology']['omega_m']['Value']
        self.baryon_fraction = self.D_param['cosmology']['baryon_fraction']['Value']
        self.n_HDF5_io_rec = self.D_param['performance']['n_HDF5_io_rec']['Value']
        #   self.sub_halo = self.D_param['model_switches']['sub_halo']['Value']


    def __str__(self):
        """ 
        Simple method to print out parameters.
            
        Returns 
        -------
        None
        """

        for item in self.D_param:
            print("{:20s}: {}".format(item,self.D_param[item]))
        return '\n'
