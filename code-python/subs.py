import h5py
import numpy as np

class C_sub:
   """
    A container for the properties needed for each subhalo.
   
    No sophisticated methods, it just truncates the GraphProperites class to 
    ensure data from the current generation is selected.

    Arrays are reference with slicing here [:], so will either be:
    * read from disk into memory, if graph contains dataset pointers;
    * a view into the corresponding graph arrays if those are already in memory.
    
    Attributes
    ----------
    graph_ID : str
        The graph_ID (from HDF5 group).
    snap_ID : int
         The snapshot ID currently being processed.
    halo_ID : int
        The halo ID of currently being processed.
    catalog_ID : int
        The ID of the halo corresponding to the original catalog.
    mass : int [Inconsistent with description below.]
        Mass of halo. Amount of dark matter particles * mass of particle.
    nprog : int
        The number of direct progenitors.
    prog_start_gid : int
        The index at which this halo's progenitors start.
    if b_HALO_FULL=True:
        prog_ids : narray of type 'int'
           Numpy array of the progenitor IDs for the halo.
        prog_mass : ndarry of type 'float'
           Numpy array of the progenitor mass contributions.
    ndesc : int
        The number of direct descendents.
    desc_start_gid : int
        The index at which this halo's descendents start.
    desc_end_gid : int
        The index at which this halo's descendents end.
    if b_HALO_FULL=True:
        desc_ids : ndarray of type 'int'
           Numpy array of the halo's descendent IDs of type.
        desc_mass : ndarray of type 'float'
           Numpy array of the descendent mass contributions.
    mass_baryon : float
        Mass of Baryons within the halo.
    mass_from_progenitors : float 
        Total mass of all the progenitor halos.
    mass_baryon_from_progenitors : float
        Total mass of all the Baryons contained within the progenitor halos.
    inclusive_contribution : float
        The amount of mass that goes 'missing' when a halo descends.
    b_done : bool
        Whether or not the halo has been processed.
 
    """   
   def __init__(self,graph_ID,snap_ID,sub_ID,graph,parameters):
      """
       Clipping graph properties to the correct generation for halo use.
    
       Parameters
       ----------
       graph_ID : str
           The graph_ID (from HDF5 group).
       snap_ID : int
           The snapshot ID currently being processed.
       sub_ID : int
           The sub ID currently being processed.
       graph : an instance of the class C_graph
           The graph contianing this halo.
       parameters : an instance of the class C_parameters
           The global parameters for this SAM run.
           
      """
      self.graph_ID = graph_ID
      self.snap_ID = snap_ID
      self.sub_ID = sub_ID
      # The following could be looked up as required but useful to define them here for quick reference
      self.host = graph.sub_host[sub_ID]
      self.n_desc = graph.sub_n_desc[sub_ID]
      self.desc_start_gid = graph.sub_desc_start_gid[sub_ID]
      self.desc_end_gid = self.desc_start_gid+self.n_desc
      self.desc_host_sid = parameters.NO_DATA_INT # descendant of host halo
      self.desc_main_sid = parameters.NO_DATA_INT

      # Intrinsic properties
      part_mass=parameters.part_mass
      self.mass = graph.sub_n_part[sub_ID]*part_mass
      self.pos = graph.sub_pos[sub_ID]
      self.vel = graph.sub_vel[sub_ID]
      # SAM properties
      self.ICM_mass = 0.
      self.hot_gas_mass = 0.
      # May be more than 1 galaxy in a subhalo (if progenitor subhalos merge):
      self.n_gal = 0
      self.gal_start_gid = parameters.NO_DATA_INT

      self.b_done = False

   def __str__(self):
      print('graph_ID =',self.graph_ID,',',end=' ')
      print('snap_ID =',self.snap_ID,',',end=' ')
      print('sub_ID =',self.sub_ID,flush=True)
      return ''

   def gal_count(self,n_gal):
      """
      Returns the current galaxy counter and updates it
      """
      gal_next = self.gal_next
      self.gal_next += n_gal
      return gal_next

   def gal_loc(self,gal_start):
      """
      Sets the location of this subhalo's galaxies in the galaxy lookup table.
      """
      # Require subhalo to have at least 1 galaxy
      self.n_gal = max(self.n_gal,1)
      self.gal_start = gal_start
      self.gal_end = self.gal_start + self.n_gal
      self.gal_next = self.gal_start # Will be used to keep track of galaxies during update_halo phase.
      return self.gal_end

class C_sub_output:
   
   """
   This class contains the attributes and methods for the subhalo output files.
   Attributes
   ----------

   Methods
   -------
   __init__
   append - add subhalos to output buffer
   close - flush io buffer then close HDF5 file
   flush - flush output buffer to HDF5 dataset
   
   """
   def __init__(self,parameters):
      """
      Opens the subhalo output file.
      Creates the subhalo output buffer.
      Creates the HDF5 halo dataset.

      Parameters:
      -----------
      parameters : obj : C_parameters
         Contains the global run paramters.
      """
      # Open file for output
      self.sub_file = h5py.File(parameters.subhalo_output_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.D_param['performance']['n_HDF5_io_rec']['Value']
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_ID',np.int32))
      dtype.append(('sub_ID',np.int32))
      dtype.append(('pos',np.float32,(3,)))
      dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass',np.float32))
      dtype.append(('ICM_mass',np.float32))
      dtype.append(('hot_gas_mass',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.sub_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and
      closes the halo output file
      """
      self.flush()
      # self.dataset.close() # There does not seem to be a need to close datasets.
      self.sub_file.close()
      return None

   def flush(self):
      """
      Writes io buffer to the HDF5 dataset and resets.
      """
      self.dataset.resize((self.dataset.shape[0]+self.i_rec,))
      self.dataset[-self.i_rec:]=self.io_buffer[:self.i_rec]
      self.i_rec=0
      return None

   def append(self,subs,parameters):
      """
      Extracts the quantities desired for halo_output and adds them to the io buffer,
      flushing if required.
      Parameters
      ----------
         halos - list of C_halo objects to be output
         parameters - C_parameters class file containing the global run parameters
      """
      for sub in subs:
         self.io_buffer[self.i_rec]['graph_ID'] = sub.graph_ID
         self.io_buffer[self.i_rec]['snap_ID'] = sub.snap_ID
         self.io_buffer[self.i_rec]['halo_ID'] = sub.sub_ID
         self.io_buffer[self.i_rec]['pos'] = sub.pos
         self.io_buffer[self.i_rec]['vel'] = sub.vel
         self.io_buffer[self.i_rec]['mass'] = sub.mass
         self.io_buffer[self.i_rec]['ICM_mass']= sub.ICM_mass
         self.io_buffer[self.i_rec]['hot_gas_mass']= sub.hot_gas_mass
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
