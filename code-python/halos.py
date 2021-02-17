import h5py
import numpy as np

class C_halo:
   """
    A container for the properties needed for each halo.
   
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
    prog_start : int
        The index at which this halo's progenitors start.
    if b_HALO_FULL=True:
        prog_ids : narray of type 'int'
           Numpy array of the progenitor IDs for the halo.
        prog_mass : ndarry of type 'float'
           Numpy array of the progenitor mass contributions.
    ndesc : int
        The number of direct descendents.
    desc_start : int
        The index at which this halo's descendents start.
    desc_end : int
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
   def __init__(self,graph_ID,snap_ID,halo_ID,graph,parameters):
      """
       Clipping graph properties to the correct generation for halo use.
    
       Parameters
       ----------
       graph_ID : str
           The graph_ID (from HDF5 group).
       snap_ID : int
           The snapshot ID currently being processed.
       halo_ID : int
           The halo ID currently being processed.
       graph : an instance of the class C_graph
           The graph contianing this halo.
       parameters : an instance of the class C_parameters
           The global parameters for this SAM run.
           
      """
      self.graph_ID = graph_ID
      self.snap_ID = snap_ID
      self.halo_ID = halo_ID
      # The following could be looked up as required but useful to define them here for quick reference
      # Note that we define links relative to the snap to enable lookup in the halo/subhalo instance lists
      halo_offset = graph.halo_start[snap_ID]
      #self.n_prog = graph.n_prog[halo_ID]
      #self.prog_start = graph.prog_start[halo_ID]
      #self.prog_start_snap = self.prog_start - halo_offset
      #self.prog_end = self.prog_start + self.n_prog
      #self.prog_end_snap = self.prog_start_snap + self.n_prog
      self.n_desc = graph.n_desc[halo_ID]
      self.desc_start = graph.desc_start[halo_ID]
      self.desc_start_snap = self.desc_start - halo_offset
      self.desc_end = self.desc_start + self.n_desc
      self.desc_end_snap = self.desc_start_snap + self.n_desc
      self.desc_main_snap = parameters.NO_DATA_INT # To be populated later
      part_mass=parameters.part_mass
      self.mass = graph.n_part[halo_ID]*part_mass
      self.pos = graph.mean_pos[halo_ID]
      self.vel = graph.mean_vel[halo_ID]
      self.rms_radius = graph.rms_radius[halo_ID]
      self.rms_speed = graph.rms_speed[halo_ID]
      # The following are properties of the SAM
      self.desc_main_snap = parameters.NO_DATA_INT  # Main descendant location in halos_this_snap
      self.mass_baryon = 0.
      self.mass_from_progenitors = 0.
      self.mass_baryon_from_progenitors = 0.
      if parameters.b_HOD==True:
         self.mass_stars = 0.
         self.mass_stars_from_progenitors = 0.
         self.star_formation_rate = 0.
      self.inclusive_contribution = 0.       
      self.b_done = False

      # Subhalos
      self.n_sub = graph.n_sub_halo[halo_ID]
      self.sub_start=parameters.NO_DATA_INT  # Updated below if n_sub>0
      # Copy in only those properties that we will use within the halo class
      if self.n_sub>0:
         sub_offset = graph.sub_start[snap_ID]
         # Many of the following could be looked up as required but useful to define them here for quick reference
         self.sub_start = graph.sub_start_halo[halo_ID]
         self.sub_start_snap = self.sub_start - sub_offset
         self.sub_end = self.sub_start+self.n_sub
         self.sub_end_snap = self.sub_end - sub_offset
         self.sub_mass = graph.sub_n_part[self.sub_start:self.sub_end]*part_mass
         self.sub_rel_pos = graph.sub_pos[self.sub_start:self.sub_end]-self.pos  # Assumes not already relative from MEGA
         self.sub_rel_vel = graph.sub_vel[self.sub_start:self.sub_end]-self.vel  #                 --"--

      # Galaxies
      self.n_gal = 0  # Total number of galaxies in halo + subhalos
      self.gal_start = parameters.NO_DATA_INT
      self.n_orphan = 0  # Galaxies not associated with a subhalo
      self.orphan_start = parameters.NO_DATA_INT

      # Identify central subhalo.
      # For now, we will assume that there IS a central subhalo; later we may relax this assumption
      if self.n_sub>0:
         metric2 = np.sum((self.sub_rel_pos/self.rms_radius)**2,1)+np.sum((self.sub_rel_vel/self.rms_speed)**2,1)
         self.sub_central = self.sub_start+np.argmin(metric2)
      else:
         self.central_sub = parameters.NO_DATA_INT

   def __str__(self):
      print('graph_ID =',self.graph_ID,',',end=' ')
      print('snap_ID =',self.snap_ID,',',end=' ')
      print('halo_ID =',self.halo_ID,flush=True)
      return ''

   def init_gal(self,gal_start0,gal_start):
      """
      Sets the location of this halo's subhalo and orphan galaxies in the galaxy lookup table.
      """
      self.gal_start = gal_start0
      self.orphan_start = gal_start
      self.orphan_next = self.orphan_start # Will be used to keep track of orphans during update_halo phase
      gal_start = gal_start + self.n_orphan
      return None

   def add_gal(self,n_orphan):
      """
      Returns the current orphan galaxy counter and updates it.
      """
      orphan_next = self.orphan_next
      self.orphan_next += n_orphan
      return orphan_next

class C_halo_output:
   
   """
   This class contains the attributes and methods for the halo output files.
   Attributes
   ----------

   Methods
   -------
   __init__
   append - add halos to output buffer
   close - flush io buffer then close HDF5 file
   flush - flush output buffer to HDF5 dataset
   
   """
   def __init__(self,parameters):
      """
      Opens the halo output file.
      Creates the halo output buffer.
      Creates the HDF5 halo dataset.

      Parameters:
      -----------
      parameters : obj : C_parameters
         Contains the gloabal run paramters.
      """
      # Open file for output
      self.halo_file = h5py.File(parameters.halo_output_file,'w')
      # Counter for total number of halos written
      self.i_halo = 0
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.D_param['performance']['n_HDF5_io_rec']['Value']
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_ID',np.int32))
      dtype.append(('mass',np.float32))
      dtype.append(('mass_baryon',np.float32))
      if parameters.b_HOD==True:
         dtype.append(('mass_stars',np.float32))
         dtype.append(('star_formation_rate',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.halo_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and
      closes the halo output file
      """
      self.flush()
      # self.dataset.close() # There does not seem to be a need to close datasets.
      self.halo_file.close()
      return None

   def flush(self):
      """
      Writes io buffer to the HDF5 dataset and resets.
      """
      self.dataset.resize((self.dataset.shape[0]+self.i_rec,))
      self.dataset[-self.i_rec:]=self.io_buffer[:self.i_rec]
      self.i_rec=0
      return None

   def append(self,halos,parameters):
      """
      Extracts the quantities desired for halo_output and adds them to the io buffer,
      flushing if required.
      Parameters
      ----------
         halos - list of C_halo objects to be output
         parameters - C_parameters class file containing the global run parameters
      """
      for halo in halos:
         self.io_buffer[self.i_rec]['graph_ID'] = halo.graph_ID
         self.io_buffer[self.i_rec]['snap_ID'] = halo.snap_ID
         self.io_buffer[self.i_rec]['halo_ID'] = halo.halo_ID
         self.io_buffer[self.i_rec]['mass'] = halo.mass
         self.io_buffer[self.i_rec]['mass_baryon']= halo.mass_baryon
         if parameters.b_HOD==True:
            self.io_buffer[self.i_rec]['mass_stars'] = halo.mass_stars
            self.io_buffer[self.i_rec]['star_formation_rate'] = halo.star_formation_rate
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
