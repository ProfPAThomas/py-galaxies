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
      self.nprog = graph.nprog[halo_ID]
      self.prog_start = graph.prog_start_index[halo_ID]
      self.ndesc = graph.ndesc[halo_ID]
      self.desc_start = graph.desc_start_index[halo_ID]
    
      part_mass=parameters.part_mass
      self.mass = graph.nparts[halo_ID]*part_mass
      self.mass_baryon = 0.
      self.mass_from_progenitors = 0.
      self.mass_baryon_from_progenitors = 0.
      if parameters.b_HOD==True:
         self.mass_stars = 0.
         self.mass_stars_from_progenitors = 0.
         self.star_formation_rate = 0.
      self.inclusive_contribution = 0.       
      self.b_done = False


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
