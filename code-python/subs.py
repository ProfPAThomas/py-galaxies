"""
Class files for storage and output of subhalo properties
"""

import ctypes
import h5py
import numpy as np

# First define the dtype for subhalo properties in a way that is compatible with ctypes
D_sub=np.dtype([
   ('mass',ctypes.c_double),
   ('pos',ctypes.c_double*3),
   ('vel',ctypes.c_double*3),
   ('rms_speed',ctypes.c_double),
   ('half_mass_radius',ctypes.c_double),
   ('half_mass_virial_speed',ctypes.c_double),
   # Derived properties
   ('temperature',ctypes.c_double),
   ('tau_dyn',ctypes.c_double),
   # SAM properties
   ('mass_baryon',ctypes.c_double),
   ('mass_gas_hot',ctypes.c_double),
   ('mass_metals_gas_hot',ctypes.c_double),
   ('mass_stars',ctypes.c_double),
   ('mass_metals_stars',ctypes.c_double)],
   align=True)

class C_sub:
   """
    A container for the properties needed for each subhalo.
   
    No sophisticated methods, it just extracted the relevant properties from the graph.

    All physical quantities are in internal code units.
    
    Attributes
    ----------
    graph_ID : int
        The graph_ID within the input file.
    snap_ID : int
         The snapshot ID currently being processed.
    halo_gid : int
        The host halo location within the graph.
    halo_sid : int
        The host halo location within the snapshot.
    sub_gid : int
        The subhalo location within the graph.
    sub_sid : int
        The subhalo location within the snapshot.
    b_done : bool
        Whether or not the subhalo has been fully processed.
    desc_end_gid : int
        The index relative to the graph at which this subhalo's descendents end.
    desc_halo_sid : int
        The index relative to the next snapshot of the descendant of the host halo of this subhalo.
    desc_start_gid : int
        The index relative to the graph at which this subhalo's descendents start.
    gal_central_sid : int
        The location in the current galaxy array of the most massive galaxy in the subhalo.
    gal_end_sid : int
        The location in the current galaxy array of the last galaxy in the subhalo (+1 because of python indexing)
    gal_next_sid : int
        Galaxy counter used when updating galaxies.
    gal_start_sid : int
        The location in the current galaxy array of the first galaxy in the subhalo
    half_mass_radius : float
        The radius containing half the total mass in the DM-only sim
    half_mass_virial_speed : float
        The circular speed at the half-mass radius
    mass : float
        The DM-only mass of the subhalo.
    mass_baryon : float
        Mass of baryons within the subhalo, inclusive of galaxies
    mass_gas_hot : float
        The mass of hot gas
    mass_metals_gas_hot : float
        The mass of metals in hot gas
    mass_metals_stars : float
        The mass of metals in stars
    mass_stars : float
        The mass of stars
    n_desc : int
        The number of direct descendants.
    n_dt : int
        Number of times that this halo has been processed this snapshot
    n_gal : int
        The number of galaxies in the subhalo.
    pos : float[3]
        The position of the subhalo.
    rms_speed : float
        The rms speed of the subhalo in the DM-only sim.
    tau_dyn : float
        The dynamical time at twice the half-mass radius.
    temperature : float
        The virial temperature as derived from the virial speed.
    vel : float[3]
        The velocity of the subhalo.
 
    """
   # Define template for new subhalo instances
   template=np.empty(1,dtype=np.dtype(D_sub,align=True))
   template['mass'] = 0.
   template['pos'] = [0.,0.,0.]
   template['vel'] = [0.,0.,0.]
   template['rms_speed'] = 0.
   template['half_mass_radius'] = 0.
   template['half_mass_virial_speed'] = 0.
   template['temperature'] = 0.
   template['tau_dyn'] = 0.
   template['mass_baryon']=0.
   template['mass_gas_hot']=0.
   template['mass_metals_gas_hot']=0.
   template['mass_stars']=0.
   template['mass_metals_stars']=0.

   def __init__(self,graph_ID,snap_ID,sub_gid,graph,parameters):
      """
       Read in the subhalo properties from the graph, including ranges for descendants;
       initialise galaxy counters.
    
       Parameters
       ----------
       graph_ID : int
           The graph_ID within the input data file.
       snap_ID : int
           The snapshot ID currently being processed.
       sub_gid : int
           The subhalo ID relative to the graph.
       graph : obj : C_graph
           The graph containing this halo.
       parameters : obj : C_parameters
           The global parameters for this SAM run.
           
      """
      self.graph_ID = graph_ID
      self.snap_ID = snap_ID
      self.halo_gid = graph.sub_host_gid[sub_gid]
      self.halo_sid = self.halo_gid - graph.snap_first_halo_gid[snap_ID]
      self.sub_gid = sub_gid
      self.sub_sid = self.sub_gid - graph.snap_first_sub_gid[snap_ID]
      self.n_desc = graph.sub_n_desc[sub_gid]
      self.desc_start_gid = graph.sub_first_desc_gid[sub_gid]
      self.desc_end_gid = self.desc_start_gid+self.n_desc
      self.desc_halo_sid = parameters.NO_DATA_INT

      # Intrinsic properties: to be passed to C routines, so store as a numpy array.
      # Create and initialise
      self.props=C_sub.template.copy() # Both the C_sub and the .copy() are required here.
      self.props['mass'] = graph.sub_mass[sub_gid]
      self.props['pos'] = graph.sub_pos[sub_gid]
      self.props['vel'] = graph.sub_vel[sub_gid]
      self.props['rms_speed'] = graph.sub_rms_speed[sub_gid]
      self.props['half_mass_radius'] = graph.sub_half_mass_radius[sub_gid]
      self.props['half_mass_virial_speed'] = (0.5*parameters.c_G*self.props['mass']/self.props['half_mass_radius'])**(0.5)
      self.props['temperature'] = self.props['half_mass_virial_speed']**2 * parameters.c_half_mass_virial_speed_to_temperature
      self.props['tau_dyn'] = 2.*self.props['half_mass_radius']/self.props['half_mass_virial_speed']

      # May be more than 1 galaxy in a subhalo (if progenitor subhalos merge):
      self.n_gal = 0
      self.gal_start_sid = parameters.NO_DATA_INT
      self.gal_end_sid = parameters.NO_DATA_INT
      self.gal_central_sid = parameters.NO_DATA_INT
      self.b_done = False # Has this subhalo been fully processed or not
      self.n_dt = 0 # Number of times that this subhalo has been processed

   def __str__(self):
      print('graph_ID =',self.graph_ID,',',end=' ')
      print('snap_ID =',self.snap_ID,',',end=' ')
      print('sub_gid =',self.sub_gid,flush=True)
      return ''

   def gal_count(self,n_gal):
      """
      Returns the current galaxy counter and updates it.

      Parameters
      ----------
      n_gal : int
         The number of galaxies being added to the subhalo

      Returns
      -------
      int
         The galaxy counter upon entry (ie before updating)
      """
      gal_next_sid = self.gal_next_sid
      self.gal_next_sid += n_gal
      return gal_next_sid

   def gal_loc(self,gal_start_sid):
      """
      Sets the location of this subhalo's galaxies in the galaxy lookup table.
      Also sets the location counter gal_next_sid to be the equal to the gal_star_sid, to initialise updating.

      Parameters
      ----------
      gal_start_sid : int
         The start location of this subhalos galaxies in the galaxy array for this graph/snapshot.

      Returns
      -------
      int
         The start location of this subhalos galaxies in the galaxy array for this graph/snapshot (+1 for python indexing).
      """
      # Require subhalo to have at least 1 galaxy
      self.n_gal = max(self.n_gal,1)
      self.gal_start_sid = gal_start_sid
      self.gal_end_sid = self.gal_start_sid + self.n_gal
      self.gal_next_sid = self.gal_start_sid # Will be used to keep track of galaxies during update_halo phase.
      return self.gal_end_sid

   def sum_mass_baryon(self,gals):
       """
       Calculates the total baryonic mass of the subhalo, including galaxies.
       Returns value rather than setting it because being used as a check on simpler method.

       Parameters
       ----------
       gals : obj : D_gal[]
          Array of records for the galaxies in this subhalo.

       Returns
       -------
       float
          The baryonic mass of the subhalo, inclusive of galaxies.
       """
       mass_baryon = self.props['mass_gas_hot'] + self.props['mass_stars'] + \
                       np.sum(gals[self.gal_start_sid:self.gal_end_sid]['mass_baryon'])
       return mass_baryon
        
class C_sub_output:
   
   """
   This class contains the attributes and methods for the subhalo output files.

   Attributes
   ----------
   sub_file : obj : File
      HDF5 file for subhalo output
   i_rec : int
      Counter for how many records have been created.
   io_buffer : obj " D_gal[n_rec]
      Storage for subhalo records prior to outputting. 
   n_rec : int
      Number records to be buffered before outputting.
   dataset : obj : HDF5 dataset
      HDF5 dataset to which the data is output. 
   """
   def __init__(self,parameters):
      """
      Opens the subhalo output file.
      Creates the subhalo output buffer.
      Creates the HDF5 halo dataset.

      Parameters
      ----------
      parameters : obj : C_parameters
         Contains the global run parameters.
      """
      # Open file for output
      self.sub_file = h5py.File(parameters.subhalo_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.n_HDF5_io_rec
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_gid',np.int32))
      dtype.append(('sub_gid',np.int32))
      dtype.append(('pos',np.float32,(3,)))
      dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass',np.float32))
      dtype.append(('rms_speed',np.float32))
      dtype.append(('half_mass_virial_speed',np.float32))
      dtype.append(('temperature',np.float32))
      dtype.append(('mass_gas_hot',np.float32))
      dtype.append(('mass_metals_gas_hot',np.float32))
      dtype.append(('mass_stars',np.float32))
      dtype.append(('mass_metals_stars',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.sub_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and closes the subhalo output file.
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
      Extracts the quantities desired for halo_output and adds them to the io buffer, flushing if required.

      Parameters
      ----------
         subs : obj : C_sub[]
            list of C_sub objects to be output.
         parameters : obj : C_parameters
            The global run parameters.
      """
      for sub in subs:
         self.io_buffer[self.i_rec]['graph_ID'] = sub.graph_ID
         self.io_buffer[self.i_rec]['snap_ID'] = sub.snap_ID
         self.io_buffer[self.i_rec]['halo_gid'] = sub.halo_gid
         self.io_buffer[self.i_rec]['sub_gid'] = sub.sub_gid
         self.io_buffer[self.i_rec]['pos'] = sub.props['pos'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['vel'] = sub.props['vel'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass'] = sub.props['mass'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['rms_speed'] = sub.props['rms_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['half_mass_virial_speed'] = sub.props['half_mass_virial_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['temperature'] = sub.props['temperature'] * parameters.temperature_internal_to_output
         self.io_buffer[self.i_rec]['mass_gas_hot']= sub.props['mass_gas_hot'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_hot']= sub.props['mass_metals_gas_hot'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_stars']= sub.props['mass_stars'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars']= sub.props['mass_metals_stars'] * parameters.mass_internal_to_output
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
