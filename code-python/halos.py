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
    halo_gid : int
        The halo ID within the graph of the halo currently being processed.
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
   def __init__(self,graph_ID,snap_ID,halo_gid,graph,parameters):
      """
       Clipping graph properties to the correct generation for halo use.
    
       Parameters
       ----------
       graph_ID : str
           The graph_ID (from HDF5 group).
       snap_ID : int
           The snapshot ID currently being processed.
       halo_gid : int
           The halo ID currently being processed, relative to the graph.
       graph : an instance of the class C_graph
           The graph containing this halo.
       parameters : an instance of the class C_parameters
           The global parameters for this SAM run.
           
      """
      # Read in halo properties from graph instance.  These should already be in internal code units.
      self.graph_ID = graph_ID
      self.snap_ID = snap_ID
      self.halo_gid = halo_gid
      self.halo_sid = self.halo_gid - graph.halo_start_gid[snap_ID]
      self.n_desc = graph.n_desc[halo_gid]
      self.desc_start_gid = graph.desc_start_gid[halo_gid]
      self.desc_end_gid = self.desc_start_gid + self.n_desc
      self.mass = graph.mass[halo_gid]
      self.pos = graph.mean_pos[halo_gid]
      self.vel = graph.mean_vel[halo_gid]
      self.half_mass_radius = graph.half_mass_radius[halo_gid]
      self.rms_radius = graph.rms_radius[halo_gid]
      self.rms_speed = graph.rms_speed[halo_gid]
      # Derived quantities
      # Using v^2=GM/r but for half mass
      self.half_mass_virial_speed = (0.5*parameters.c_G*self.mass/self.half_mass_radius)**(0.5)
      self.temperature = self.half_mass_virial_speed**2 * parameters.c_half_mass_virial_speed_to_temperature
      self.tau_dyn = 2.*self.half_mass_radius/self.half_mass_virial_speed
      # The following are properties of the SAM
      self.desc_main_sid = parameters.NO_DATA_INT  # Main descendant location in halos_this_snap
      self.mass_baryon = 0.
      self.mass_baryon_simple = 0.
      self.mass_from_progenitors = 0.
      self.mass_baryon_from_progenitors = 0.
      self.mass_gas_hot = 0.
      self.mass_metals_gas_hot = 0.
      self.mass_stars = 0.
      self.mass_metals_stars = 0.
      if parameters.b_HOD == True:
         self.star_formation_rate = 0.
      self.inclusive_contribution = 0.       
      self.n_dt = 0 # Number of times that this halo has been processed
      self.b_done = False # Has this halo been fully processed or not.
      self.b_desc_exists = True  # Needed to prevent trying to pass properties on to non-existent halo!

      # Subhalos
      self.n_sub = graph.n_sub_halo[halo_gid]
      self.sub_start_gid=parameters.NO_DATA_INT  # Updated below if n_sub>0
      # Copy in only those properties that we will use within the halo class
      if self.n_sub>0:
         sub_offset = graph.sub_start_gid[snap_ID]
         # Many of the following could be looked up as required but useful to define them here for quick reference
         self.sub_start_gid = graph.sub_start_halo_gid[halo_gid]
         self.sub_start_sid = self.sub_start_gid-sub_offset
         self.sub_end_gid = self.sub_start_gid+self.n_sub
         self.sub_end_sid = self.sub_start_sid+self.n_sub
         self.sub_mass = graph.sub_mass[self.sub_start_gid:self.sub_end_gid]
         self.sub_rel_pos = graph.sub_pos[self.sub_start_gid:self.sub_end_gid]-self.pos  # Assumes not already relative from MEGA
         self.sub_rel_vel = graph.sub_vel[self.sub_start_gid:self.sub_end_gid]-self.vel  #                 --"--

      # Galaxies
      self.n_gal = 0  # Total number of galaxies in halo + subhalos
      #self.gal_start_gid = parameters.NO_DATA_INT
      self.n_orphan = 0  # Galaxies not associated with a subhalo
      self.orphan_start_sid = parameters.NO_DATA_INT

      # Identify central subhalo.
      # For now, we will assume that there IS a central subhalo; later we may relax this assumption.
      # The appropriate metric for defining distance from the centre of phase-space also needs exploring.
      if self.n_sub>0:
         if parameters.b_lgalaxies:
            # In L-galaxies mode the most massive subhalo is assigned as the central subhalo.
            metric = self.sub_mass
            self.sub_central_gid = self.sub_start_gid+np.argmax(metric)
         else:
            # For now set equal to the minimum displacement from the phase-space ellipsoid.
            metric2 = np.sum((self.sub_rel_pos/self.rms_radius)**2,1)+np.sum((self.sub_rel_vel/self.rms_speed)**2,1)
            self.sub_central_gid = self.sub_start_gid+np.argmin(metric2)
         self.sub_central_sid = self.sub_central_gid - sub_offset
      else:
         self.sub_central_gid = parameters.NO_DATA_INT
         self.sub_central_sid = parameters.NO_DATA_INT

   def __str__(self):
      print('graph_ID =',self.graph_ID,',',end=' ')
      print('snap_ID =',self.snap_ID,',',end=' ')
      print('halo_gid =',self.halo_gid,flush=True)
      return ''

   def orphan_count(self,n_orphan):
      """
      Returns the current orphan galaxy location counter and then updates it.
      """
      orphan_next_sid = self.orphan_next_sid
      self.orphan_next_sid += n_orphan
      return orphan_next_sid

   def gal_loc(self,gal_start0_sid,gal_start_sid):
      """
      Sets the location of this halo's subhalo and orphan galaxies in the galaxy lookup table for this snapshot.
      """
      self.gal_start_sid = gal_start0_sid
      self.orphan_start_sid = gal_start_sid
      self.orphan_next_sid = self.orphan_start_sid # Will be used to keep track of orphans during update_halo phase
      return gal_start_sid+self.n_orphan
    
   def set_mass_baryon(self,subs,gals):
     """
     Calculates the total baryonic mass of the subhalo, including galaxies
     """
     self.mass_baryon = self.mass_gas_hot + self.mass_stars
     for i_sub in range(self.n_sub): 
        self.mass_baryon += subs[self.sub_start_sid+i_sub].mass_baryon
     # The orphan galaxies are not included in the subhalo baryon count, so add them in here
     if self.n_orphan >0:
        self.mass_baryon += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_gas_cold'])
        self.mass_baryon += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_stars_bulge'])
        self.mass_baryon += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_stars_disc'])
     return None


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
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.D_param['performance']['n_HDF5_io_rec']['Value']
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_gid',np.int32))
      dtype.append(('pos',np.float32,(3,)))
      dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass',np.float32))
      dtype.append(('temperature',np.float32))
      dtype.append(('rms_speed',np.float32))
      dtype.append(('half_mass_virial_speed',np.float32))
      dtype.append(('mass_baryon',np.float32))
      dtype.append(('mass_baryon_simple',np.float32))
      dtype.append(('mass_gas_hot',np.float32))
      dtype.append(('mass_metals_gas_hot',np.float32))
      dtype.append(('mass_stars',np.float32))
      dtype.append(('mass_metals_stars',np.float32))
      if parameters.b_HOD==True:
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
         self.io_buffer[self.i_rec]['halo_gid'] = halo.halo_gid
         self.io_buffer[self.i_rec]['pos'] = halo.pos * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['vel'] = halo.vel * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass'] = halo.mass * parameters.mass_internal_to_output      
         self.io_buffer[self.i_rec]['temperature'] = halo.temperature * parameters.temperature_internal_to_output
         self.io_buffer[self.i_rec]['rms_speed'] = halo.rms_speed * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['half_mass_virial_speed'] = halo.half_mass_virial_speed * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass_baryon']= halo.mass_baryon  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_baryon_simple']= halo.mass_baryon_simple  * parameters.mass_internal_to_output         
         self.io_buffer[self.i_rec]['mass_gas_hot'] = halo.mass_gas_hot  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_hot'] = halo.mass_metals_gas_hot  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_stars'] = halo.mass_stars  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars'] = halo.mass_metals_stars  * parameters.mass_internal_to_output
         if parameters.b_HOD==True:
            self.io_buffer[self.i_rec]['star_formation_rate'] = halo.star_formation_rate
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
