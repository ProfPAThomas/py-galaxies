"""
Class files for storage and output of halo properties
"""

import ctypes
import h5py
import numpy as np

# First define the dtype for halo properties in a way that is compatible with ctypes
D_halo=np.dtype([
   ('mass',ctypes.c_double),
   ('mass_from_progenitors',ctypes.c_double),
   ('pos',ctypes.c_double*3),
   ('vel',ctypes.c_double*3),
   ('rms_speed',ctypes.c_double),
   ('half_mass_radius',ctypes.c_double),
   ('half_mass_virial_speed',ctypes.c_double),
   ('rms_radius',ctypes.c_double),
   # Derived properties
   ('temperature',ctypes.c_double),
   ('tau_dyn',ctypes.c_double),
   # SAM properties
   ('mass_baryon',ctypes.c_double),
   ('mass_baryon_from_progenitors',ctypes.c_double),
   ('mass_gas_eject',ctypes.c_double),
   ('mass_metals_gas_eject',ctypes.c_double),
   ('mass_gas_hot',ctypes.c_double),
   ('mass_metals_gas_hot',ctypes.c_double),
   ('mass_stars',ctypes.c_double),
   ('mass_metals_stars',ctypes.c_double)],
   align=True)

class C_halo:
   """
   A container for the properties needed for each halo.
  
   No sophisticated methods, it just truncates the GraphProperites class to 
   ensure data from the current generation is selected.

   All physical quantities are in internal code units.
    
   Attributes
   ----------
   graph_ID : str
       The graph_ID (from HDF5 group).
   snap_ID : int
       The snapshot ID currently being processed.
   halo_gid : int
       The halo location within the graph of the halo currently being processed.
   halo_sid : int
       The halo location within the snapshot of the halo currently being processed.
   b_done : bool
       Whether or not the halo has been fully processed.
   delta_baryon_dthalo : float
       The accretion needed per halo timestep to bring the baryon content up to the universal mean.
   desc_end_gid : int
       The index at which this halo's descendents end (+1 as is usual for python).
   desc_main_sid : int
       The main descendant location in this snapshot (ie in halos_this_snap).
   desc_start_gid : int
       The index at which this halo's descendents start.
   half_mass_radius : float
       The radius containing half the total mass of the halo in the DM-only sim.
   half_mass_virial_speed : float
       The circular speed at the half-mass radius.
   mass : float
       The DM-only mass of this halo.
   mass_baryon : float
       Mass of baryons within the halo, inclusive of subhalos and galaxies.
   mass_baryon_from_progenitors : float
       Total mass of all the baryons contained within the progenitor halos.
   mass_from_progenitors : float 
       Total DM-only mass of all the progenitor halos.
   mass_gas_eject : float
       The mass of ejected gas
   mass_gas_hot : float
       The mass of hot gas in the halo, exclusive of subhalos.
   mass_stars : float
       The mass of stars in the halo, exclusive of subhalos.
   mass_metals_gas_eject : float
       The mass of metals in ejected gas.
   mass_metals_gas_hot : float
       The mass of metals in hot gas in the halo, exclusive of subhalos.
   mass_metals_stars : float
       The mass of metals in stars in the halo, exclusive of subhalos.
   n_desc : int
       The number of direct descendants of this halo
   n_dt : int
       Number of times that this halo has been processed this snapshot
   n_gal : int
       The number of galaxies in the halo, inclusive of subhalos
   n_orphan : int
       The number of orphan galaxies (i.e. galaxies not contained in subhalos)
   n_sub : int
       The number of subhalos.
   orphan_start_sid : int
       The location of the first orphan galaxy within the current snapshot galaxy array
   orphan_next_sid : int
       To track the processing of orphan galaxies during the halo_update phase.
   pos : float[3]
       The position of the halo.
   rms_radius : float
       The rms radius of the halo particles in the DM-only sim.
   rms_speed : float
       The rms speed of the halo particles in the DM-only sim.
   sub_central_gid : int
       The location in the graph of the subhalo at the centre of the halo (if any)
   sub_central_sid : int
       The location in the snapshot of the subhalo at the centre of the halo (if any)
   sub_end_gid : int
       The location in the graph of the last subhalo (+1 because of python indexing)
   sub_end_sid : int
       The location in the snapshot of the last subhalo (+1 because of python indexing)
   sub_mass : float[n_sub]
       The DM-only masses of the subhalos
   sub_rel_pos : float[n_sub,3]
       The positions of the subhalos relative to that of the halo
   sub_rel_vel : float[n_sub,3]
       The velocities of the subhalos relative to that of the halo
   sub_start_gid : int
       The location in the graph of the first subhalo
   sub_start_sid : int
       The location in the snapshot of the first subhalo
   tau_dyn : float
       Twice the dynamical time at the half-mass radius (= dynamical time at twice the half mass radius for isothermal sphere).
   temperature : float
       The temperature as derived from the virial speed.
   vel : float[3]
       Velocity of the halo
   """
   # Define template for new halo instances
   template=np.empty(1,dtype=np.dtype(D_halo,align=True))
   template['mass'] = 0.
   template['mass_from_progenitors'] = 0.
   template['pos'] = [0.,0.,0.]
   template['vel'] = [0.,0.,0.]
   template['rms_speed'] = 0.
   template['half_mass_radius'] = 0.
   template['half_mass_virial_speed'] = 0.
   template['rms_radius'] = 0.
   template['temperature'] = 0.
   template['tau_dyn'] = 0.
   template['mass_baryon']=0.
   template['mass_baryon_from_progenitors']=0.
   template['mass_gas_eject']=0.
   template['mass_metals_gas_eject']=0.
   template['mass_gas_hot']=0.
   template['mass_metals_gas_hot']=0.
   template['mass_stars']=0.
   template['mass_metals_stars']=0.
   
   def __init__(self,graph_ID,snap_ID,halo_gid,graph,parameters):
      """
      Read in the halos properties from the graph, including ranges for decendants, subhalos and galaxies.
    
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
      self.halo_sid = self.halo_gid - graph.snap_first_halo_gid[snap_ID]
      self.n_desc = graph.halo_n_desc[halo_gid]
      self.desc_start_gid = graph.halo_first_desc_gid[halo_gid]
      self.desc_end_gid = self.desc_start_gid + self.n_desc
      self.desc_main_sid = parameters.NO_DATA_INT  # Main descendant location in halos_this_snap
      
      # Intrinsic properties: to be passed to C routines, so store as a numpy array.
      # Create and initialise
      self.props=C_halo.template.copy() # Both the C_halo and the .copy() are required here.
      self.props['mass'] = graph.halo_mass[halo_gid]
      self.props['pos'] = graph.halo_mean_pos[halo_gid]
      self.props['vel'] = graph.halo_mean_vel[halo_gid]
      self.props['half_mass_radius'] = graph.halo_half_mass_radius[halo_gid]
      self.props['rms_radius'] = graph.halo_rms_radius[halo_gid]
      self.props['rms_speed'] = graph.halo_rms_speed[halo_gid]
      # Derived quantities
      # Using v^2=GM/r but for half mass
      self.props['half_mass_virial_speed'] = (0.5*parameters.c_G*self.props['mass']/self.props['half_mass_radius'])**(0.5)
      self.props['temperature'] = self.props['half_mass_virial_speed']**2 * parameters.c_half_mass_virial_speed_to_temperature
      self.props['tau_dyn'] = 2.*self.props['half_mass_radius']/self.props['half_mass_virial_speed']
      
      self.n_dt = 0 # Number of times that this halo has been processed
      self.b_done = False # Has this halo been fully processed or not.

      # Subhalos
      self.n_sub = graph.halo_n_sub[halo_gid]
      self.sub_start_gid=parameters.NO_DATA_INT  # Updated below if n_sub>0
      # Copy in only those properties that we will use within the halo class
      if self.n_sub>0:
         sub_offset = graph.snap_first_sub_gid[snap_ID]
         # Many of the following could be looked up as required but useful to define them here for quick reference
         self.sub_start_gid = graph.halo_first_sub_gid[halo_gid]
         self.sub_start_sid = self.sub_start_gid-sub_offset
         self.sub_end_gid = self.sub_start_gid+self.n_sub
         self.sub_end_sid = self.sub_start_sid+self.n_sub
         self.sub_mass = graph.sub_mass[self.sub_start_gid:self.sub_end_gid]
         self.sub_rel_pos = graph.sub_pos[self.sub_start_gid:self.sub_end_gid]-self.props['pos']  # Assumes not already relative from MEGA
         self.sub_rel_vel = graph.sub_vel[self.sub_start_gid:self.sub_end_gid]-self.props['vel']  #                 --"--

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

   def accrete_primordial_gas(self,base_metallicity):
      """
      Updates the baryon content as calcuated previously in the halo set_mass_baryon method.
      Any excess baryons arrive in the form of base_metallicity hot gas.

      Arguments
      ---------
      base_metallicity : float
         The metallicity of pristine, infalling gas.

      Returns
      -------
      None
      """
      delta_mass=self.delta_baryon_dthalo
      self.props['mass_baryon']+=delta_mass
      self.props['mass_gas_hot']+=delta_mass
      self.props['mass_metals_gas_hot']+=delta_mass*base_metallicity

      return None

   def gal_loc(self,gal_start0_sid,gal_start_sid):
      """
      Sets the location of this halo's subhalo and orphan galaxies in the galaxy lookup table for this snapshot.

      Parameters
      ----------
      gal_start0_sid : int 
         The location of the first galaxy in the halo, inclusive of subhalos, in the galaxy array for this snapshot.
      gal_start_sid : int 
         The location of the first orphan galaxy in the halo, exlcusive of subhalos, in the galaxy array for this snapshot.

      Returns
      -------
      int
         The updated galaxy pointer (ie location of the last galaxy in this halo +1 for python indexing)
      """
      self.gal_start_sid = gal_start0_sid
      self.orphan_start_sid = gal_start_sid
      self.orphan_next_sid = self.orphan_start_sid # Will be used to keep track of orphans during update_halo phase
      return gal_start_sid+self.n_orphan
    
   def orphan_count(self,n_orphan):
      """
      Returns the current orphan galaxy location counter and then updates it.

      Parameters
      ----------
      n_orphan : int
         The number of orphans to add to the halo

      Returns
      -------
      int
         The current orphan count for this snapshot of the graph
      """
      orphan_next_sid = self.orphan_next_sid
      self.orphan_next_sid += n_orphan
      return orphan_next_sid

   def set_mass_baryon(self,subs,gals,baryon_fraction,n_dt_halo):
      """
      Calculates and updates the total baryonic mass of the subhalo, inclusive of subhalos and galaxies.
      Also determines the accretion required to bring the halo up to the universal mean baryon fraction, 
      or the sum of the baryon content from the progenitors, whichever is larger (so that baryons are not lost).

      Parameters
      ----------
      subs : obj: C_sub[n_sub]
        The subhalo instances contained within this halo
      gals : obj: D_gal[n_gal]
         The galaxies  contained within this halo, inclusive of subhalos
      baryon_fraction : float
         The universal baryon fraction.
      n_dt_halo : int
         The number of halo timesteps in this snapshot interval.

      Returns
      -------
      None
      """
      self.props['mass_baryon'] = self.props['mass_gas_hot'] + self.props['mass_gas_eject'] + self.props['mass_stars']
      for i_sub in range(self.n_sub): 
         self.props['mass_baryon'] += subs[self.sub_start_sid+i_sub].props['mass_baryon']
      # The orphan galaxies are not included in the subhalo baryon count, so add them in here
      if self.n_orphan >0:
         self.props['mass_baryon'] += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_gas_cold'])
         self.props['mass_baryon'] += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_stars_bulge'])
         self.props['mass_baryon'] += np.sum(gals[self.orphan_start_sid:self.orphan_start_sid+self.n_orphan]['mass_stars_disc'])
         
      # The amount of accretion needed per halo timestep to bring the baryon content up to the universal mean over a snapshot interval.
      self.delta_baryon_dthalo=max(0.,baryon_fraction*max(self.props['mass'],self.props['mass_from_progenitors'])-self.props['mass_baryon'])/float(n_dt_halo)

      return None


class C_halo_output:
   
   """
   This class contains the attributes and methods for the halo output files. 

   Attributes
   ----------
   halo_file : obj : File
      HDF5 file for halo output
   i_rec : int
      Counter for how many records have been created.
   io_buffer : obj : D_gal[n_rec]
      Storage for halo records prior to outputting. 
   n_rec : int
      Number records to be buffered before outputting.
   dataset : obj : HDF5 dataset
      HDF5 dataset to which the data is output.
   """
   def __init__(self,parameters):
      """
      Opens the halo output file.
      Creates the halo output buffer.
      Creates the HDF5 halo dataset.

      Parameters
      ----------
      parameters : obj : C_parameters
         Contains the global run parameters.
      Returns
      -------
      None
      """
      # Open file for output
      self.halo_file = h5py.File(parameters.halo_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.n_HDF5_io_rec
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
      dtype.append(('mass_gas_hot',np.float32))
      dtype.append(('mass_metals_gas_hot',np.float32))
      dtype.append(('mass_gas_eject',np.float32))
      dtype.append(('mass_metals_gas_eject',np.float32))
      dtype.append(('mass_stars',np.float32))
      dtype.append(('mass_metals_stars',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.halo_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and
      closes the halo output file.
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
         halos : obj : C_halo[]
            list of C_halo objects to be output.
         parameters : obj : C_parameters
            The global run parameters.
      """
      for halo in halos:
         self.io_buffer[self.i_rec]['graph_ID'] = halo.graph_ID
         self.io_buffer[self.i_rec]['snap_ID'] = halo.snap_ID
         self.io_buffer[self.i_rec]['halo_gid'] = halo.halo_gid
         self.io_buffer[self.i_rec]['pos'] = halo.props['pos'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['vel'] = halo.props['vel'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass'] = halo.props['mass'] * parameters.mass_internal_to_output      
         self.io_buffer[self.i_rec]['temperature'] = halo.props['temperature'] * parameters.temperature_internal_to_output
         self.io_buffer[self.i_rec]['rms_speed'] = halo.props['rms_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['half_mass_virial_speed'] = halo.props['half_mass_virial_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass_baryon']= halo.props['mass_baryon']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_gas_hot'] = halo.props['mass_gas_hot']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_hot'] = halo.props['mass_metals_gas_hot']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_gas_eject'] = halo.props['mass_gas_eject']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_eject'] = halo.props['mass_metals_gas_eject']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_stars'] = halo.props['mass_stars']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars'] = halo.props['mass_metals_stars']  * parameters.mass_internal_to_output
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
