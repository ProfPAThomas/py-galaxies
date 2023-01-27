import numpy as np
    
class C_graph:
   """
   A container for all data contained within a single graph.
    
   Data attributes:
   ----------------
   graph_ID : int
      The ID of this graph
   n_halo : int
      Number of halos in the graph
   root_mass : float
      Mass of the root halo
   # Properties of graph: halos per snaphot (generation)
   generation_id : ndarray : int32
      ID of this snapshot (NO_DATA if no halos); redundant
   generation_length : ndarray : float32
      Number of halos in each generation(snapshot)
   generation_start_index : ndarry : float32
      First generation (snaphot) with halos in it
   # Halo properties, fixed length arrays
   desc_start_index : ndarray : int32
      First entry in descendant arrays
   half_mass_radius : ndarray : float32
      Half mass radius of halo (Units?)
   half_mass_speed : ndarray : float32
      Half mass speed of halo (Units?)
   halo_catalog_halo_ids : ndarray : int64
      IDs in original halo catalogue from N-body simulation
   mean_pos : ndarray : float32[3]
      Position of CofM of halo (Units? Presumably comoving?)
   mean_vel : ndarray : float32[3]
      Velocity of CofM of halo (Units?)
   ndesc : ndarray : int32
      Number of descendant halos
   nparts : ndarray : int32
      Number of particles in this halo
   nprog : ndarray : int32
      Number of progenitors of this halo
   prog_start_index : ndarray : int32
      First entry in progenitor arrays
   redshifts : ndarray : float32
      Redshift of halo (redundant)
   rms_radius : ndarray : float32
      RMS size of halo (Units?)
   rms_speed : ndarray : float32
      RMS speed of halo particles (3-D velocity dispersion) (Units?)
   snapshots : ndarray : int32
      Snapshot of halo (redundant)
   v_max : ndarray : float32
      Maximum rotation speed of halo (Units?)      
   # Halo properties, variable length arrays (because of possible graph branching)
   direct_desc_contribution : ndarray : int32
      Number of particles contributed to descendant
   direct_desc_ids : ndarry : int32
      The IDs of the direct descendants
   direct_prog_contribution : ndarray : int32
      Number of particles contributed by progenitor
   direct_prog_ids : ndarry : int32
      The IDs of the direct progenitors
   ... and more for subhalos ...

    Authors: Andrew Bowell & Peter Thomas
        
   """
   
   def __init__(self,graph_ID,open_graph_file,parameters):
      """ 
      Opening HDF datasets for a graph and saving them to the class.
        
      Parameters
      ----------
      graph_ID : int
         The ID of the graph to be opened.
      open_graph_file : obj : 'File'
         Open HDF5 file that is to be read from.
      parameters : obj : 'Class'
         parameters class object.  
         Included in case graph file has parameters to be read.
      
      """

      """ Note:
      * Unclear whether better to read in array data here (the [:] forces that), or
      at point of use.  The latter is more memory efficient but may be slower.
      * Many variables are renamed upon read to simplify and impose uniform convention throughout the rest of the code.
      * Uncomment variables as required.  Conditional inclusion depending upon run-time parameter choices is fine.
      * Indices are relative to the graph (with 0 in the past) unless labelled otherwise.
      * Similarly, counts, etc refer to the whole graph unless labelled otherwise.
      """
      self.graph_ID = graph_ID
      graph = open_graph_file[str(graph_ID)]
      # Attributes of graph
      self.n_halo = graph.attrs['nhalos_in_graph']            # Total number of halos in the graph
      try:
         self.n_sub = graph.attrs['sub_nhalos_in_graph']      # Total number of subhalos in the graph
      except:
         self.n_sub = 0
      self.root_mass=graph.attrs['root_mass']
      # Properties of graph: halos per snaphot (generation)
      self.snap_ID = graph['generation_id'][:] # Contains no data flag if there are no halos.
      # Note that the following is set to the no data flag, not 0, if there are no halos - need to correct
      self.n_halo_snap = graph['generation_length'][:]        # Number of halos in each snapshot
      self.n_halo_snap=np.where(self.n_halo_snap == parameters.NO_DATA_INT, 0, self.n_halo_snap)
      self.halo_start_gid = graph['generation_start_index'][:]    # First halo in each snapshot
      # Halo properties, fixed length arrays
      self.desc_start_gid = graph['desc_start_index'][:]
      #self.half_mass_radius = graph['half_mass_radius'][:]
      #self.half_mass_speed = graph['half_mass_velocity_radius'][:]
      #self.catalog_halo_ids = graph['halo_catalog_halo_ids'][:]
      self.mean_pos = graph['mean_pos'][:]
      self.mean_vel = graph['mean_vel'][:]
      self.n_desc = graph['ndesc'][:]
      self.n_part = graph['nparts'][:]
      self.n_prog = graph['nprog'][:]
      self.prog_start_gid = graph['prog_start_index'][:]
      #self.redshifts = graph['redshifts'][:]
      self.rms_radius = graph['rms_radius'][:]
      self.rms_speed = graph['3D_velocity_dispersion'][:]
      #self.snapshots = graph['snapshots'][:]
      #self.v_max = graph['v_max'][:]
      # Halo properties, variable length arrays (because of possible graph branching)
      self.desc_contribution = graph['direct_desc_contribution'][:]
      self.desc_IDs_gid = graph['direct_desc_ids'][:]
      #self.prog_contribution = graph['direct_prog_contribution'][:]
      #self.prog_IDs = graph['direct_prog_ids'][:]      
      # Subhalos
      if self.n_sub == 0:
         self.n_sub_halo = np.zeros(self.n_halo,dtype=np.int32)
      else:
         self.n_sub_halo = graph['nsubhalos'][:]               # Number of subhalos in each halo
         #self.sub_rms_speed = graph['sub_3D_velocity_dispersion'][:]
         #self.sub_catalog_halo_ids = graph['subhalo_catalog_halo_ids'][:]
         self.sub_desc_start_gid = graph['sub_desc_start_index'][:]
         self.sub_desc_contribution = graph['sub_direct_desc_contribution'][:]
         self.sub_desc_IDs_gid = graph['sub_direct_desc_ids'][:]
         #self.sub_direct_prog_contribution = graph['sub_direct_prog_contribution'][:]
         #self.sub_direct_prog_ids = graph['sub_direct_prog_ids'][:]
         #self.sub_generation_id = graph['sub_generation_id'][:]
         self.n_sub_snap = graph['sub_generation_length'][:]       # Number of subhalos in this snapshot
         self.n_sub_snap=np.where(self.n_sub_snap == parameters.NO_DATA_INT, 0, self.n_sub_snap)
         self.sub_start_gid = graph['sub_generation_start_index'][:]   # First subhalo in each snapshot
         #self.sub_half_mass_radius = graph['sub_half_mass_radius'][:]
         #self.sub_half_mass_speed = graph['sub_half_mass_velocity_radius'][:]
         self.sub_host_gid = graph['host_halos'][:]
         self.sub_pos = graph['sub_mean_pos'][:]
         self.sub_vel = graph['sub_mean_vel'][:]
         self.sub_n_desc = graph['sub_ndesc'][:]
         self.sub_n_part = graph['sub_nparts'][:]
         #self.sub_n_prog = graph['sub_nprog'][:]
         #self.sub_prog_start_gid = graph['sub_prog_start_index'][:]
         #self.sub_redshifts = graph['sub_redshifts'][:]
         #self.sub_rms_radius = graph['sub_rms_radius'][:]
         #self.sub_snapshots = graph['sub_snapshots'][:]
         self.sub_start_halo_gid = graph['subhalo_start_index'][:]      # First subhalo in each halo
         #self.sub_v_max = graph['sub_v_max']
      # Galaxies
      self.n_gal=0
