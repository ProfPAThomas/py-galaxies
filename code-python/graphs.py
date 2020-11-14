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
      Unclear whether better to read in array data here (the [:] forces that), or
      at point of use.  The latter is more memory efficient but may be slower. """
      
      self.graph_ID = graph_ID
      graph = open_graph_file[str(graph_ID)]
      # Attributes of graph
      self.n_halo=graph.attrs['nhalos_in_graph']
      self.root_mass=graph.attrs['root_mass']
      # Properties of graph: halos per snaphot (generation)
      self.generation_id = graph['generation_id'][:]
      self.generation_length = graph['generation_length'][:]
      self.generation_start_index = graph['generation_start_index'][:]
      # Halo properties, fixed length arrays
      self.desc_start_index = graph['desc_start_index'][:]
      self.half_mass_radius = graph['half_mass_radius'][:]
      self.half_mass_speed = graph['half_mass_velocity_radius'][:]
      self.halo_catalog_halo_ids = graph['halo_catalog_halo_ids'][:]
      self.mean_pos = graph['mean_pos'][:]
      self.mean_vel = graph['mean_vel'][:]
      self.ndesc = graph['ndesc'][:]
      self.nparts = graph['nparts'][:]
      self.nprog = graph['nprog'][:]
      self.prog_start_index = graph['prog_start_index'][:]
      self.redshifts = graph['redshifts'][:]
      self.rms_radius = graph['rms_radius'][:]
      self.rms_speed = graph['3D_velocity_dispersion'][:]
      self.snapshots = graph['snapshots'][:]
      self.v_max = graph['v_max'][:]
      # Halo properties, variable length arrays (because of possible graph branching)
      self.direct_desc_contribution = graph['direct_desc_contribution'][:]
      self.direct_desc_ids = graph['direct_desc_ids'][:]
      self.direct_prog_contribution = graph['direct_prog_contribution'][:]
      self.direct_prog_ids = graph['direct_prog_ids'][:]      
      # Subhalos (galaxies)
      try:
         self.nsubhalos = graph['nsubhalos'][:]
         self.sub_desc_start_index = graph['sub_desc_start_index'][:]
         self.sub_direct_desc_contribution = graph['sub_direct_desc_contribution'][:]
         self.sub_direct_desc_ids = graph['sub_direct_desc_ids'][:]
         self.sub_direct_prog_contribution = graph['sub_direct_prog_contribution'][:]
         self.sub_direct_prog_ids = graph['sub_direct_prog_ids'][:]
         self.sub_generation_id = graph['sub_generation_id'][:]
         self.sub_generation_length = graph['sub_generation_length'][:]
         self.sub_generation_start_index = graph['sub_generation_start_index'][:]
         self.sub_mean_pos = graph['sub_mean_pos'][:]
         self.sub_ndesc = graph['sub_ndesc'][:]
         self.sub_nparts = graph['sub_nparts'][:]
         self.sub_nprog = graph['sub_nprog'][:]
         self.sub_prog_start_index = graph['sub_prog_start_index'][:]
         self.sub_redshifts = graph['sub_redshifts'][:]
         self.sub_snapshots = graph['sub_snapshots'][:]
         self.subhalo_catalog_halo_ids = graph['subhalo_catalog_halo_ids'][:]
      except:
         self.nsubhalos=np.zeros(self.n_halo)
         
