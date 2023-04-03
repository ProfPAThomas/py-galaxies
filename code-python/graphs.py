"""
Class file for reading and storage of graph properties
"""

import numpy as np
    
class C_graph:
   """
   A container for all data contained within a single graph.

   Note that all attributes are in internal units.
   In the following attributes the sizes of the arrays are specified in terms of the total number of snapshots, halos, subhalos, progenitors and descendants within the graph.
   
   Attributes:
      # Graph and snapshot properties:
      graph_ID (int) : graph_ID.
      halo_start_gid (int[n_snap]) : the first halo in each snapshot.
      n_gal (int) : running total of the number of galaxies created by py-galaxies when processing the graph
      n_halo (int) : total number of halos in the graph.
      n_halo_snap (int[n_snap]) : number of halos in each snapshots.
      n_sub (int) : total number of subhalos in the graph.
      n_sub_snap (int[n_snap]) : number of subhalos in each snapshots.
      root_mass (float) : mass of the root halo in the graph.
      snap_ID (int[n_snap]) : list of snapshots with halos: contains no data flag if there are no halos.
      sub_start_gid (int[n_snap]) : = first subhalo in each snapshot.

      # Halo properties:
      desc_contribution (float[n_desc_halo]) : list of all particle contributions to descendant halos
      desc_IDs_gid (int[n_desc_halo]) : list of all descendant halos
      desc_start_gid (int[n_halo]) : location in graph of first descendant halo of each halo.
      half_mass_radius (float[n_halo]) : radius containing half of the halo particles
      mass  (float[n_halo]) : mass of each halo.
      mean_pos (float[n_halo,3]) : position of each halo.
      mean_vel (float[n_halo,3]) : velocity of each halo.
      n_desc (float[n_halo]) : number of descendants of each halo.
      n_prog  (float[n_halo]) : number of progenitors of each halo.
      n_sub_halo (float[n_halo]) : number of subhalos for each halo
      prog_start_gid (int) : location in graph of first progenitor halo of each halo.
      rms_radius (float[n_halo]) : rms radius of halo particles
      rms_speed (float[n_halo]) : rms speed of halo particles
      sub_start_halo_gid (float[n_halo]) : the first subhalo in each halo

      # Subhalo properties:
      sub_desc_contribution (float[n_desc_sub]) : list of all particle contributions to descendant halos
      sub_desc_IDs_gid  (int[n_desc_sub]) : list of all descendant subhalos
      sub_desc_start_gid (int[n_sub]): location in graph of first descendant halo of each halo.
      sub_half_mass_radius (float[n_sub]) : radius containing half of the subhalo particles
      sub_host_gid (int[n_sub]) : the host halo of each subhalo
      sub_mass  (float[n_sub]) : mass of each subhalo.
      sub_n_desc (float[n_sub]) : number of descendants of each subhalo.
      sub_pos (float[n_sub,3]) : position of each subhalo.
      sub_rms_speed (float[n_sub]) : rms speed of subhalo particles
      sub_vel (float[n_sub,3]) : velocity of each subhalo.        
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
         Contains parameters describing the simulation.
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
      part_mass_internal=parameters.part_mass*parameters.mass_input_to_internal
      # Attributes of graph
      self.n_halo = graph.attrs['nhalos_in_graph']            # Total number of halos in the graph
      try:
         self.n_sub = graph.attrs['sub_nhalos_in_graph']      # Total number of subhalos in the graph
      except:
         self.n_sub = 0
      self.root_mass=graph.attrs['root_mass'] * part_mass_internal 
      # Properties of graph: halos per snaphot (generation)
      self.snap_ID = graph['generation_id'][:] # Contains no data flag if there are no halos.
      # Note that the following is set to the no data flag, not 0, if there are no halos - need to correct
      self.n_halo_snap = graph['generation_length'][:]        # Number of halos in each snapshot
      self.n_halo_snap=np.where(self.n_halo_snap == parameters.NO_DATA_INT, 0, self.n_halo_snap)
      self.halo_start_gid = graph['generation_start_index'][:]    # First halo in each snapshot
      # Halo properties, fixed length arrays
      self.desc_start_gid = graph['desc_start_index'][:]
      self.mean_pos = graph['mean_pos'][:] * parameters.length_input_to_internal
      self.mean_vel = graph['mean_vel'][:] * parameters.speed_input_to_internal
      self.n_desc = graph['ndesc'][:]
      self.mass = graph['nparts'][:] * part_mass_internal
      self.n_prog = graph['nprog'][:]
      self.prog_start_gid = graph['prog_start_index'][:]
      self.rms_radius = graph['rms_radius'][:] * parameters.length_input_to_internal
      self.rms_speed = graph['3D_velocity_dispersion'][:] * parameters.speed_input_to_internal
      self.half_mass_radius = graph['half_mass_radius'][:] * parameters.length_input_to_internal
      # Halo properties, variable length arrays (because of possible graph branching)
      self.desc_contribution = graph['direct_desc_contribution'][:]
      self.desc_IDs_gid = graph['direct_desc_ids'][:]

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
         self.sub_host_gid = graph['host_halos'][:]
         self.sub_pos = graph['sub_mean_pos'][:] * parameters.length_input_to_internal
         self.sub_vel = graph['sub_mean_vel'][:]* parameters.speed_input_to_internal
         self.sub_n_desc = graph['sub_ndesc'][:]
         self.sub_mass = graph['sub_nparts'][:] * part_mass_internal
         self.sub_rms_speed = graph['sub_3D_velocity_dispersion'][:] * parameters.speed_input_to_internal
         self.sub_half_mass_radius = graph['sub_half_mass_radius'][:] * parameters.length_input_to_internal
         self.sub_start_halo_gid = graph['subhalo_start_index'][:]      # First subhalo in each halo
         #self.sub_v_max = graph['sub_v_max']
      # Galaxies
      self.n_gal=0
