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
      n_gal (int) : running total of the number of galaxies created by py-galaxies when processing the graph
      n_halo (int) : total number of halos in the graph.
      snap_n_halo (int[n_snap]) : number of halos in each snapshots.
      snap_first_halo_gid (int[n_snap]) : the first halo in each snapshot.
      n_sub (int) : total number of subhalos in the graph.
      snap_n_sub (int[n_snap]) : number of subhalos in each snapshots.
      snap_first_sub_gid (int[n_snap]) : = first subhalo in each snapshot.
      # root_mass (float) : mass of the root halo in the graph.

      # Halo properties:
      halo_desc_contribution (float[halo_n_desc]) : list of all particle contributions to descendant halos
      halo_desc_IDs_gid (int[halo_n_desc]) : list of all descendant halos
      halo_first_desc_gid (int[n_halo]) : location in graph of first descendant halo of each halo.
      half_mass_radius (float[n_halo]) : radius containing half of the halo particles
      mass  (float[n_halo]) : mass of each halo.
      mean_pos (float[n_halo,3]) : position of each halo.
      mean_vel (float[n_halo,3]) : velocity of each halo.
      halo_n_desc (float[n_halo]) : number of descendants of each halo.
      halo_n_sub (float[n_halo]) : number of subhalos for each halo
      halo_rms_radius (float[n_halo]) : rms radius of halo particles
      halo_rms_speed (float[n_halo]) : rms speed of halo particles
      halo_first_sub_gid (float[n_halo]) : the first subhalo in each halo

      # Subhalo properties:
      sub_desc_contribution (float[n_desc_sub]) : list of all particle contributions to descendant halos
      sub_desc_IDs_gid  (int[n_desc_sub]) : list of all descendant subhalos
      sub_first_desc_gid (int[n_sub]): location in graph of first descendant halo of each halo.
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
      parameters : obj : 'C_parameters class'
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
      graph = open_graph_file['graph_'+str(graph_ID)]
      # Attributes of graph
      for key,value in graph.attrs.items():
         exec('self.'+key+'=value')
      # Properties of graph: halos per snaphot (generation)
      self.snap_n_halo = graph['snap_n_halo'][:]           # Number of halos in each snapshot
      self.snap_first_halo_gid = graph['snap_first_halo'][:]    # First halo in each snapshot
      # Halo properties, fixed length arrays
      self.n_desc = graph['halo_n_desc'][:]
      self.halo_first_desc_gid = graph['halo_first_desc'][:]
      self.halo_mean_pos = graph['halo_pos'][:] * parameters.length_input_to_internal
      self.halo_mean_vel = graph['halo_vel'][:] * parameters.speed_input_to_internal
      self.halo_mass = graph['halo_mass'][:] * parameters.mass_input_to_internal
      self.halo_rms_radius = graph['halo_rms_radius'][:] * parameters.length_input_to_internal
      self.halo_rms_speed = graph['halo_rms_speed'][:] * parameters.speed_input_to_internal
      self.halo_half_mass_radius = graph['halo_half_mass_radius'][:] * parameters.length_input_to_internal
      # Halo properties, variable length arrays (because of possible graph branching)
      self.halo_desc_contribution = graph['halo_desc_contribution'][:]
      self.halo_desc_IDs_gid = graph['halo_desc_halo'][:]

      # Subhalos
      self.snap_n_sub = graph['snap_n_sub'][:]               # Number of subhalos in each snapshot
      self.snap_first_sub_gid = graph['snap_first_sub'][:]   # First subhalo in each snapshot
      self.halo_n_sub = graph['halo_n_sub'][:]               # Number of subhalos in each halo
      self.halo_first_sub_gid = graph['halo_first_sub'][:]   # First subhalo in each halo
      self.sub_first_desc_gid = graph['sub_first_desc'][:]
      self.sub_desc_contribution = graph['sub_desc_contribution'][:]
      self.sub_desc_IDs_gid = graph['sub_desc_sub'][:]
      # In case host halo not already set (it is tautologous)
      try:
         self.sub_host_gid = graph['sub_host'][:]
      except:
         self.sub_host_gid = np.full(self.n_sub,parameters.NO_DATA_INT,int)
         for i_halo in range(self.n_halo):
            for i_sub in self.halo_first_sub_gid[i_halo]+range(self.halo_n_sub[i_halo]):
               self.sub_host_gid[i_sub]=i_halo
      self.sub_pos = graph['sub_pos'][:] * parameters.length_input_to_internal
      self.sub_vel = graph['sub_vel'][:]* parameters.speed_input_to_internal
      self.sub_n_desc = graph['sub_n_desc'][:]
      self.sub_mass = graph['sub_mass'][:] * parameters.mass_input_to_internal
      self.sub_rms_speed = graph['sub_rms_speed'][:] * parameters.speed_input_to_internal
      self.sub_half_mass_radius = graph['sub_half_mass_radius'][:] * parameters.length_input_to_internal
      # Galaxies
      self.n_gal=0

