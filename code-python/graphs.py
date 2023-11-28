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
      graph_ID (int) : The unique identifier of the graph within the file.
      n_gal (int) : Running total of the number of galaxies created by py-galaxies when processing the graph
      snap_n_halo (int[n_snap]) : The number of halos in each snapshot.
      snap_first_halo_gid (int[n_snap]) : The first halo in each snapshot, relative to the graph.
      snap_n_sub (int[n_snap]) : The number of subhalos in each snapshot.
      snap_first_sub_gid (int[n_snap]) : The first subhalo in each snapshot, relative to the graph.

      # Halo properties:
      halo_desc_contribution (float[halo_n_desc]) : The particle contributions to descendant halos.
      halo_desc_IDs_gid (int[halo_n_desc]) : The locations in the graph of all descendant halos of each halo.
      halo_first_desc_gid (int[n_halo]) : The locations in the graph of the first descendant halo of each halo.
      halo_first_sub_gid (float[n_halo]) : The locations in the graph of first subhalo in each halo.
      half_mass_radius (float[n_halo]) : The radii containing half of the halo particles.
      halo_mass  (float[n_halo]) : The mass of each halo.
      halo_mean_pos (float[n_halo,3]) : The position of each halo.
      halo_mean_vel (float[n_halo,3]) : The velocity of each halo.
      halo_n_desc (float[n_halo]) : The number of descendants of each halo.
      halo_n_sub (float[n_halo]) : The number of subhalos of each halo
      halo_rms_radius (float[n_halo]) : The rms radii of halo particles
      halo_rms_speed (float[n_halo]) : The rms speeds of halo particles

      # Subhalo properties:
      sub_desc_contribution (float[n_desc_sub]) : The particle contributions to descendant subhalos.
      sub_desc_IDs_gid  (int[n_desc_sub]) : The locations in the graph of all descendant subhalos.
      sub_first_desc_gid (int[n_sub]): The locations in the graph of the first descendant subhalo of each subhalo.
      sub_half_mass_radius (float[n_sub]) : The radii containing half of the subhalo particles.
      sub_host_gid (int[n_sub]) : The host halo of each subhalo.
      sub_mass  (float[n_sub]) : The mass of each subhalo.
      sub_n_desc (float[n_sub]) : The number of descendants of each subhalo.
      sub_pos (float[n_sub,3]) : The position of each subhalo.
      sub_rms_speed (float[n_sub]) : The rms speed of subhalo particles within each subhalo.
      sub_vel (float[n_sub,3]) : The velocity of each subhalo.        
   """
   
   def __init__(self,graph_ID,open_graph_file,parameters):
      """ 
      Opening HDF datasets for a graph and saving them to the class.
        
      Parameters
      ----------
      graph_ID : int
         The ID within the file of the graph to be opened.
      open_graph_file : obj : 'File'
         Open HDF5 file that is to be read from.
      parameters : obj : 'C_parameters class'
         Contains parameters describing the simulation.
      """

      """ Note:
      * Unclear whether better to read in array data here (the [:] forces that), or
      at point of use.  The latter is more memory efficient but may be slower.
      * Conditional inclusion of variables depending upon run-time parameter choices is fine.
      * Indices are relative to the graph (with 0 in the past) as denoted by the _gid postfix.
      """
      self.graph_ID = graph_ID
      graph = open_graph_file['graph_'+str(graph_ID)]

      # Attributes of graph
      for key,value in graph.attrs.items():
         exec('self.'+key+'=value')

      # Properties of graph: halos and subhalos per snaphot
      self.snap_n_halo = graph['snap_n_halo'][:]
      self.snap_first_halo_gid = graph['snap_first_halo'][:]
      self.snap_n_sub = graph['snap_n_sub'][:]
      self.snap_first_sub_gid = graph['snap_first_sub'][:]

      # Halo properties, fixed length arrays
      self.halo_first_desc_gid = graph['halo_first_desc'][:]
      self.halo_first_sub_gid = graph['halo_first_sub'][:]
      self.halo_half_mass_radius = graph['halo_half_mass_radius'][:] * parameters.length_input_to_internal
      self.halo_mean_pos = graph['halo_pos'][:] * parameters.length_input_to_internal
      self.halo_mean_vel = graph['halo_vel'][:] * parameters.speed_input_to_internal
      self.halo_mass = graph['halo_mass'][:] * parameters.mass_input_to_internal
      self.halo_n_desc = graph['halo_n_desc'][:]
      self.halo_n_sub = graph['halo_n_sub'][:]
      self.halo_rms_radius = graph['halo_rms_radius'][:] * parameters.length_input_to_internal
      self.halo_rms_speed = graph['halo_rms_speed'][:] * parameters.speed_input_to_internal
      # Halo properties, variable length arrays (because of possible graph branching)
      self.halo_desc_contribution = graph['halo_desc_contribution'][:]
      self.halo_desc_IDs_gid = graph['halo_desc_halo'][:]

      # Subhalo properties
      self.sub_desc_contribution = graph['sub_desc_contribution'][:]
      self.sub_desc_IDs_gid = graph['sub_desc_sub'][:]
      self.sub_first_desc_gid = graph['sub_first_desc'][:]
      self.sub_half_mass_radius = graph['sub_half_mass_radius'][:] * parameters.length_input_to_internal
      # In case host halo not already set (it is tautologous)
      try:
         self.sub_host_gid = graph['sub_host'][:]
      except:
         self.sub_host_gid = np.full(self.n_sub,parameters.NO_DATA_INT,int)
         for i_halo in range(self.n_halo):
            for i_sub in self.halo_first_sub_gid[i_halo]+range(self.halo_n_sub[i_halo]):
               self.sub_host_gid[i_sub]=i_halo
      self.sub_mass = graph['sub_mass'][:] * parameters.mass_input_to_internal
      self.sub_n_desc = graph['sub_n_desc'][:]
      self.sub_pos = graph['sub_pos'][:] * parameters.length_input_to_internal
      self.sub_vel = graph['sub_vel'][:]* parameters.speed_input_to_internal
      self.sub_rms_speed = graph['sub_rms_speed'][:] * parameters.speed_input_to_internal
      # Galaxies
      self.n_gal=0

