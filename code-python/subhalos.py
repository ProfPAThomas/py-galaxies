import h5py
import numpy as np

class C_subhalo:
   """
    A container for the properties needed for each subhalo.
   
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
   def __init__(self,graph_ID,snap_ID,subhalo_ID,graph,parameters):
      """
       Clipping graph properties to the correct generation for halo use.
    
       Parameters
       ----------
       graph_ID : str
           The graph_ID (from HDF5 group).
       snap_ID : int
           The snapshot ID currently being processed.
       subhalo_ID : int
           The subhalo ID currently being processed.
       graph : an instance of the class C_graph
           The graph contianing this halo.
       parameters : an instance of the class C_parameters
           The global parameters for this SAM run.
           
      """
      self.graph_ID = graph_ID
      self.snap_ID = snap_ID
      self.subhalo_ID = subhalo_ID
      # The following could be looked up as required but useful to define them here for quick reference
      self.host = graph.sub_host[subhalo_ID]
      self.n_desc = graph.sub_n_desc[subhalo_ID]
      self.desc_start = graph.sub_desc_start[subhalo_ID]
      self.desc_host_snap = parameters.NO_DATA_INT # descendant of host halo
      self.desc_main_snap = parameters.NO_DATA_INT

      # Intrinsic properties
      part_mass=parameters.part_mass
      self.mass = graph.sub_n_part[subhalo_ID]*part_mass
      self.pos = graph.sub_pos[subhalo_ID]
      self.vel = graph.sub_vel[subhalo_ID]
      # SAM properties
      self.ICM_mass = 0.
      self.hot_gas_mass = 0.
      # May be more than 1 galaxy in a subhalo (if progenitor subhalos merge):
      self.n_gal = 0
      self.gal_start = parameters.NO_DATA_INT

      self.b_done = False

   def init_gal(self,gal_start):
      """
      Sets the location of this subhalo's galaxies in the galaxy lookup table.
      """
      # Require subhalo to have at least 1 galaxy
      self.n_gal = max(self.n_gal,1)
      self.gal_start = gal_start
      self.gal_next = self.gal_start # WIll be used to keep track of galaxies during update_halo phase.
      gal_start = gal_start + self.n_gal
      return None

   def add_gal(self,n_gal):
      """
      Returns the current galaxy counter and updates it.
      """
      gal_next = self.gal_next
      self.gal_next += n_gal
      return gal_next
