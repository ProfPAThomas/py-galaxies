class C_Halo:
    """
    A container for the properties needed for each halo.
   
    No sophisticated methods, it just truncates the GraphProperites class to 
    ensure data from the current generation is selected.
    
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
           The halo ID of currently being processed.
       graph : an instance of the class C_graph
           The graph contianing this halo.
       parameters : an instance of the class C_parameters
           The global parameters for this SAM run.
           
       """
       part_mass=parameters.part_mass
       self.graph_ID = graph_ID
       self.snap_ID = snap_ID
       self.halo_ID = halo_ID
       self.catalog_ID = graph.halo_catalog_halo_ids[halo_ID]
       self.mass = graph.nparts[halo_ID]*part_mass
       
       self.nprog = graph.nprog[halo_ID]
       self.prog_start = graph.prog_start_index[halo_ID]
       if b_HALO_FULL:
          self.prog_ids = graph.direct_prog_ids[self.prog_start:self.prog_start+self.nprog]
          self.prog_mass =  graph.direct_prog_contribution[self.prog_start:self.prog_start+self.nprog]*part_mass
           
       self.ndesc = graph.ndesc[halo_ID]
       self.desc_start = graph.desc_start_index[halo_ID]
       self.desc_end = self.desc_start + self.ndesc
       if b_HALO_FILL:
          self.desc_ids = graph.direct_desc_ids[self.desc_start:self.desc_start+self.ndesc]
       self.desc_mass = graph.direct_desc_contribution[self.desc_start:self.desc_start+self.ndesc]*part_mass
    
       self.mass_baryon = 0.
       self.mass_from_progenitors = 0.
       self.mass_baryon_from_progenitors = 0.
       self.inclusive_contribution = 0.
       self.b_done = False
       
      try:
         if parameters.model_switches['HOD']==True:
            self.mass_stars = 0.
            self.mass_stars_from_progenitors = 0.
      except:
         pass
