import numpy as np
    
class C_graph:
   """A container for all data contained within a single graph
    
    This class consists of data gathered from a single graph. This documentat-
    tion will essentially be copied from Will's. 
    
    
    Attributes
    ----------
    desc_start_index : ndarray of type 'int'
        The starting index (pointer) for each halo’sentries in all descendant 
        halo arrays (i.e.direct_desc_ids,direct_desc_contribution,​ etc.). 
        Entries containing 2**30 have no descendants.
    direct_desc_contribution : ndarray of type 'int'
        The number of dark matter particles contributed ​to​ each direct
        descendent ​from the halo.
    direct_desc_ids : ndarray of type 'int'
        The descendent halo IDs, extracted usingdesc_start_index​ and ​ndesc​
    direct_prog_contribution : ndarray of type 'int'
        The number of dark matter particles contributed ​by​ each direct 
        progenitor ​to​ thehalo.
    direct_prog_ids : ndarray of type 'int'
        The progenitor halo IDs, extracted using prog_start_index​ and ​nprog​.
    generation_id : ndarray of type 'int'
        The ID (or number) associated with each generation. Counting starts
        from the earliest snapshot.
    generation_length : ndarray of type 'int'
        The number of halos in each generation
    generation_start_index : ndarray of type 'int'
        The starting index (pointer) for each host halo generation.
    halo_catalog_halo_ids : ndarray of type 'int'
        The halo catalog ID assigned to each halo.
    mean_pos : ndarray of type 'float'
        The mean position of the particles in the halo.
    ndesc : ndarray of type 'int'
        The number of descendants for each halo.
    nparts : ndarray of type 'int'
        The number of dark matter particles in each halo. 
    nprog : ndarray of type 'int'
        The number of progenitors for each halo.
    prog_start_index : ndarray of type 'int'
        The starting index (pointer) for each halo’s entries in all progenitor
        halo arrays (i.e. direct_prog_ids, direct_prog_contribution,​ etc.).
        Entries containing 2**30 have nodescendants.
    redshifts : ndarray of type 'float'
        The redshift for each halo in the graph.
    snapshots : ndarray of type 'int'
        The index of the snapshot in the snapshottext file dictated in the 
        param file, for each halo.
    sub_desc_start_index : ndarray of type 'int'
        The starting index (pointer) for each subhalo’s entries in all
        descendant subhalo arrays (i.e. ​sub_direct_desc_ids,
        sub_direct_desc_contribution,​ etc.). Entries containing 2**30 have 
        no descendants.
    sub_direct_desc_contribution : ndarray of type 'int'
        The number of dark matter particles contributed ​to​ each direct
        descendent ​from the subhalo.
    sub_direct_desc_ids : ndarray of type 'int'
        The descendent subhalo IDs, extracted using ​sub_desc_start_index​ 
        and sub_ndesc​.
    sub_direct_prog_contribution : ndarray of type 'int'
        The number of dark matter particles contributed ​by​ each direct
        progenitor ​to​ the subhalo.
    sub_direct_prog_ids : ndarray of type 'int'
        The progenitor subhalo IDs, extracted using sub_prog_start_index​ and 
        ​sub_nprog​.
    sub_generation_id : ndarray of type 'int'
        The ID (or number) associated with each generation. Counting starts
        from the earliest snapshot.
    sub_generation_length : ndarray of type 'int'
        The number of subhalos in each generation.
    sub_generation_start_index : indarray of type 'int'nt
        The starting index (pointer) for each subhalo generation.
    sub_mean_pos : ndarray of type 'float'
        The mean position of the particles in the subhalo.
    sub_ndesc : ndarray of type 'int'
        The number of descendents for eachsubhalo.
    sub_nparts : ndarray of type 'int'
        The number of dark matter particles in eachhalo.
    sub_nprog : ndarray of type 'int'
        The number of progenitors for each halo.
    sub_prog_start_index : ndarray of type 'int'
        The starting index (pointer) for eachsubhalo’s entries in all 
        progenitor subhaloarrays (i.e. ​sub_direct_prog_ids,
        sub_direct_prog_contribution,​ etc.).Entries containing 2**30 have
        no descendants.
    sub_redshifts : ndarray of type 'float'
        The redshift for each subhalo in the graph.
    sub_snapshots : ndarray of type 'int'
        The index of the snapshot in the snapshottext file dictated in the
        param file, for each subhalo.
    subhalo_catalog_halo_ids : ndarray of type 'int'
        The subhalo catalog ID assigned to each subhalo.

    Authors: Andrew Bowell & Peter Thomas
    Latest edit: 9-Nov-20
        
   """
   
   def __init__(self,graph_ID,open_graph_file,parameters):
      """ 
      Opening HDF datasets for a graph and saving them to the class.
        
      Parameters
      ----------
      graph_ID : int
         The ID of the graph to be opened.
      open_graph_file : :obj: 'File'
         Open HDF5 file that is to be read from.
      model_params : obj: 'Class'
         parameters class object.  
         Included in case graph file has parameters to be read.
      
      """
    
      self.graph_ID = graph_ID

      graph = open_graph_file[str(graph_ID)]
      # Properties of graph: halos per snaphot (generation)
      self.generation_id = graph['generation_id'][:]
      self.generation_length = graph['generation_length'][:]
      self.generation_start_index = graph['generation_start_index'][:]
      # Variable length arrays (because of possible graph branching)
      self.direct_desc_contribution = graph['direct_desc_contribution'][:]
      self.direct_desc_ids = graph['direct_desc_ids'][:]
      self.direct_prog_contribution = graph['direct_prog_contribution'][:]
      self.direct_prog_ids = graph['direct_prog_ids'][:]
      # Halo properties
      self.rms_speed = graph['3D_velocity_dispersion'][:]
      self.desc_start_index = graph['desc_start_index'][:]
      self.half_mass_radius = graph['half_mass_radius'][:]
      self.half_mass_speed = graph['half_mass_velocity_radius'][:]
      self.halo_catalog_halo_ids = graph['halo_catalog_halo_ids'][:]
      self.mean_pos = graph['mean_pos'][:]
      self.ndesc = graph['ndesc'][:]
      self.nparts = graph['nparts'][:]
      self.nprog = graph['nprog'][:]
      self.prog_start_index = graph['prog_start_index'][:]
      self.redshifts = graph['redshifts'][:]
      self.rms_radius = graph['rms_radius'][:]
      self.snapshots = graph['snapshots'][:]
      self.v_max = graph['v_max'][:]
      self.n_halos = len(self.nparts)
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
         self.nsubhalos=np.zeros(self.n_halos)
         
