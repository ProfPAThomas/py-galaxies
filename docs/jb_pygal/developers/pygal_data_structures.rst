py-galaxies data structures: halos, subhalos and galaxies
=========================================================

Halos, subhalos and galaxies are the main data structures within py-galaxies.  Halos and subhalos are defined as an instance of the corresponding class.  Each carries a large number of attributes and methods which can be viewed in the API interface.  Galaxies are stored as a numpy array, so as to allow easy and rapid interface to the C-code astrophysics functions.

Some terminology
----------------

py-galaxies builds galaxies within distinct **merger graphs** whose nodes are located at a sequence of **snapshots**.  Unlike merger trees, where nodes have a single descendent, graphs can branch and lead to multiple descendants.

The nodes are constructed from friends of friends halos using two different linking lengths, corresponding to two different overdensities.  This leads to a structure in which one set of nodes is spatially nested within the other.
The large linking length, (relatively) low density nodes are **halos**: these correspond to virialised, collapsed regions (groups, clusters).  The short linking length, high density nodes are **subhalos**: these are the locations where galaxies may form.  Ideally, each halo in the input graph should have at least one subhalo; however the code provides the option to create dummy subhalos with properties identical to that of the halo, should any subhalos be missing, by setting :code:`b_create_missing_subhalos` to :code:`True` in the input yaml file; if this input parameter is absent or set to :code:`False` then the code will otherwise abort.  The halo within which each subhalo resides is the **host halo**; within any given subhalo this can simply be referred to as 'the' halo.  Usually, but not always, there will be a subhalo close to the dynamical centre of each halo: this is the **central subhalo**; the others are **satellite subhalos**.

**Galaxies** form initially within subhalos, but may become **orphans** if their subhalo is tidally destroyed in the subsequent evolution; if that happens, they are passed on to the relevant host halo.

Through merging, it is possible that subhalos may contain more than one galaxy.  One of these will be deemed to be the **central galaxy**; the others are **satellite galaxies**.

Data structures
---------------

The input data for graphs and their enclosed halos and subhalos is designed to match that of MEGA.  Helper scripts are provided in the folder code-helper to convert other tree data into the required format.  The data format is described in the section on :doc:`mega`.

Each graph is handled separately.

Processing loops over snapshots, from past to present.

Within each snapshot, halos and subhalos are stored as python lists of halo and subhalo class instances.  When passing forward information from one snapshot to the next, two copies of these lists are required: one for the progenitor and one for the new snapshot; for subsequent processing, only the current snapshot is required.

halos and subhalos are therefore stored in two different kinds of structures within which they are indexed differently:

* MEGA input files and graphs: ``_gid`` is the location relative to the graph datasets.
* snapshot lists: ``_sid`` is the location relative to the first halo/subhalo in the list.

These are related by constant offsets: ``_gid = graph.snap_first_(sub)_gid[snap_ID] + _sid``

Galaxies are stored in numpy arrays, with again 2 copies being required.  Unlike (sub)halos, their number can vary, and they can move between subhalos and halos.  Therefore, the ordering of galaxies from one snapshot to the next may change, and some book-keeping is required to keep track of which galaxies belong to which (sub)halo.

The API section contains full details of the graph, halo, subhalo and galaxy data structures.  Here we just outline the main variables used in book-keeping.

Graphs
^^^^^^

* Graph and snapshot properties:
  
  *  graph_ID (int) : The unique identifier of the graph within the file.
  *  n_gal (int) : Running total of the number of galaxies created by py-galaxies when processing the graph.
  *  snap_n_halo (int[n_snap]) : The number of halos in each snapshot.
  *  snap_first_halo_gid (int[n_snap]) : The first halo in each snapshot, relative to the graph.
  *  snap_n_sub (int[n_snap]) : The number of subhalos in each snapshot.
  *  snap_first_sub_gid (int[n_snap]) : The first subhalo in each snapshot, relative to the graph.

* Halo properties:
  
  *  halo_desc_contribution (float[halo_n_desc]) : The particle contributions to descendant halos.
  *  halo_desc_IDs_gid (int[halo_n_desc]) : The locations in the graph of all descendant halos of each halo.
  *  halo_first_desc_gid (int[n_halo]) : The locations in the graph of the first descendant halo of each halo.
  *  halo_first_sub_gid (float[n_halo]) : The locations in the graph of first subhalo in each halo.
  *  halo_n_desc (float[n_halo]) : The number of descendants of each halo.
  *  halo_n_sub (float[n_halo]) : The number of subhalos of each halo

* Subhalo properties:
  
  *  sub_desc_contribution (float[n_desc_sub]) : The particle contributions to descendant subhalos.
  *  sub_desc_IDs_gid  (int[n_desc_sub]) : The locations in the graph of all descendant subhalos.
  *  sub_first_desc_gid (int[n_sub]): The locations in the graph of the first descendant subhalo of each subhalo.
  *  sub_host_gid (int[n_sub]) : The host halo of each subhalo.
  *  sub_n_desc (float[n_sub]) : The number of descendants of each subhalo.

Halos
^^^^^

Halos are instances of the C_halo class, contained in a python list ``halos_this_snap``.  At the end of processing each snapshot, this is copied into ``halos_last_snap`` to enable forward propagation between snapshots.  Thus a maximum of two snapshots are required at any given time.

We mention here only those attributes that are related to traversal of the resulting subhalo and galaxy structures.

Galaxies create the most book-keeping.  Within the galaxy array, we first store galaxies associated with each subhalo of a given halo, then follow that by the orphan halos associated with that halo.  Orphans can be passed on from progenitors, but also from subhalos in the previous snapshot that have no descendants.

* Halo identifiers:

  - graph_ID (str) : The graph_ID (from HDF5 group).
  - snap_ID (int) : The snapshot ID currently being processed.
  - halo_gid (int) : The halo location within the graph of the halo currently being processed.
  - halo_sid (int) : The halo location within the snapshot of the halo currently being processed.

* Subhalos:

  - n_sub (int) : The number of subhalos of this halo.
  - sub_start_(gid|sid) : The location of the first subhalo of this halo.
  - sub_end_(gid|sid) : The location of the last subhalo of this halo +1 (python indexing).
  - sub_central_(gid|sid) : The central subhalo of this halo, if any.

* Descendants:

  - n_desc (int) : The number of descendant halos.
  - desc_start_gid (int) : The location of the first descendant halo in the descendant arrays
  - desc_end_gid (int) : The location of the last descendant halo in the descendant arrays +1 (python indexing)
  - desc_main_sid (int) : The location of the main (most massive) descendant halo in the following snapshot halo list, for galaxy propagation.
    
* Galaxies:

  - n_gal (int) : The number of galaxies in this halo + all subhalos.
  - gal_start_(sid|gid) (int) : The location of the first galaxy in this halo relative to the snapshot and to the graph.
  - n_orphan (int) : the number of galaxies not associated with subhalos.
  - orphan_start_sid (int) : The location of the first orphan galaxy in the galaxies array.
  - orphan_next_sid (int) : The next available location to store orphans in the galaxy array.

 
Subhalos
^^^^^^^^

Subhalos are instances of the C_sub class, contained in a python list ``subs_this_snap``.  At the end of processing each snapshot, this is copied into ``subs_last_snap`` to enable forward propagation between snapshots.

We mention here only those attributes that are related to traversal of the halo and galaxy structures.

* Subhalo identifiers:

  - graph_ID (str) : The graph_ID (from HDF5 group).
  - snap_ID (int) : The snapshot ID currently being processed.
  - halo_gid (int) : The location within the graph of the host halo.
  - halo_sid (int) : The location within the snapshot of the host halo.
  - sub_gid (int) : The subhalo location within the graph.
  - sub_sid (int) : The subhalo location within the snapshot.

* Descendants:

  - n_desc (int) : The number of descendant halos.
  - desc_start_gid (int) : The location of the first descendant subhalo in the descendant arrays
  - desc_end_gid (int) : The location of the last descendant subhalo in the descendant arrays +1 (python indexing)
  - desc_halo_sid (int) : The location of the descendant halo of the host halo in the following snapshot halo list, for galaxy propagation.
    
* Galaxies:

  - n_gal (int) : The number of galaxies in this halo + all subhalos.
  - gal_start_sid (int) : The location of the first galaxy in this halo relative to the snapshot.
  - gal_end_sid (int) : The location of the last galaxy in this halo relative to the snapshot +1 (python indexing).
  - gal_next_sid (int) : Galaxy counter used when updating galaxies.
  - gal_central_sid (iint) : The location in the current galaxy array of the most massive galaxy in the subhalo.

Galaxies
^^^^^^^^

Galaxies are stored in a numpy array for each snapshot.  That array contains the following entries to relate them to halos and subhalos

* Graph, halos & subhalos

  - graph_ID (str) : The graph_ID (from HDF5 group).
  - snap_ID (int) : The snapshot ID currently being processed.
  - halo_gid (int) : The location within the graph of the host halo.
  - halo_sid (int) : The location within the snapshot of the host halo.
  - sub_gid (int) : The subhalo location within the graph.
  - sub_sid (int) : The subhalo location within the snapshot.

We also give galaxies unique labels within the graph, and track their merging tree (as galaxies cannot split, we do not need a graph).  As we have a tree, we can use the usual pointers, stored within the galaxy array.
    
* Galaxies
  
  - gal_gid (int) : The location of this galaxy within all galaxies for this graph.
  - desc_gid (int) : The location of the descendant galaxy within all galaxies for this graph.
  - first_prog_gid (int) : The location within all galaxies for this graph of the progenitor galaxy in the previous snapshot from which this galaxy derives.
  - next_prog_gid (int) : The location within all galaxies for this graph of the next galaxy in this snapshot that merges with the same descendant.


Structure of the output files
-----------------------------

Information is stored in separate HDF5 files each containing a numpy structured array for halos, subhalos and galaxies.
Note that the units of quantities are those specified in the input.yml file.

Halos
^^^^^

This dataset is labelled "Halos".

========================  ================  =======================================================================
   Name                     Type                  Description
========================  ================  =======================================================================
graph_ID                    int32            The graph in which this halo resides
snap_ID                     int32            The snapshot in which this halo resides
halo_gid                    int32            The number/location of the halo within the graph
pos                        (float32,(3,))    The location of the halo
vel                        (float32,(3,))    The mean velocity of the halo
mass                        float32          The dark matter mass of the halo
temperature                 float32          The temperature of the halo assuming an SIS halo model
rms_speed                   float32          The 3-D rms speed of DM particles in the halo
half_mass_virial_speed      float32          The circular speed at the half mass radius assuming an SIS halo model
mass_baryon                 float32          The baryon mass in the halo, including all subhalos and galaxies
mass_gas_hot                float32          The hot gas mass in the halo, excluding subhalos
mass_metals_gas_hot         float32          The hot gas metal mass in the halo, excluding subhalos
mass_gas_eject              float32          The mass of gas ejected from the halo
mass_metals_gas_eject       float32          The moss of metals in the gas ejected from the halo
mass_stars                  float32          The (initial) stellar mass in the halo, excluding subhalos
mass_metals_stars           float32          The (initial) stellar metals mass in the halo, excluding subhalos
========================  ================  =======================================================================

Subhalos
^^^^^^^^

This dataset is labelled "Subhalos"

========================  ================  =======================================================================
   Name                     Type                  Description
========================  ================  =======================================================================
graph_ID                    int32            The graph in which this halo resides
snap_ID                     int32            The snapshot in which this halo resides
halo_gid                    int32            The number/location of the host halo within the graph
sub_gid                     int32            The number/location of the subhalo within the graph
pos                        (float32,(3,))    The location of the halo
vel                        (float32,(3,))    The mean velocity of the halo
mass                        float32          The dark matter mass of the halo
temperature                 float32          The temperature of the halo assuming an SIS halo model
rms_speed                   float32          The 3-D rms speed of DM particles in the halo
half_mass_virial_speed      float32          The circular speed at the half mass radius assuming an SIS halo model
mass_gas_hot                float32          The hot gas mass in the halo, excluding subhalos
mass_metals_gas_hot         float32          The hot gas metal mass in the halo, excluding subhalos
mass_stars                  float32          The (initial) stellar mass in the halo, excluding subhalos
mass_metals_stars           float32          The (initial) stellar metals mass in the halo, excluding subhalos
========================  ================  =======================================================================

Galaxies
^^^^^^^^

This dataset is labelled "Galaxies"

.. Comments don't work in the table :-(

========================  ================  =======================================================================
   Name                     Type                  Description
========================  ================  =======================================================================
graph_ID                    int32            The graph in which this halo resides
snap_ID                     int32            The snapshot in which this halo resides
halo_gid                    int32            The number/location of the host halo within the graph
sub_gid                     int32            The number/location of the host subhalo within the graph
gal_gid                     int32            The number/location of the galaxy within the graph
desc_gid                    int32            The number/location of the descendant galaxy within the graph
first_prog_gid              int32            The number/location of the first progenitor galaxy within the graph
next_prog_gid               int32            If looping over progenitors, this is a pointer to the next one
b_exists                    bool             Galaxy exists (otherwise has merged and should be ignored
#pos                        (float32,(3,))    The location of the galaxy
#vel                        (float32,(3,))    The mean velocity of the galaxy
mass_stars_bulge            float32          The (initial) stellar mass in the bulge
mass_metals_stars_bulge     float32          The (initial) stellar metal mass in the bulge
mass_stars_disc             float32          The (initial) stellar mass in the disc
mass_metals_stars_disc      float32          The (initial) stellar metals mass in the disc
mass_gas_cold               float32          The mass in the cold gas (ie ISM)
mass_metals_gas_cold        float32          The metals mass in the cold gas (ie ISM)
mass_BH                     float32          The mass of the central  black hole
radius_gas_cold             float32          The exponential disc scale length of the cold gas
radius_stars_disc           float32          The exponential disc scale length of the stellar disc
radius_stars_bulge          float32          The half-mass radius of the stellar bulge
SFR_dt                      float32          The star formation rate in the last galaxy timestep
SFR_snap                    float32          The star formation rate averaged over the last snapshot
========================  ================  =======================================================================

The pointers to locations in the galaxy table refer to the current graph.  Therefore we also need a record of where each graph starts.  The relevant dataset is labelled "Graph_start_locations"; it is a 1-D numpy array.

========================  ================  =======================================================================
   Name                     Type                  Description
========================  ================  =======================================================================
<None>                      int32            The location within the Galaxies dataset where each graph starts
========================  ================  =======================================================================
