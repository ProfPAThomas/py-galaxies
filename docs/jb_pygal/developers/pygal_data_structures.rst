py-galaxies data structures: halos, subhalos and galaxies
=========================================================

Halos, subhalos and galaxies are the main data structures within py-galaxies.  Halos and subhalos are defined as an instance of the corresponding class.  Each carries a large number of attributes and methods which can be viewed in the API interface.  Galaxies are stored as a numpy array, so as to allow easy and rapid interface to the C-code astrophysics functions.

Some terminology
----------------

py-galaxies builds galaxies within distinct **merger graphs** whose nodes are located at a sequence of **snapshots**.  Unlike merger trees, where nodes have a single descendent, graphs can branch and lead to multiple descendants.

The nodes are constructed from friends of friends halos using two different linking lengths, corresponding to two different overdensities.  This leads to a structure in which one set of nodes is spatially nested within the other.

The large linking length, (relatively) low density nodes are **halos**: these correspond to virialised, collapsed regions (groups, clusters).  The short linking length, high density nodes are **subhalos**: these are the locations where galaxies may form.  The halo within which each subhalo resides is the **host halo**; within any given subhalo this can simply be referred to as 'the' halo.  Usually, but not always, there will be a subhalo close to the dynamical centre of each halo: this is the **central subhalo**; the others are **satellite subhalos**.

**Galaxies** form initially within subhalos, but may become orphans if their subhalo is tidally destroyed in the subsequent evolution; if that happens, they are passed on to the relevant host halo.

Through merging, it is possible that subhalos may contain more than one galaxy.  One of these will be deemed to be the **central galaxy**; the others are **satellite galaxies**.

Data structures
---------------

The input data for graphs and their enclosed halos and subhalos comes from MEGA.  The datasets are described in the subsection on Graphs below.

Each graph is handled separately.

Processing loops over snapshots, from past to present.

Within each snapshot, halos, subhalo and galaxies are stored as python lists.  When passing forward information from one snapshot to the next, two copies of these lists are required: one for the progenitor and one for the new snapshot; for subsequent processing, only the current snapshot is required.

halos and subhalos are therefore stored in two different kinds of structures within which they are indexed differently:

* MEGA input files: ``_gid`` is the location relative to the graph datasets.
* snapshot lists: `_sid`` is the location relative to the first halo/subhalo in the list

These are related by constant offsets:

* ``halo_offset = graph.halo_start_gid[snap_ID]``
* ``sub_offset = graph.sub_start_gid[snap_ID]``

Galaxies are only stored in numpy arrays.  Moreover, their number can vary, and they can move between subhalos and halos.  Therefore, the ordering of galaxies from one snapshot to the next may change, and some book-keeping is required to keep track of which galaxies belong to which (sub)halo.

Graphs
^^^^^^

The MEGA output stores data within graphs.  Each graph then contains a (large) number of datasets.  We list here only the ones required for traversal of the resulting halo and subhalo structures.

* Header:
  These are not currently read in but could be used to trim the set of graphs that are processed.
  
  - graph_lengths: the number of shapshots of each graph that contain haloss.
  - root_nparts: the number of particles in the root (final snapshot) of each graph.
  - nhalos_in_graph: total number of halos in each graph.
  - sub_graph_lengths: the number of shapshots of each graph that contain subhalos.
  - sub_nhalos_in_graph: total number of subhalos in each graph.
    
* Graphs [listed as numbers]
  Note that the subhalo datasets exist only if there are subhalos present within the graph.  We list here the *name within the MEGA file* and the **internal name within the graph class** in py-galaxies.
  
  - attributes:
    
    + *length* - number of snapshots with halos in.
    + *nhalos_in_graph* | **n_halo** - the number of halos in the graph.
    + *root_mass* | **root_mass** - the number of particles in the most massive halo in the final snapshot.
    + *sub_length* - the number of snapshots that contain subhalos.
    + *sub_nhalos_in_graph* | **n_sub** - the number of subhalos in the graph.
    + *sub_root_mass* - the number of particles in the most massive subhalo in the final snapshot.

  - *generation_length* | **n_halo_snap** - the number of halos in each snapshot (set to NO_DATA_INT if no halos).
  - *generation_start_index* | **halo_start_gid** - the first halo in each snapshot.
  - *ndesc* | **n_desc** - the number of descendents of each halo.
  - *desc_start_index* | **desc_start_gid** - the first descendent halo (within the following dataset - the descendent halos are not guaranteed to be in halo order, or even unique)).
  - *direct_desc_ids* | **desc_IDs_gid** - a list of descendant halos for each halo.
  - *nprog* | **n_prog** - the number of progenitors of each halo.
  - *prog_start_index* | **prog_start_gid** - the first progenitor halo.
  - *sub_generation_length* | **n_sub_snap** - the number of subhalos in each snapshot.
  - *sub_generation_start_index* | **sub_start_gid** - the first subhalo for each snapshot.
  - *nsubhalo* | **n_sub_halo** - the number of subhalos within each halo.
  - *subhalo_start_index* | **sub_start_gid** - the first subhalo within each halo.
  - *host_halo* | **sub_host_gid** - the host halo of each subhalo.
  - *sub_ndesc* | **n_desc** - the number of descendants of each subhalo.
  - *sub_desc_start_index* | **sub_desc_start_gid** - the first descendant subhalo for each subhalo (within the following dataset - the descendent subhalos are not guaranteed to be in subhalo order, or even unique).
  - *sub_direct_desc_ids* | **sub_desc_IDs_gid** - a list of descendent subhalos for each subhalo.

Halos
^^^^^

Halos are instances of the halo class, contained in a python list ``halos_this_snap``.  At the end of processing each snapshot, this is copied into ``halos_last_snap`` to enable forward propagation between snapshots.  The relationship between location in the graph files ``_gid`` and the snapshot lists ``_sid`` is a simple offset ``halo_offset = graph.halo_start_gid[snap_ID]``.

We mention here only those attributes that are related to traversal of the resulting subhalo and galaxy structures.

Galaxies create the most book-keeping.  Within the galaxy array, we first store galaxies associated with each subhalo of a given halo, then follow that by the orphan halos associated with that halo.  Orphans can be passed on from progenitors, but also from subhalos in the previous snapshot that have no descendants.

* Halos:

  - halo_offset - to relate location in graph to that in snapshot
  - n_desc - number of descendant halos
  - desc_start_gid - first descendant halo
  - desc_end_gid - last descendant halo +1 (python indexing)
  - desc_main_sid - main (most massive) descendant halo, for galaxy propagation

* Subhalos:

  - sub_offset - relate location in graph to that in snapshot
  - n_sub - number of subhalos
  - sub_start_(gid|sid) - first subhalo
  - sub_end_(gid|sid) - last subhalo +1 (python indexing)
  - sub_central_(gid|sid) - the central subhalo of this halo, if any

* Galaxies:

  - n_gal - the number of galaxies in this halo + all subhalos.
  - gal_start_(sid|gid) - location of the first galaxy in this halo in the galaxies (snapshot list | graph).
  - n_orphan - the number of galaxies not associated with subhalos.
  - orphan_start_sid - the location of the first orphan galaxy in the galaxies snapshot list.
  - orphan_next_sid - keeps track of the next available location to store orphans in the galaxy array.

 
Subhalos
^^^^^^^^

Subhalos are instances of the sub class, contained in a python list ``subs_this_snap``.  At the end of processing each snapshot, this is copied into ``subs_last_snap`` to enable forward propagation between snapshots.  The relationship between location in the graph files ``_gid`` and the snapshot lists ``_sid`` is a simple offset ``sub_offset = graph.sub_start_gid[snap_ID]``.

We mention here only those attributes that are related to traversal of the halo and galaxy structures.

Galaxies create the most book-keeping.  Within the galaxy array, we first store galaxies associated with each subhalo of a given halo, then follow that by the orphan halos associated with that halo.  Orphans can be passed on from progenitors, but also from subhalos in the previous snapshot that have no descendants.

* Halos:

  - host_ID_gid - the host halo of this subhalo
  - desc_host_sid - main descendant of host halo; needed if no descendant of subhalo

* Subhalos
    
  - n_desc - the number of descendant subhalos.
  - desc_start_gid - first descendant subhalo
  - desc_end_gid - last descendant subhalo +1 (python indexing)
  - desc_host_sid - the descendant halo of the host which will receive this subhalos contents, should this subhalo have no descendant.

* Galaxies

  - n_gal - the number of galaxies in this halo + all subhalos.
  - gal_start_sid - location of the first galaxy in this subhalo in the galaxies snapshot list.
  - gal_next_sid - keeps track of the next available location to store inherited galaxies in the galaxy array.

Galaxies
^^^^^^^^

Galaxies are stored in a numpy array for each snapshot.  That array contains the following entries to relate them to halos and subhalos

* Graph, halos & subhalos

  - halo_(gid|sid) - the location of the host halo in the graph and snapshot datasets.
  - sub_(gid|sid) - the location of the host subhalo in the graph and snapshot datasets.

We also give galaxies unique labels within the graph, and track their merging tree (as galaxies cannot split, we do not need a graph).  As we have a tree, we can use the usual pointers, stored within the galaxy array.
    
* Galaxies
  
  - gal_gid - overall counter of the id of this galaxy within all galaxies for this graph.
  - desc_gid - the descendant galaxy in the next snapshot.
  - first_prog_gid - the galaxy in the previous snapshot from which this galaxy derives.
  - next_prog_gid - the next galaxy in this snapshot that merges with the same descendant.


