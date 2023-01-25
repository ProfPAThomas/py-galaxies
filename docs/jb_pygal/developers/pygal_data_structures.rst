py-galaxies data structures: halos, subhalos and galaxies
=========================================================

Some terminology
----------------

py-galaxies builds galaxies within distinct **merger graphs** whose nodes are located at a sequence of **snapshots**.  Unlike merger trees, where nodes have a single descendent, graphs can branch and lead to multiple descendants.

The nodes are constructed from friends of friends halos using two different linking lengths, corresponding to two different overdensities.  This leads to a structure in which one set of nodes is spatially nested within the other.

The large linking length, (relatively) low density nodes are **halos**: these correspond to virialised, collapsed regions (groups, clusters).  The short linking length, high density nodes are **subhalos**: these are the locations where galaxies may form.  The halo within which each subhalo resides is  the **host halo**.  Usually, but not always, there will be a subhalo close to the dynamical centre of each halo: this is the **central subhalo**; the others are **satellite subhalos**.

**Galaxies** form initially within subhalos, but may become orphans if their subhalo is tidally destroyed in the subsequent evolution; if that happens, they are passed on to the relevant host halo.

Through merging, it is possible that subhalos may contain more than one galaxy.  One of these will be deemed to be the **central galaxy**; the others are **satellite galaxies**.

Data structures
---------------

The input data for graphs and their enclosed halos and subhalos comes from MEGA.  The datasets are described in the subsection on Graphs below.

Each graph is handled separately.

Processing loops over snapshots, from past to present.

Within each snapshot, halos, subhalo and galaxies are stored as python lists.  When passing forward information from one snapshot to the next, two copies of these lists are required: one for the progenitor and one for the new snapshot; for subsequent processing, only the current snapshot is required.

halos and subhalos are therefore stored in two different kinds of structures within which they are indexed differently:

* MEGA input files: _gid is the location relative to the graph datasets.
* snapshot lists: _sid is the location relative to the first halo/subhalo in the list

These are related by constant offsets:

* halo_offset = graph.halo_start_gid[snap_ID]
* sub_offset = graph.sub_start_gid[snap_ID]

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


