.. _GraphFileStructure:

Input graph file structure
==========================

MEGA is the merger graph algorithm: https://github.com/WillJRoper/mega

The precise form of the MEGA HDF5 output is still evolving, but for the current version of py-galaxies, we require the following as a minimum:

Top level (/)
-------------

Attributes
^^^^^^^^^^

* :code:`n_graph`: the number of graphs contained in the file.

* :code:`n_snap`: the (maximum) number of snapshots for which there is graph data.

Ideally, the input file should also contain these other attributes, although these can be set within the input parameter file if they are missing or need to be over-ridden (do this with care!)

* Cosmological parameters:

  * :code:`Hubble_h`: dimensionless Hubble parameter at :code:`z=0`.

  * :code:`Omega_m`: dimensionless matter density relative to critical at :code:`z=0`.

  * :code:`baryon_fraction`: the fraction of matter contained in baryons at :code:`z=0`.

* Units.  Note that there is no facility yet to read these in automatically, so one needs to check that the input parameter file matches the values stored in the graph file.

  * :code:`unit_length`: units of length used in the input graph file.

  * :code:`unit_mass`: units of mass used in the input graph file.

  * :code:`unit_speed`: units of speed used in the input graph file.

  * :code:`unit_time`: units of time used in the input graph file.

Datasets
^^^^^^^^

* :code:`snap_table`: a table of properties for each snapshot containing at a minimum a column labelled `time_in_years` that contains the cosmic time of each snapshot as measured from the big bang.

Groups
^^^^^^

A separate group for each graph with a name of the form `graph_#`, where `#` is a non-negative integer in the range :math:`0-n_{graph}`.

Graph_#/
--------

This is a list of the required attributes and datasets for each graph.

Attributes
^^^^^^^^^^

* :code:`n_halo`: the number of halos in the graph.

* :code:`n_halo_desc`: the number of entries in the descendant halo arrays.  (May not be required, not sure!)

* :code:`n_sub`: the number of subhalos in the graph.

* :code:`n_sub_desc`: the number of entries in the descendant subhalo arrays.  (May not be required, not sure!)

Datasets
^^^^^^^^

* :code:`snap_n_halo (int[n_snap])`: The number of halos in each snapshot.

* :code:`snap_first_halo_gid (int[n_snap])`: The first halo in each snapshot, relative to the graph.

* :code:`snap_n_sub (int[n_snap])`: The number of subhalos in each snapshot.

* :code:`snap_first_sub_gid (int[n_snap])`: The first subhalo in each snapshot, relative to the graph.

* :code:`halo_desc_contribution (float[halo_n_desc])`: The particle contributions to descendant halos.

* :code:`halo_desc_IDs_gid (int[halo_n_desc])`: The locations in the graph of all descendant halos of each halo.

* :code:`halo_first_desc_gid (int[n_halo])`: The locations in the graph of the first descendant halo of each halo.

* :code:`halo_first_sub_gid (float[n_halo])`: The locations in the graph of first subhalo in each halo.

* :code:`half_mass_radius (float[n_halo])`: The radii containing half of the halo particles.

* :code:`halo_mass  (float[n_halo])`: The mass of each halo.

* :code:`halo_mean_pos (float[n_halo,3])`: The position of each halo.

* :code:`halo_mean_vel (float[n_halo,3])`: The velocity of each halo.

* :code:`halo_n_desc (float[n_halo])`: The number of descendants of each halo.

* :code:`halo_n_sub (float[n_halo])`: The number of subhalos of each halo

* :code:`halo_rms_radius (float[n_halo])`: The rms radii of halo particles

* :code:`halo_rms_speed (float[n_halo])`: The rms speeds of halo particles

* :code:`sub_desc_contribution (float[n_desc_sub])`: The particle contributions to descendant subhalos.

* :code:`sub_desc_IDs_gid  (int[n_desc_sub])`: The locations in the graph of all descendant subhalos.

* :code:`sub_first_desc_gid (int[n_sub]`: The locations in the graph of the first descendant subhalo of each

* :code:`sub_half_mass_radius (float[n_sub])`: The radii containing half of the subhalo particles.

* :code:`sub_host_gid (int[n_sub])`: The host halo of each subhalo.

* :code:`sub_mass  (float[n_sub])`: The mass of each subhalo.

* :code:`sub_n_desc (float[n_sub])`: The number of descendants of each subhalo.

* :code:`sub_pos (float[n_sub,3])`: The position of each subhalo.

* :code:`sub_rms_speed (float[n_sub])`: The rms speed of subhalo particles within each subhalo.

* :code:`sub_vel (float[n_sub,3])`: The velocity of each subhalo.






