To-do
=====

Optimisation
------------

This is a priority for **v0.4** of the code prior to wider distribution.

* Convert push_snap to C.  This will require switching to using numpy structured arrays for graphs.
* Eliminate classes from halo, subhalo and galaxy output buffers: convert routines to C.
* Better timing tests needed to determine the efficiency of the code, ideally by comparison with the release version of L-Galaxies.
  

Major structural development
----------------------------

* Implement MCMC to optimise parameters
* Implement multiple metallicities
* Implement resolved discs

Pre-processing development
--------------------------

* Write helper routines to convert from other tree formats to py-gal input.  A version now exists for Millennium trees (which does, however, still generate a very few broken graphs).

Minor structural development
----------------------------

* (Possibly) cool from subhalos onto the central galaxy every galaxy timestep (currently done every halo timestep).

* Option to pre-process merger graph to convert it to a tree:
  
  - for comparison with other models, specifically L-Galaxies
  - could possibly be a run-time option within L-Galaxies (give everything to main descendant) but should work out of the box with tree as input.
  
* Code to locate the central subhalo (if any) in a halo.

* Need angular momentum of halos from merger graph:

  - follow angular momentum of gas and stellar discs.
  - investigate using shrinking of disc in mergers to trigger starburst.

Galaxy physics modules to import from L-galaxies
------------------------------------------------

* Reincorporate gas onto halos:
  Basic version completed
  
  - implement infall from Ejected phase using new prescription.

* Hot gas cooling onto galaxy:
  Basic SIS version completed and tested.

  - need better version of merger graphs so can use angular momentum of halo.
  - implement a beta model for the hot gas distribution.
  - add resolved galactic discs
  
* Star formation and feedback:
  Unresolved version implemented
  
  - add resolved galactic discs with inflow of gas

* AGN accretion and feedback:
  BH growth implemented
  
  - Add in calculation of quasar luminosity.
  - Implement feedback models
  
* Galaxy merging:
  Basic version with instantaneous merging implemented
  
  - ideally need a way of determining merger time.
  - triggering starburst (does this need to be explicit, or can it arise naturally from contraction of disc?).

* Stripping:

  - Need to implement stripping of non-central galaxies (and subhalos).

Plotting developments
---------------------

* Generic interface to add observations to plots
* Galaxy stellar mass function
* Luminosity functions (requires code to generate SEDs from star formation history)
* Stellar to halo mass ratio

Documentation
-------------

.. * Eliminate warnings and errors in documentation generation.

* Add in API for C code.
* Check through documentation before public release as some may have become outdated.

Pedagogy
--------

* Write a galaxy formation primer that gradually switches on the astrophysics, one step at a time, and illustrates the effect on the galaxy population.

Testing
-------

* Determine the cause of halos with excessively large baryon fractions -- should probably wait until we have better input graphs.

Known issues
------------

* code-helper/Millennium_to_pygal.py

  - produces the very occasional halo that is not in the final snapshot yet has halo.desc_start_gid=-1 (these should have been eliminated).  A specific example is Millennium input tree #5, graph_ID = 122 , snap_ID = 47 , halo_gid = 945.

  - produces many halos with baryon fractions exceeding the universal value, including one as large as 15 -- review the baryon fraction calculation.
