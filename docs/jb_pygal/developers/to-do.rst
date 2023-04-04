To-do
=====

Major structural development
----------------------------

* Add in performance monitoring & test timing.
* Work out how to interface to C-routines
* Implement MCMC to optimise parameters
* Implement multiple metallicities
* Implement resolved discs
* Implement star formation and metal enrichment histories
* Write a galaxy formation primer that gradually switches on the astrophysics, one step at a time, and illustrates the effect on the galaxy population.

Minor structural development
----------------------------

* Pre-processing merger graph to:
  
  - ensure every halo has a subhalo (introduce dummy subhalos as necessary)
  - eliminate halos that have no descendants

* Option to pre-process merger graph to convert it to a tree:
  
  - for comparison with other models, specifically L-Galaxies
  - could possibly be a run-time option within L-Galaxies (give everything to main descendant) but should work out of the box with tree as input.
  
* Code to locate the central subhalo (if any) in a halo.

* (Probably) code to add a dummy subhalo if running in L-Galaxies model (or even if not, might allow higher resolution).

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
  - cool in stages and check that this better aligns the instantaneous and averaged cooling rates.
  
* Star formation and feedback:
  
  - implement a version with no resolved discs:

    + basic implementation done
    + need to add in calculation of SFR
    + need to test
    
  - add resolved galactic discs with inflow of gas

* AGN accretion and feedback:
  
  - BH growth from mergers: add in calculation of quasar luminosity.
  
* Galaxy merging:
  Basic version with instantaneous merging implemented
  
  - need better merger graphs so 
  - triggering starburst (does this need to be explicit, or can it arise naturally from contraction of disc?)

* Stripping:

  - Need to implement stripping of non-central galaxies (and subhalos).

Plotting developments
---------------------

* Galaxy stellar mass function
* Luminosity functions (requires code to generate SEDs from star formation history)

Testing
-------

* Check where the crazily large gas disc sizes come from.
