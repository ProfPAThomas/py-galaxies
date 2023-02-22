To-do
=====

Major structural development
----------------------------

* Work out how to interface to C-routines
* Implement MCMC to optimise parameters
* Implement star formation and metal enrichment histories
* Write a galaxy formation primer that gradually switches on the astrophysics, one step at a time, and illustrates the effect on the galaxy population.

Minor structural development
----------------------------

* Go through routines and document headers in numpy format; create API interface in documentation.

* Test that:
  
  - Where galaxies have no subhlao, it is because they are orphans.
  - Where galaxies have no halo, it is because the are orphans for which no descendant halo exists.

* Implement mini-steps (might be a major development)

* Pre-processing merger graph to:
  
  - ensure every halo has a subhalo (introduce dummy subhalos as necessary)
  - eliminate halos that have no descendants

* Option to pre-process merger graph to convert it to a tree:
  
  - for comparison with other models, specifically L-Galaxies
  - could possibly be a run-time option within L-Galaxies (give everything to main descendant) but should work out of the box with tree as input.
  
* Option to alter structure to mimic L-Galaxies as closely as possible (might be a major development).

* Code to locate the central subhalo (if any) in a halo.

* (Probably) code to add a dummy subhalo if running in L-Galaxies model (or even if not, might allow higher resolution).

Galaxy physics modules to import from L-galaxies
------------------------------------------------

* Reincorporate gas onto halos:
  Basic version completed
  
  - decide whether we need an **Ejected** phase for halos; I suspect that we do.

* Hot gas cooling onto galaxy:
  Basic SIS version completed but not tested.
  
  - check that cooling rates have been correctly modelled - this will involve generating a test halo.
  - implement a beta model for the hot gas distribution.
  - add resolved galactic discs
  
* Star formation and feedback:
  
  - implement a version with no resolved discs
  - add resolved galactic discs with inflow of gas

* AGN accretion and feedback:
  
  - BH growth from mergers: quasar luminosity?
  - Radio mode accretion + feedback (L-Galaxies model is feeble in this regard, but sort of works)
  
* Galaxy merging:
  
  - redistributing gas and stars
  - triggering starburst (does this need to be explicit, or can it arise naturally from contraction of disc?)

Plotting developments
---------------------

* Galaxy stellar mass function
* Luminosity functions (requires code to generate SEDs from star formation history)
