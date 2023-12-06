Change log
==========

Log of major changes to the code.

28-Nov-23
---------

Have t​​wo working versions in separate branches:

* main – works with existing MEGA input.  That is deficient in many ways and the input format needs to be changed.  This has been done in...
* new_input_format – required entries in the graph input file are much better defined with cleaner names.  A helper function has been written to convert Millennium output tree files into that format (but without the ability to trace orphaned halos).
  
The basic code contains the following astrophysics:

* Pushing of material from one snapshot  to the next: halos, subhalos, galaxies.
* Accretion of baryons onto halos.
* Reaccretion of ejected gas.
* Cooling of gas from halos onto subhalos and subhalos onto galaxies.
* Merging of galaxies within subhalos, with possible starburst.
* Star formation and feedback.
* AGN feedback.
  
The code currently has the following complexity:

* Star formation histories.
* **Not** resolved discs.
* **Not** complex metallicities, just a single overall metallicity.
  
Preliminary routines have been written to analyse Halo, subhalo and galaxy output.  These show that the code performs broadly as expected but with some outliers: I suspect that this is a problem with the input graphs/trees from Millennium.

Routines have been implemented to monitor cpu memory usage (see separate document).  It is fair to say that both of these seem very poor at the moment compared to a pure C implementation and the next step is to improve this to demonstrate that the code will have a practical use.
