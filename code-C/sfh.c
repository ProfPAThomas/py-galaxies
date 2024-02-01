/* sfh.c

Functions to merger SFH bins if required.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "halos.h"
#include "subs.h"
#include "gals.h"
#include "sfh.h"
#include "parameters.h"
// proto.h has to come last in order not to generate warnings about multiple struct definitions
// Could probably get around that by using #ifndef.
#include "proto.h"

void F_sfh_update_bins(struct struct_gal gals[], int n_gal, int i_dt) {
    /*
      Merges galaxy bins if required.

      Arguments
      ---------
      gals : obj : np.array[n_gal]
         Structured array of galaxies
      n_gal : int
         Size of galaxy array
      i_dt : int
         ministep ID BEFORE updating
    */

    int i_level;
    int i_bin, j_bin, k_bin, i_gal;    
    int n_bin_in_level[n_level];
    int level[n_bin+1];

    // Number of bins used in this timestep
    i_bin=i_bin_all[i_dt];
    // Set number of bins at each level in the merging hierarchy
    for (i_level=0; i_level<n_level; i_level++) n_bin_in_level[i_level]=n_bin_in_level_all[i_dt][i_level];
    // Set level in merging hierarchy of each SFH bin
    for (j_bin=0; j_bin<n_bin; j_bin++) level[j_bin]=level_all[i_dt][j_bin];
    level[n_bin]=parameters.NO_DATA_INT;
    // Create new bin
    level[i_bin]=0;
    n_bin_in_level[0]+=1;
    // Run through levels merging as required
    j_bin=i_bin+1;
    for (i_level=0; i_level<n_level; i_level++) {
	j_bin-=n_bin_in_level[i_level];
	if (n_bin_in_level[i_level] == n_merge) {
	    // Need to merge bins j_bin and j_bin+1
	    // This essential here means first resetting the times and counts of bins at each level,
	    // then shuffling bins downward.
	    i_bin-=1;
	    n_bin_in_level[i_level+1]+=1;
	    n_bin_in_level[i_level]-=2;
	    level[j_bin]+=1;
	    for (k_bin=j_bin+1; k_bin<n_bin; k_bin++) level[k_bin]=level[k_bin+1];
	    level[n_bin]=parameters.NO_DATA_INT;
	    // Now combine the data.
	    // Would this be faster (and make the code look simpler) if all the SFH data was a sub-array?
	    for (i_gal=0; i_gal<n_gal; i_gal++) {
		(gals[i_gal]).mass_stars_bulge_sfh[j_bin]+=(gals[i_gal]).mass_stars_bulge_sfh[j_bin+1];
		(gals[i_gal]).mass_metals_stars_bulge_sfh[j_bin]+=(gals[i_gal]).mass_metals_stars_bulge_sfh[j_bin+1];
		(gals[i_gal]).mass_stars_disc_sfh[j_bin]+=(gals[i_gal]).mass_stars_disc_sfh[j_bin+1];
		(gals[i_gal]).mass_metals_stars_disc_sfh[j_bin]+=(gals[i_gal]).mass_metals_stars_disc_sfh[j_bin+1];
		for (k_bin=j_bin+1; k_bin<n_bin; k_bin++) {
		    (gals[i_gal]).mass_stars_bulge_sfh[k_bin]=(gals[i_gal]).mass_stars_bulge_sfh[k_bin+1];
		    (gals[i_gal]).mass_metals_stars_bulge_sfh[k_bin]=(gals[i_gal]).mass_metals_stars_bulge_sfh[k_bin+1];
		    (gals[i_gal]).mass_stars_disc_sfh[k_bin]=(gals[i_gal]).mass_stars_disc_sfh[k_bin+1];
		    (gals[i_gal]).mass_metals_stars_disc_sfh[k_bin]=(gals[i_gal]).mass_metals_stars_disc_sfh[k_bin+1];
		}
		(gals[i_gal]).mass_stars_bulge_sfh[n_bin]=0.;
		(gals[i_gal]).mass_metals_stars_bulge_sfh[n_bin]=0.;
		(gals[i_gal]).mass_stars_disc_sfh[n_bin]=0.;
		(gals[i_gal]).mass_metals_stars_disc_sfh[n_bin]=0.;
	    }
	    j_bin+=1;      // First bin at this level has been pushed upwards by one slot
	}
    }
    return;
}
      
