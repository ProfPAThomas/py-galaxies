/*
F_process_snap - perform the galaxy astrophysics

*/

#include "all_headers.h"

//---------------------------------------------------------------------------------------------------

void F_process_snap(struct struct_halo halos[], int n_halo,
			    struct struct_sub subs[], int n_sub,
			    struct struct_gal gals[], int n_gal,
			    struct struct_var variables) {
    /*
    This is the controlling routine for halo processing.

    Arguments
    ---------
    halos : obj : D_halo[n_halo]
       The halos currently being processed in this graph/snapshot
    n_halo : int
       The number of such halos.
    subs : obj : D_sub[n_sub]
       The subhalos currently being processed in this graph/snapshot
    n_sub : int
       The number of such subhalos.  (Not used when subhalos done as inner loop within halos.)
    gals : obj : D_gal[n_gal]
       The galaxies currently being processed in this graph/snapshot.
    n_gal : int
       The number of such galaxies.
    variables : obj : C_variables
       Instance of class containing global variables structure

    Returns
    -------
    None
    */

    double dt_halo, mass_stars;
    int i_gal, i_gal_central, i_halo, i_sub, i_sub_central;
    int i_dt_halo, n_dt_halo, i_dt_gal, n_dt_gal;
    int sub_sid;
#ifdef SFH
    int i_dt_sfh;
#endif

    // Load timestep information
    dt_halo=variables.dt_halo;
    n_dt_halo=variables.n_dt_halo;
    n_dt_gal=variables.n_dt_gal;
#ifdef SFH
    i_dt_sfh=variables.i_dt_sfh;
#endif

    /*
      Note about the structure.
      Because subhalos contain links to their containing halo, the structure could be altered to process first all halos, 
      then all subhalos.   I have instead embedded the subhalo loop within the halo loop as that better reflects the dependencies.
      It is not clear to me which is better, if either.
      Because galaxies are updated on a different timestep, and I want to process all relevant objects within each separate increment
      of the time counters, galaxy processing takes place within the halo timestep but after all halo and subhalo processing.
    */

    // Set the accretion rate required to bring the halos up to the universal mean over the course of a snapshot interval.
    F_halos_set_mass_baryon(halos, n_halo, subs, gals, n_dt_halo);

    // Loop over halo timestep
    for (i_dt_halo=0; i_dt_halo<n_dt_halo; i_dt_halo++) {
    
	// Loop over halos
	for (i_halo=0; i_halo<n_halo; i_halo++) {

	    // Sanity check - halos should not be processed at this point
	    if (halos[i_halo].b_done==true) {
		printf("i_dt_halo = %d, i_halo = %d: halo %d in graph %d already processed\n",
		       i_dt_halo,i_halo,halos[i_halo].halo_gid,halos[i_halo].graph_ID);
		fflush(stdout);
		exit(1);
	    }
	    
	    // Accretion onto halos.
	    F_halos_accrete_primordial_gas(&halos[i_halo]);
	    // Reincorporation of ejected gas
	    if (halos[i_halo].mass_gas_eject > parameters.mass_minimum_internal) F_halos_reincorporation(&halos[i_halo],variables);
	    // Cooling of gas from halo onto central subhalo (or, in L-Galaxies mode, the most massive subhalo)
	    // Cooling occurs only if a central subhalo exists.
	    if (halos[i_halo].sub_central_sid != parameters.NO_DATA_INT) {
		i_sub_central=halos[i_halo].sub_central_sid;
		F_cooling_halo(&halos[i_halo],&subs[i_sub_central],dt_halo);
		// In l-galaxies mode the virial velocity of the subhalo may have changed, so need to reset that of the central galaxy also
		if (parameters.b_lgalaxies) {
		    i_gal_central=subs[i_sub_central].gal_central_sid;
		    if (i_gal_central != parameters.NO_DATA_INT) {
			gals[i_gal_central].v_vir=subs[i_sub_central].half_mass_virial_speed;
		    }
		}
	    }
	    
	    // Cooling of subhalos should probably be done on the galaxy timestep, in tandem with galaxy processing:
	    // in a sense the characteristic time of a subhalo is the same as that of a galaxy, as they represent galaxy halos.
	    // However, change in SFR over snap suggests that it would make very little difference and would slow down code.
	    // Loop over subhalos
	    for (i_sub=halos[i_halo].sub_start_sid; i_sub<halos[i_halo].sub_end_sid; i_sub++) {
		
		// Sanity check - subhalos should not be processed at this point
		if (subs[i_sub].b_done==true) {
		    printf("subhalo %d in graph %d already processed\n",subs[i_sub].sub_gid,halos[i_halo].graph_ID);
		    exit(1);
		}
		
		// Merge galaxies
		if (subs[i_sub].n_gal>1) {
		    // Initially assume instantaneous merging of galaxies in subhalos
		    subs[i_sub].gal_central_sid=F_mergers_merge_gals(&halos[i_halo],&subs[i_sub],gals,subs[i_sub].n_gal,variables);
		}
		
		// Cooling.  This also includes radio mode BH growth and feedback.
		// Not all subhalos may have hot gas
		if (subs[i_sub].mass_gas_hot > parameters.mass_minimum_internal) {
		    F_cooling_sub(&subs[i_sub],&gals[subs[i_sub].gal_central_sid],dt_halo);
		}

		subs[i_sub].n_dt+=1;
		if (subs[i_sub].n_dt==n_dt_halo) subs[i_sub].b_done=true;
	    }
	    // End loop over subhalos
	    
	    halos[i_halo].n_dt+=1;
	    if (halos[i_halo].n_dt==n_dt_halo) halos[i_halo].b_done=true;

	    // Loop over galaxy timestep
	    for (i_dt_gal=0; i_dt_gal<n_dt_gal; i_dt_gal++) {
	    
		// Loop over galaxies 
		for (i_gal=halos[i_halo].gal_start_sid; i_gal<halos[i_halo].gal_end_sid; i_gal++) {
		    if (!gals[i_gal].b_exists) continue; // Galaxies may have merged
		    gals[i_gal].SFR_dt = 0.; // This will fail to capture mergers (done above in subs loop) until we have a proper merger time for them.
		    if (gals[i_gal].mass_gas_cold > parameters.mass_minimum_internal) {
			mass_stars = F_SFF_gal_form_stars(&gals[i_gal],variables);
			// If subhalo does not exist, use halo as proxy.  This will work here as only need access to hot gas phase.
			sub_sid=gals[i_gal].sub_sid;
			if (sub_sid==parameters.NO_DATA_INT) {
			    F_SFF_orphan_SN_feedback(mass_stars,&gals[i_gal],&halos[i_halo]);
			} else {
			    F_SFF_gal_SN_feedback(mass_stars,&gals[i_gal],&subs[sub_sid],&halos[i_halo]);
			}
		    }
		}
		// End loop over galaxies
		
#ifdef SFH
		// Update SFH bins
		F_sfh_update_bins(gals,n_gal,i_dt_sfh);
		i_dt_sfh+=1;
#endif
		// Some SFR loggers
		if (i_dt_gal==0 && i_dt_halo==0) {
		    for (i_gal=0; i_gal<n_gal; i_gal++) {
			gals[i_gal].SFR_dt_start=gals[i_gal].SFR_dt; // This will contain the mergers
		    }
		}

	    }
	    //End loop over galaxy timestep
	
	}
	// End loop over halos
	
    }
    // End loop over halo timestep
    
    return;
}
