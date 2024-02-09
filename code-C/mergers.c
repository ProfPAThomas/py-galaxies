/* mergers.c

Functions related to galaxy mergers.

*/

#include "all_headers.h"

/* Currently dt and dt_snap not used (only needed to update SFRs).
   Could eliminate the need for multiple function parameter lists, and excessive passing of variables,
   if combined everything currently in common into a struct - 
   should do that as it would be in the spirit of making the code much simpler to understand.
*/
/* Could avoid a return value below if included gal_central_sid in struct_sub.
   Not doing that at the moment as it is a book-keeping parameter rather than a physical property.
   However, this is maybe a moot point as struct_gal contains some book-keeping parameters.
*/

/*
  Use first form of call if you need to change entries in the struct variables, and the second if you do not.
  Entries are accessed by (*variables).name in the first instance and variables.name in the latter.
  Will be interesting to see if the two differ in timing test.
*/
// int F_mergers_merge_gals(struct struct_halo *halo, struct struct_sub *sub, struct struct_gal gals[],
// 			 int n_gal, struct struct_var *variables) {
int F_mergers_merge_gals(struct struct_halo *halo, struct struct_sub *sub, struct struct_gal gals[],
			 int n_gal, struct struct_var variables) {
    /*
    Merges galaxies within subhalo.

    This currently merges all galaxies.  
    In the future, it should check to see whether galaxies should be merged yet or not.

    Arguments
    ---------
    halo : obj : C_halo
       The host halo of the subhalo currently being processed.
    sub : obj : C_sub
       The subhalo currently being processed.
    parameters : obj : C_parameters
       Instance of class containing global parameters
    gals : obj : D_gal[n_gal]
       The galaxies in the subhalo currently being processed.
    n_gal : int
       The number of galaxies in the subhalo.

    Returns
    -------
    gal_central_sid : int
       The sid of the galaxy at the centre of the subhalo
    */

    int gal_central_sid, i_gal, i_main;
    double ang_mom_gas_cold, dm_BH, dm_metals_BH, gal_main_mass_baryon;
    double mass_gas_cold_before_starburst, mass_bulge_main, mass_bulge_sat, mass_ratio, mass_starburst;
    double max_mass, PE_bulge, radius_half_main, radius_half_sat, radius_sum;
    struct struct_gal *gal_main, *gal_sat;
    
    // Find the most massive galaxy: we will take this to be the one onto which everything accretes
    // Note: this is a pain in C; would it be better to pass in from python? i_main=np.argmax(gals['mass_baryon'])
    i_main=0;
    max_mass=(gals[0]).mass_baryon;
    for (i_gal=1; i_gal<n_gal; i_gal++) {
	if ((gals[i_gal]).mass_baryon > max_mass) {
	    i_main=i_gal;
	    max_mass=(gals[i_gal]).mass_baryon;
	}
    }
    gal_main=&gals[i_main];

    // The central galaxy mass may change as a result of mergers, so use the initial one to separate minor/major
    gal_main_mass_baryon=(*gal_main).mass_baryon;
    // At this point rename the central galaxy to be the most massive one.
    // Note that the need to return this value to the python calling routine could be another reason to identify
    // the main galaxy within the python code instead.
    gal_central_sid = (*gal_main).gal_sid;

    // If no galaxies have any mass then can quit early
    if (gal_main_mass_baryon <= parameters.mass_minimum_internal) {
	// Newly created galaxies sharing a subhalo: nothing to do except set b_exists status
	for (i_gal=0; i_gal<n_gal; i_gal++) {
	    (gals[i_gal]).b_exists=false;
	}
	(*gal_main).b_exists=true;
	return gal_central_sid;
    }
    
    // Loop over each of the other galaxies in turn, triggering the mergers
    for (i_gal=0; i_gal<n_gal; i_gal++) {
        gal_sat = &gals[i_gal];

	// Don't merge with yourself or non-existent galaxies
        if ((i_gal==i_main) || ((*gal_sat).b_exists==false)) continue;  
	
        // If satellite galaxy has no mass then can skip
        if ((*gal_sat).mass_baryon < parameters.mass_minimum_internal) {
	    (*gal_sat).b_exists=false;
	    continue;
	}
        
        // Store quantities needed to determine the final gas disc and stellar bulge sizes
	mass_ratio = (*gal_sat).mass_baryon/gal_main_mass_baryon;
        ang_mom_gas_cold = 2. * ((*gal_main).mass_gas_cold*(*gal_main).v_vir*(*gal_main).radius_gas_cold +
                                 (*gal_sat).mass_gas_cold*(*gal_sat).v_vir*(*gal_sat).radius_gas_cold);
        // Mass-weighted half mass radii (for exp disc 1.68*R_d; for Jaffe R)
        mass_bulge_sat = (*gal_sat).mass_stars_disc+(*gal_sat).mass_stars_bulge;
        // Not sure how there can be no stars in the satellite, but there can be
        if (mass_bulge_sat > parameters.mass_minimum_internal) {
	    radius_half_sat = (1.68*(*gal_sat).mass_stars_disc*(*gal_sat).radius_stars_disc +
			       (*gal_sat).mass_stars_bulge*(*gal_sat).radius_stars_bulge) / mass_bulge_sat;
	} else {
	    radius_half_sat = 0.;
	}
        // Minor/major mergers differ in contributions to size calculation.  Equation S34 of Hen15.
        if (mass_ratio < parameters.major_merger_fraction) {
            // Minor merger.  Only bulge of main ends up in final bulge (excepting starburst contribution, added later)
            mass_bulge_main = (*gal_main).mass_stars_bulge;
            radius_half_main = (*gal_main).radius_stars_bulge;
        } else {
            // Major merger.  All stars in main galaxy end up in bulge
            mass_bulge_main = (*gal_main).mass_stars_disc+(*gal_main).mass_stars_bulge;
            if (mass_bulge_main >  parameters.mass_minimum_internal) {
                radius_half_main = (1.68*(*gal_main).mass_stars_disc*(*gal_main).radius_stars_disc +
				    (*gal_main).mass_stars_bulge*(*gal_main).radius_stars_bulge) / mass_bulge_main;
	    } else {
                radius_half_main = 0.;
	    }
	}
        PE_bulge = 0.;
        if (radius_half_main > parameters.length_minimum_internal) {
            PE_bulge += pow(mass_bulge_main,2)/radius_half_main;
	}
        if (radius_half_sat > parameters.length_minimum_internal) {
            PE_bulge += pow(mass_bulge_sat,2)/radius_half_sat;
	}
        // The following is the orbital contribution
        radius_sum = radius_half_main +radius_half_sat;
        if (radius_sum > parameters.length_minimum_internal) {
            PE_bulge += (mass_bulge_main+mass_bulge_sat)/radius_sum;
	}
                
        // Transfer mass onto central gal.  
        // Note that the disc of satellite goes into bulge of central galaxy: if there is a major merger
        // the disc of the main galaxy will get transferred to the bulge also.
        (*gal_main).mass_gas_cold += (*gal_sat).mass_gas_cold;
        (*gal_main).mass_metals_gas_cold += (*gal_sat).mass_metals_gas_cold;
        (*gal_main).mass_stars_bulge += (*gal_sat).mass_stars_bulge + (*gal_sat).mass_stars_disc;
        (*gal_main).mass_metals_stars_bulge += (*gal_sat).mass_metals_stars_bulge + (*gal_sat).mass_metals_stars_disc;
#ifdef SFH
	// Would be better to amend gals.h to include n_SFH rather than calculating it here.
	int i_SFH, n_SFH;
	n_SFH=sizeof((*gal_main).mass_stars_bulge_sfh)/sizeof((*gal_main).mass_stars_bulge_sfh[0]);
	for (i_SFH=0; i_SFH<n_SFH; i_SFH++) {
	    (*gal_main).mass_stars_bulge_sfh[i_SFH] += (*gal_sat).mass_stars_bulge_sfh[i_SFH] +
		(*gal_sat).mass_stars_disc_sfh[i_SFH];
            (*gal_main).mass_metals_stars_bulge_sfh[i_SFH] += (*gal_sat).mass_metals_stars_bulge_sfh[i_SFH] +
                (*gal_sat).mass_metals_stars_disc_sfh[i_SFH];
	}
#endif
	(*gal_main).mass_BH += (*gal_sat).mass_BH;
        (*gal_main).mass_metals_BH += (*gal_sat).mass_metals_BH;
        (*gal_main).mass_baryon += (*gal_sat).mass_baryon;
        // Not strictly necessary as will nullify galaxy, but tidier/safer this way
        (*gal_sat).mass_gas_cold=0.;
        (*gal_sat).mass_metals_gas_cold=0.;
        (*gal_sat).mass_stars_bulge=0.;
        (*gal_sat).mass_metals_stars_bulge=0.;
        (*gal_sat).mass_stars_disc=0.;
        (*gal_sat).mass_metals_stars_disc=0.;
#ifdef SFH
	for (i_SFH=0; i_SFH<n_SFH; i_SFH++) {
            (*gal_sat).mass_stars_bulge_sfh[i_SFH]=0.;
            (*gal_sat).mass_metals_stars_bulge_sfh[i_SFH]=0.;
            (*gal_sat).mass_stars_disc_sfh[i_SFH]=0.;
            (*gal_sat).mass_metals_stars_disc_sfh[i_SFH]=0.;
	}
#endif
        (*gal_sat).mass_BH=0.;
        (*gal_sat).mass_metals_BH=0.;
        (*gal_sat).mass_baryon=0.;
        (*gal_sat).b_exists=false;

        // Perform starburst
        if (((*gal_main).mass_gas_cold>=parameters.mass_minimum_internal) &&
	    (mass_ratio >= parameters.major_merger_fraction)) { 
            // Performs starburst and transfers stellar disc to bulge.
            // Note that the starburst will likely be compact but in the first instance we assume
            // that it follows the profile of the disc, so the angular momentum is reduced in proportion
            mass_gas_cold_before_starburst=(*gal_main).mass_gas_cold;
	    mass_starburst = F_mergers_starburst(mass_ratio,gal_main,variables);
            // Feedback associated with starburst
            F_SFF_gal_SN_feedback(mass_starburst,gal_main,sub,halo);
	    ang_mom_gas_cold *= (*gal_main).mass_gas_cold/mass_gas_cold_before_starburst;
	}
        
        // BH growth: could be absorbed into above but cleaner to do it here
        // There is no BH feedback associated with mergers
        if ((*gal_main).mass_gas_cold>=parameters.mass_minimum_internal) {
	    dm_BH = F_BH_growth_quasar(mass_ratio,(*gal_main).mass_gas_cold,(*gal_main).v_vir,
				       parameters.BH_f_q,parameters.BH_v_q_internal);
            dm_metals_BH = dm_BH * (*gal_main).mass_metals_gas_cold/(*gal_main).mass_gas_cold;
            (*gal_main).mass_BH += dm_BH;
            (*gal_main).mass_metals_BH += dm_metals_BH;
            (*gal_main).mass_gas_cold -= dm_BH;
            (*gal_main).mass_metals_gas_cold -= dm_metals_BH;
	}
            
        // Set new sizes for gas disc and stellar bulge
        // Note that the starburst is assumed compact and transfers ZERO angular momentum
        if ((*gal_main).mass_gas_cold>parameters.mass_minimum_internal) {
            (*gal_main).radius_gas_cold=ang_mom_gas_cold/(2. * (*gal_main).mass_gas_cold * (*gal_main).v_vir);
	}
        if ((*gal_main).mass_stars_bulge>parameters.mass_minimum_internal) {
            if (PE_bulge>0) {
                (*gal_main).radius_stars_bulge=pow((*gal_main).mass_stars_bulge,2)/PE_bulge;
	    } else {
                (*gal_main).radius_stars_bulge=0.;
            }
	}
    }	
    return gal_central_sid;
}

/*
  Use first form of call if you need to change entries in the struct variables, and the second if you do not.
  Entries are accessed by (*variables).name in the first instance and variables.name in the latter.
  Will be interesting to see if the two differ in a timing test.
*/
//double F_mergers_starburst(double mass_ratio, struct struct_gal *gal, struct struct_var *variables) {
double F_mergers_starburst(double mass_ratio, struct struct_gal *gal, struct struct_var variables) {
    /*
      Major mergers of galaxies trigger a burst of star formation that goes into the bulge of the remnant.

      The starburst will likely be compact but in the first instance we take it to follow the profile of the disc.

      Arguments
      ---------
      mass_ratio : float
         The ratio of the baryonic mass of the smaller galaxy to the larger one.
      gal : obj : D_gal
         The galaxy that represents the merger product.
      variables : obj : struct struct_var
         variables structure containing counters and timing information

      Returns
      -------
      double
         The mass of stars formed in the starburst (before recycling).
    */

    double mass_stars, mass_stars_imf, mass_metals_stars;

    // Extract the variables that we need
    //Uncomment the following if needed for SFR below
    //double dt, dt_snap;
    //dt=variables.dt_gal;
    //dt_snap=variables.dt_snap;
#ifdef SFH
    int i_bin_sfh;
    // i_bin_sfh=(*variables).i_bin_sfh;
    i_bin_sfh=variables.i_bin_sfh;
#endif

    if (mass_ratio < parameters.major_merger_fraction) {
	printf("F_mergers_starburst: mass_ratio < parameters.major_merger_fraction\n");
	exit(1);
    }
    if (mass_ratio >1.) {
	printf("F_mergers_starburst: mass_ratio >1.\n");
	exit(1);
    }
    mass_stars_imf = parameters.merger_f_burst * pow(mass_ratio,parameters.merger_beta_burst) * (*gal).mass_gas_cold;

    if (mass_stars_imf < parameters.mass_minimum_internal) return 0.;

    // Record star formation rates
    // Not sure why this is commented out, although SFR for merger should involve sub-grid modelling
    // (*gal).SFR_dt += mass_stars_imf/dt;
    // (*gal).SFR_snap += mass_stars_imf/dt_snap;
   
    // For now assume instantaneous recycling back into the cold gas
    // Then the mass stored in stars is that AFTER recycling, not the initial mass
    mass_stars=(1.-parameters.sfr_recycle_fraction)*mass_stars_imf;

    // Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
    mass_metals_stars = mass_stars * ((*gal).mass_metals_gas_cold / (*gal).mass_gas_cold);
    (*gal).mass_gas_cold -= mass_stars;
    (*gal).mass_metals_gas_cold -= mass_metals_stars;
    (*gal).mass_stars_bulge += mass_stars;
    (*gal).mass_metals_stars_bulge += mass_metals_stars;
#ifdef SFH
    (*gal).mass_stars_bulge_sfh[i_bin_sfh-1] += mass_stars;
    (*gal).mass_metals_stars_bulge_sfh[i_bin_sfh-1] += mass_metals_stars;
#endif
    // Now we enrich our surroundings with prompt metals returned from the newly formed stars
    // To begin with we will assume that all metal enrichment is prompt and everything goes to cold gas.
    // (The return of gas mass to the cold gas has been handled by the recycling fraction above.)
    // Subsequent feedback will then transfer some of those metals to the hot gas.
    // Later versions of L-Galaxies allowed some fraction of the metals to go straight to the hot gas phase
    (*gal).mass_metals_gas_cold += parameters.sfr_yield * mass_stars_imf;
    if ((*gal).mass_metals_gas_cold>(*gal).mass_gas_cold*0.2) {
        printf("Warning, high Z for cold gas: mass, Z = %.3g, %.3g\n",
	       (*gal).mass_gas_cold,(*gal).mass_metals_gas_cold/(*gal).mass_gas_cold);
    }

    return mass_stars_imf;
}

