/* star_formation_and_feedback.c

Functions to make stars from galaxies and to provide feedback of cold gas from SN.

*/

#include "all_headers.h"

struct struct_SN_feedback{double mass_reheat; double mass_eject_max;};

//-------------------------------------------------------------------------------------------------------------

void F_SFF_gal_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_sub *sub, struct struct_halo *halo){
    /*
    Implements feedback from SN after star formation.

    For now, assume prompt feedback.

    Arguments
    ---------
    mass_stars : float
       The initial mass of stars formed, before recycling.
    gal : obj : D_gal
       The galaxy currently being processed.
    sub : obj : D_sub
       The physical properties of the host subhalo of this galaxy.
    halo : obj : D_halo
       The physical properties of the host halo of this galaxy.
   
    Returns
    -------
    None
    */

    double mass_eject, mass_eject_max, mass_metals_eject, mass_metals_reheat, mass_reheat, Zmet;
    struct struct_SN_feedback SN_feedback;
    
    if (strcmp(parameters.snr_model,"Hen15")==0) {
	SN_feedback = F_SFF_SN_feedback_Hen15(mass_stars,(*sub).half_mass_virial_speed,(*gal).mass_gas_cold);
	mass_reheat=SN_feedback.mass_reheat;
	mass_eject_max=SN_feedback.mass_eject_max;
    } else {
        printf("star_formation_and_feedback.F_gal_SN_feedback: snr model %s not implemented.\n",parameters.snr_model);
	exit(1);
    }

    // Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
    if (mass_reheat > parameters.mass_minimum_internal) {
	Zmet = (*gal).mass_metals_gas_cold / (*gal).mass_gas_cold;
	mass_metals_reheat = mass_reheat * Zmet;
	(*gal).mass_gas_cold -= mass_reheat;
        (*gal).mass_metals_gas_cold -= mass_metals_reheat;
        (*sub).mass_gas_hot += mass_reheat;
	(*sub).mass_metals_gas_hot += mass_metals_reheat;
        (*gal).mass_baryon -= mass_reheat;
    }
    // The following assumes full mixing of metals before subsequent ejection
    mass_eject=fmin(mass_eject_max,(*sub).mass_gas_hot);
    if (mass_eject > parameters.mass_minimum_internal) {
	Zmet = (*sub).mass_metals_gas_hot / (*sub).mass_gas_hot;
	mass_metals_eject = mass_eject * Zmet;
	(*sub).mass_gas_hot -= mass_eject;
        (*sub).mass_metals_gas_hot -= mass_metals_eject;
        (*halo).mass_gas_eject += mass_eject;
        (*halo).mass_metals_gas_eject += mass_metals_eject;
    }
    // Ejected gas leaves the subhalo
    (*sub).mass_baryon -= mass_eject;
    return;
}

//-------------------------------------------------------------------------------------------------------------

void F_SFF_orphan_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_halo *halo){
    /*
    Implements feedback from SN after star formation for galaxies which are orphans (ie no subhalo).
    During code development, treat the halo as the subhalo: later we should probably think more carefully about
    what the differences would be.
    Note that the python version of the code made no distinction between F_gal_SN_feedback and F_orphan_SN_feedback,
    but in the C version of the code they have to be distinct (no overloading of functions).

    For now, assume prompt feedback.

    Arguments
    ---------
    mass_stars : float
       The initial mass of stars formed, before recycling.
    gal : obj : D_gal
       The galaxy currently being processed.
    halo : obj : D_halo
       The physical properties of the host halo of this galaxy.
   
    Returns
    -------
    None
    */

    double mass_eject, mass_eject_max, mass_metals_eject, mass_metals_reheat, mass_reheat, Zmet;
    struct struct_SN_feedback SN_feedback;
    
    if (strcmp(parameters.snr_model,"Hen15")==0) {
	SN_feedback = F_SFF_SN_feedback_Hen15(mass_stars,(*halo).half_mass_virial_speed,(*gal).mass_gas_cold);
	mass_reheat=SN_feedback.mass_reheat;
	mass_eject_max=SN_feedback.mass_eject_max;
    } else {
        printf("star_formation_and_feedback.F_gal_SN_feedback: snr model %s not implemented.\n",parameters.snr_model);
	exit(1);
    }

    // Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
    if (mass_reheat > parameters.mass_minimum_internal) {
	Zmet = (*gal).mass_metals_gas_cold / (*gal).mass_gas_cold;
	mass_metals_reheat = mass_reheat * Zmet;
	(*gal).mass_gas_cold -= mass_reheat;
        (*gal).mass_metals_gas_cold -= mass_metals_reheat;
        (*halo).mass_gas_hot += mass_reheat;
	(*halo).mass_metals_gas_hot += mass_metals_reheat;
        (*gal).mass_baryon -= mass_reheat;
    }
    // The following assumes full mixing of metals before subsequent ejection
    mass_eject=fmin(mass_eject_max,(*halo).mass_gas_hot);
    if (mass_eject > parameters.mass_minimum_internal) {
	Zmet = (*halo).mass_metals_gas_hot / (*halo).mass_gas_hot;
	mass_metals_eject = mass_eject * Zmet;
	(*halo).mass_gas_hot -= mass_eject;
        (*halo).mass_metals_gas_hot -= mass_metals_eject;
        (*halo).mass_gas_eject += mass_eject;
        (*halo).mass_metals_gas_eject += mass_metals_eject;
    }
    return;
}

//------------------------------------------------------------------------------------------------------------

struct struct_SN_feedback F_SFF_SN_feedback_Hen15(double mass_stars, double v_vir, double mass_gas_cold) {
    /*
    Evaluates the mass returned to the corona and ejected from the halo in the Hen15 model.
	
    Arguments
    ---------
    mass_stars : double
       The mass of stars formed, prior to recycling.
    v_vir : double
       The circular speed of the host (sub)halo of the galaxy.
    mass_gas_cold : double
       The amount of cold gas in the galaxy prior to feedback.

    Returns
    -------
    SN_feedback : obj : struct_SN_feedback
       Struct containing:
          mass_reheat : float
             The mass reheated into the hot gas (coronal) phase.
          mass_eject_max : float
       The (maximum) mass ejected from the host halo.
    */

    double DE_feedback, DM_eject_max, DM_feedback, DM_reheat, DM_reheat_max, epsilon_feedback, mu_reheat_max;
    // Eq S16 & S17 SN feedback energy (note: have absorbed minus sign into ratio by swapping num. & denom.)
    epsilon_feedback = parameters.Hen15_eta*(0.5+pow(parameters.Hen15_v_eject_internal/v_vir,parameters.Hen15_beta2));
    DE_feedback = parameters.c_Hen15_S16*epsilon_feedback*mass_stars;
    // Amount of mass that this energy can raise to the virial temperature (well, actually, 0.5*v_vir^2)
    DM_feedback = DE_feedback / (0.5 * pow(v_vir,2));
    // Eq S18 & S19 SN reheated mass (note: have absorbed minus sign into ratio by swapping num. & denom.)
    mu_reheat_max = parameters.Hen15_epsilon*(0.5+pow(parameters.Hen15_v_reheat_internal/v_vir,parameters.Hen15_beta1));
    DM_reheat_max = mu_reheat_max*mass_stars;
    // Can't reheat more gas than exists, or that we have energy for
    DM_reheat = fmin(fmin(DM_feedback,DM_reheat_max),mass_gas_cold);
    // Use excess energy to eject gas from the (sub)halo
    DM_eject_max = DM_feedback-DM_reheat;

    struct struct_SN_feedback SN_feedback;
    SN_feedback.mass_reheat=DM_reheat;
    SN_feedback.mass_eject_max=DM_eject_max;
    return SN_feedback;
}

//------------------------------------------------------------------------------------------------------------

double F_SFF_star_formation_unresolved(double mass_gas, double R_d, double v_vir, double dt) {
    /*
      Implements star formation assuming gas disk unresolved, using the model from Hen15 (arXiv:1410.0365) S1.6 S14/15.

      Arguments
      ---------
      mass_gas : double
         The mass of the cold gas disc.
      R_d : double
         The cold gas disc exponential scale length.
      v_vir : double
         The circular speed of the halo that is associated with the galaxy
      dt : double
         The (galaxy) timestep.

      Returns
      -------
      double
         The mass of stars formed (before recycling).
    */

    double mass_excess, mass_stars, sfr_Mcrit, t_dyn;
    t_dyn = R_d / v_vir;
    sfr_Mcrit = parameters.c_sfr_Mcrit * v_vir * R_d;
    mass_excess = mass_gas-sfr_Mcrit;
    // L-Galaxies assumes that dt/t_dyn <= 1 here.
    // It seems more correct to take an exponential decline (requires an exp but saves a check to see if mass goes negative)
    if (mass_excess>0.) {
	mass_stars = mass_excess * (1.-exp(-parameters.sfr_efficiency*dt/t_dyn));
    } else {
	mass_stars = 0.;
    }
    return mass_stars;
}

//------------------------------------------------------------------------------------------------------------

/*
  Use first form of call if you need to change entries in the struct variables, and the second if you do not.
  Entries are accessed by (*variables).name in the first instance and variables.name in the latter.
  Will be interesting to see if the two differ in timing test.
*/
//double F_SFF_gal_form_stars(struct struct_gal *gal, struct struct_var *variables) {
double F_SFF_gal_form_stars(struct struct_gal *gal, struct struct_var variables) {
    /*
      Creates stars from gas in the cold gas disc.

      Arguments
      ---------
      gal : obj : D_gal
	 The galaxy that represents the merger product.
      dt : double
         The galaxy timestep.
      i_bin_sfh : int
         [Only if -DSFH] The SFH bin currently in use.

      Returns
      -------
      double
	 The mass of stars formed (before recycling).
    */

    char* sfr_model;
    double gal_mass_gas_cold, gal_mass_metals_gas_cold, gal_mass_stars_disc, gal_mass_metals_stars_disc;
    double gal_radius_gas_cold, gal_radius_stars_disc, gal_v_vir;
    double ang_mom_stars_disc, mass_stars, mass_stars_imf, mass_metals_stars;
    
    // Extract the variables that we need
    double dt, dt_snap;
    //dt=(*variables).dt_gal;
    dt=variables.dt_gal;
    //dt_snap=(*variables).dt_snap;
    dt_snap=variables.dt_snap;
#ifdef SFH
    int i_bin_sfh;
    //i_bin_sfh=(*variables).i_bin_sfh;
    i_bin_sfh=variables.i_bin_sfh;
#endif

    gal_mass_gas_cold=(*gal).mass_gas_cold;
    gal_mass_metals_gas_cold=(*gal).mass_metals_gas_cold;
    gal_mass_stars_disc=(*gal).mass_stars_disc;
    gal_mass_metals_stars_disc=(*gal).mass_metals_stars_disc;
    gal_radius_gas_cold=(*gal).radius_gas_cold;
    gal_radius_stars_disc=(*gal).radius_stars_disc;
    gal_v_vir=(*gal).v_vir;
   
    // We will need this later to set the new disc scale length
    ang_mom_stars_disc = 2. * gal_mass_stars_disc * gal_v_vir * gal_radius_stars_disc;

    // Determine the mass of stars formed before (mass_stars_imf) and after (mass_stars) prompt recycling
    sfr_model=parameters.sfr_model;
    if (strcmp(sfr_model,"Unresolved")==0) {
	mass_stars_imf=F_SFF_star_formation_unresolved(gal_mass_gas_cold,gal_radius_gas_cold,gal_v_vir, dt);
    } else {
	printf("sfr model %s not yet implemented\n",sfr_model);
	exit(1);
    }
    if (mass_stars_imf < parameters.mass_minimum_internal) return 0.;

    // Record star formation rates
    (*gal).SFR_dt += mass_stars_imf/dt;        // zeroed at start of timestep
    (*gal).SFR_snap += mass_stars_imf/dt_snap; // This one is cumulative over the snapshot
   
    // For now assume instantaneous recycling back into the cold gas
    // Then the mass stored in stars is that AFTER recycling, not the initial mass
    mass_stars=(1.-parameters.sfr_recycle_fraction)*mass_stars_imf;

    // Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
    mass_metals_stars = mass_stars * ((*gal).mass_metals_gas_cold / (*gal).mass_gas_cold);
    gal_mass_gas_cold -= mass_stars;
    gal_mass_metals_gas_cold -= mass_metals_stars;
    gal_mass_stars_disc += mass_stars;
    gal_mass_metals_stars_disc += mass_metals_stars;
#ifdef SFH
    (*gal).mass_stars_disc_sfh[i_bin_sfh-1] += mass_stars;
    (*gal).mass_metals_stars_disc_sfh[i_bin_sfh-1] += mass_metals_stars;
#endif

    // Now we enrich our surroundings with prompt metals returned from the newly formed stars
    // To begin with we will assume that all metal enrichment is prompt and everything goes to cold gas.
    // (The return of gas mass to the cold gas has been handled by the recycling fraction above.)
    gal_mass_metals_gas_cold += parameters.sfr_yield * mass_stars_imf;
    if (gal_mass_metals_gas_cold>gal_mass_gas_cold*0.2) {
	printf("Warning, high Z for cold gas: mass, Z = %.3g, %.3g\n", gal_mass_gas_cold,gal_mass_metals_gas_cold/gal_mass_gas_cold);
    }
    // We want to model the size of the stellar disc.  We do this using angular momentum (assuming exponential in shape)
    ang_mom_stars_disc += mass_stars * gal_v_vir * gal_radius_gas_cold;
    gal_radius_stars_disc = ang_mom_stars_disc / (2 * gal_mass_stars_disc * gal_v_vir);

    // Now finally update everything; could do this as we went along but danger of confusing initial and final values
    (*gal).mass_gas_cold=gal_mass_gas_cold;
    (*gal).mass_metals_gas_cold=gal_mass_metals_gas_cold;
    (*gal).mass_stars_disc=gal_mass_stars_disc;
    (*gal).mass_metals_stars_disc=gal_mass_metals_stars_disc;
    (*gal).radius_gas_cold=gal_radius_gas_cold;
    (*gal).radius_stars_disc=gal_radius_stars_disc;

    return mass_stars_imf;
}
