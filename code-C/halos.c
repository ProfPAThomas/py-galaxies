/*
halos.c - functions acting on halos
*/

#include "all_headers.h"

//----------------------------------------------------------------------------------------------------------------

void F_halos_accrete_primordial_gas(struct struct_halo *halo) {
    /*
      Updates the baryon content as calcuated previously in the halo set_mass_baryon method.
      Any excess baryons arrive in the form of base_metallicity hot gas.
    */
    
    double delta_mass;
    
    delta_mass=(*halo).mass_baryon_delta_dthalo;
    (*halo).mass_baryon+=delta_mass;
    (*halo).mass_gas_hot+=delta_mass;
    (*halo).mass_metals_gas_hot+=delta_mass*parameters.base_metallicity;
    return;
}

//----------------------------------------------------------------------------------------------------------------

void F_halos_central_subhalo(struct struct_halo halos[], int n_halo, struct struct_sub subs[]){
    /*
      Identify the central subhalo within each halo.
      For now, we will assume that there IS a central subhalo; later we may relax this assumption.
      The appropriate metric for defining distance from the centre of phase-space also needs exploring.
      Could do this one halo at a time and return the result, but probably more efficient to do all at once.
    */
    int i_halo, i_sub, i_sub_central, n_sub, sub_start_sid;
    double metric, metric2, metric2_min;

    for (i_halo=0; i_halo<n_halo; i_halo++) {
	i_sub_central=parameters.NO_DATA_INT;    
	n_sub=halos[i_halo].n_sub;
	if (n_sub>0) {
	    sub_start_sid=halos[i_halo].sub_start_sid;
	    if (parameters.b_lgalaxies) {
		// In L-galaxies mode the most massive subhalo is assigned as the central subhalo.
		metric=0.;
		for (i_sub=sub_start_sid; i_sub<sub_start_sid+n_sub; i_sub++) {
		    if (subs[i_sub].mass>metric) {
			metric=subs[i_sub].mass;
			i_sub_central=i_sub;
		    }
		}
	    } else {
		metric2=INFINITY;
		for (i_sub=sub_start_sid; i_sub<sub_start_sid+n_sub; i_sub++) {
		    metric2=(pow(subs[i_sub].pos[0]-halos[i_halo].pos[0],2)+pow(subs[i_sub].pos[1]-halos[i_halo].pos[1],2)+
			     pow(subs[i_sub].pos[2]-halos[i_halo].pos[2],2))/pow(halos[i_halo].rms_radius,2)+
			    (pow(subs[i_sub].vel[0]-halos[i_halo].vel[0],2)+pow(subs[i_sub].vel[1]-halos[i_halo].vel[1],2)+
			     pow(subs[i_sub].vel[2]-halos[i_halo].vel[2],2))/pow(halos[i_halo].rms_speed,2);
		    if (metric2<metric2_min) {
			metric2_min=metric2;
			i_sub_central=i_sub;
		    }
		}
	    }
	    halos[i_halo].sub_central_sid=i_sub_central;
	}
    }
    return;
}

//----------------------------------------------------------------------------------------------------------------
       
void F_halos_reincorporation(struct struct_halo *halo, struct struct_var variables) {
    /*
      Reincorporation of ejected gas.

      Currently just assumes Hen15 model.
      Might be better to pass dt_halo and c_reinc explicitly.
      Assumes that the minimum reincorporation timescale is the dynamical timescale, tau_dyn.

      Arguments
      ---------
      halo : obj : struct struct_halo
         Properties of the halo.
      variables : obj: struct struct_var
         Struct containing runtime variables.

      Returns
      -------
      None
    */

    double dt_halo, mass_reinc, mass_metals_reinc, t_reinc;
    
    dt_halo=variables.dt_halo;
    
    t_reinc = fmax(parameters.c_Hen15_reinc/(*halo).mass,(*halo).tau_dyn);
    mass_reinc = (*halo).mass_gas_eject * (1.-exp(-dt_halo/t_reinc));
    mass_metals_reinc = mass_reinc * ((*halo).mass_metals_gas_eject/(*halo).mass_gas_eject);
    (*halo).mass_gas_eject -= mass_reinc;
    (*halo).mass_metals_gas_eject -= mass_metals_reinc;
    (*halo).mass_gas_hot += mass_reinc;
    (*halo).mass_metals_gas_hot += mass_metals_reinc;
    return;
}

//------------------------------------------------------------------------------------------------------------

void F_halos_set_mass_baryon(struct struct_halo halos[], int n_halo, struct struct_sub subs[], struct struct_gal gals[], int n_dt_halo) {
    /*
      Calculates and updates the total baryonic mass of each halo, inclusive of subhalos and galaxies.
      Also determines the accretion required to bring the halo up to the universal mean baryon fraction, 
      or the sum of the baryon content from the progenitors, whichever is larger (so that baryons are not lost).

      Parameters
      ----------
      halos : obj : D_halo[n_halo]
        The halos in the snapshot
      n_halo : int
        The number of halos in this snapshot
      subs : obj: D_sub[n_sub]
        The subhalos in this snapshot
      gals : obj: D_gal[n_gal]
         The galaxies  contained within this halo, inclusive of subhalos
      n_dt_halo : int
         The number of halo timesteps in this snapshot interval.

      Returns
      -------
      None
    */

    int i_halo, i_orphan, i_sub, orphan_end_sid;
    double mass_baryon_delta;

    for (i_halo=0; i_halo<n_halo; i_halo++) {
	halos[i_halo].mass_baryon = halos[i_halo].mass_gas_hot + halos[i_halo].mass_gas_eject + halos[i_halo].mass_stars;
	for (i_sub=halos[i_halo].sub_start_sid; i_sub<halos[i_halo].sub_end_sid; i_sub++) {
	    halos[i_halo].mass_baryon += subs[i_sub].mass_baryon;
	}
        // The orphan galaxies are not included in the subhalo baryon count, so add them in here
	if (halos[i_halo].n_orphan >0) {
	    orphan_end_sid=halos[i_halo].orphan_start_sid+halos[i_halo].n_orphan;
	    for (i_orphan=halos[i_halo].orphan_start_sid; i_orphan<orphan_end_sid; i_orphan++) {
		halos[i_halo].mass_baryon += gals[i_sub].mass_gas_cold+gals[i_sub].mass_stars_bulge+gals[i_sub].mass_stars_disc;
	    }
	}
         
	// The amount of accretion needed per halo timestep to bring the baryon content up to the universal mean over a snapshot interval.
	mass_baryon_delta = parameters.baryon_fraction*fmax(halos[i_halo].mass,halos[i_halo].mass_from_progenitors) - halos[i_halo].mass_baryon;
	halos[i_halo].mass_baryon_delta_dthalo=fmax(0.,mass_baryon_delta)/(float)n_dt_halo;
    }
    return;
}
