/* cooling.c

   Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy.
   The cooling table itself is generated in the python module cooling.py and stored in the header file cooling.h

*/

#include "all_headers.h"
#include "cooling.h"

//------------------------------------------------------------------------------------------------------------

double F_cooling_get_metaldependent_cooling_rate(double log10_T, double log10_Z) {
    /*
    Returns the cooling function, ie the cooling rate per unit density of electrons and ions.
    Assumes that the cooling function is tabulated in code units.

    Arguments
    ---------
    log10_T : double
       Temperature of cooling gas
    log10_Z : double
       Metallicities associated with cooling table entries.

    Returns
    -------
    double
       Value of cooling function in units such as to simplify cooling formula in code.
    */

    int i_T, i_Z;
    double fracT, fracZ;
    double log10_Lambda, log10_Lambda0, log10_Lambda1;

    // This should make i_T and i_T+1 bracket the input temperature.
    // We do NOT allow extrapolation.
    if (log10_T>=log10_T_table[n_T-2]) {
	i_T=n_T-2;
	log10_T=log10_T_table[n_T-1];
    }
    else {
	i_T = 0;
	while (log10_T>log10_T_table[i_T+1]) i_T++;
    }
    fracT=(log10_T-log10_T_table[i_T])/(log10_T_table[i_T+1]-log10_T_table[i_T]);
    if (fracT<0. || fracT>1.) {
	printf("fracT = %f\n",fracT);
	printf("log10_T = %f\n",log10_T);
	printf("i_T =%d\n",i_T);
	printf("log10_T_table[i_T-1] =%f\n",log10_T_table[i_T-1]);
	printf("log10_T_table[i_T] =%f\n",log10_T_table[i_T]);
	exit(EXIT_FAILURE);
    }
    // Ditto for metallicity
    if (log10_Z>=log10_Z_table[n_Z-2]) {
	i_Z=n_Z-2;
	log10_Z=log10_Z_table[n_Z-1];
    }
    else {
	i_Z = 0;
	while (log10_Z>log10_Z_table[i_Z+1]) i_Z++;
    }
    fracZ=(log10_Z-log10_Z_table[i_Z])/(log10_Z_table[i_Z+1]-log10_Z_table[i_Z]);
    if (fracZ<0. || fracZ>1.) {
	printf("fracZ = %f\n",fracZ);
	exit(EXIT_FAILURE);
    }

    // Interpolate in 2-d
    log10_Lambda0 = fracT*log10_Lambda_table[i_Z][i_T+1]+(1-fracT)*log10_Lambda_table[i_Z][i_T];
    log10_Lambda1 = fracT*log10_Lambda_table[i_Z+1][i_T+1]+(1-fracT)*log10_Lambda_table[i_Z+1][i_T];
    log10_Lambda = fracZ*log10_Lambda1+(1-fracZ)*log10_Lambda0;
    return pow(10.,log10_Lambda);
}

//------------------------------------------------------------------------------------------------------------

double F_cooling_SIS(double mass, double tau_dyn, double half_mass_radius, double mass_gas, double mass_metals_gas, double temp_start, double temp_end, double dt) {
    /*
    Implements the isothermal cooling model as used in L-Galaxies and many other SAMs.

    Can be used to determine the amount of cooled gas either in halos or subhalos.

    Arguments
    ---------
    mass : double
       The mass of the halo/subhalo.
    tau_dyn : double
       The dynamical time of the halo/subhalo.
    half_mass_radius : double
       The half_mass_radius of the halo/subhalo.
    mass_gas : double
       The mass of hot gas in the halo/subhalo.
    mass_metals_gas : double
       The mass of metals in the hot gas in the halo/subhalo.
    temp_start : double
       The starting temperature of the cooling gas.
    temp_end : double
       The final temperature of the cooling gas.
    dt : double
       The time for which the gas cools.
    cooling_table : obj : np.array[,]
       Cooling table in units such as to simplify cooling formula in code.

    Returns
    -------
    double
       Mass of cooled gas.
    */

    double dt_ratio, fg0, fg, log10_Z, Lambda, tau_cool, tau_ratio, teq_ratio;

    
    // Not sure if subhalo virial temperature can ever exceed that of the halo that it is in.
    // If it can, trap out before call to this subroutine, so raise error here.
    if (temp_end >= temp_start) {
        printf("cooling:F_cooling_SIS: mass, temp_start, temp_end = %f, %f, %f",mass,temp_start,temp_end);
        return mass_gas;
    }
    
    // Determine cooling function
    log10_Z = log10(mass_metals_gas/mass_gas);
    // Cooling rate per unit density of electrons & ions
    Lambda = F_cooling_get_metaldependent_cooling_rate(log10(temp_start),log10_Z);
    // Could save a little time in defining a conversion factor for half_mass_radius**3/mass in halos/subhalos.
    // (Because we execute the cooling every mini-step).
    tau_cool = pow(half_mass_radius,3)*(temp_start-temp_end)/(mass*Lambda);
    
    // Cooling at constant temperature, but allowing density to vary: see documentation.
    // The gas fraction here is relative to the halo mass (= 200 times the critical density in Millennium/SIS)
    fg0 = mass_gas/mass;
    dt_ratio=dt/tau_dyn;
    tau_ratio=tau_dyn*fg0/tau_cool;
    if (tau_ratio <=1.)
        fg=fg0/pow(1+0.5*sqrt(tau_ratio)*dt_ratio,2);
    else {
        teq_ratio=log(tau_ratio);
        if (dt_ratio <= teq_ratio)
            fg=fg0*exp(-dt_ratio);
        else
            fg=fg0/(tau_ratio*pow(1+0.5*(dt_ratio-teq_ratio),2));
    }
            
    return (fg0-fg)*mass;
}

//------------------------------------------------------------------------------------------------------------

void F_cooling_halo(struct struct_halo *halo,struct struct_sub *sub, double dt) {
    /*
    Cooling of halo onto subhalo.

    Arguments
    ---------
    halo : obj : struct struct_halo
       Properties of the halo
    sub : obj : struct struct_sub
       Properties of the subhalo
    dt : double
       The timestep over which cooling takes place.
    */

    double mass_cooled, mass_metals_cooled;
    
    if (parameters.b_lgalaxies) {
	/*
        Simple model to mimic that of L-Galaxies.  
        The hot gas in halos and subhalos is regarded as a single entity with temperature equal to that of the halo.
        So all we do here is give the hot gas to the subhalo, and
        set the virial speed and temperature of the subhalo to equal that of the halo.
        We will also need to update the virial speed of the galaxies upon return from this function.
        */
        (*sub).mass_gas_hot += (*halo).mass_gas_hot;
	(*sub).mass_metals_gas_hot += (*halo).mass_metals_gas_hot;
        (*sub).mass_baryon += (*halo).mass_gas_hot;
        (*halo).mass_gas_hot = parameters.mass_minimum_internal; // Just so it appears in the plots as a non-zero value.
	(*halo).mass_metals_gas_hot = parameters.mass_minimum_internal*parameters.base_metallicity;
	// Set subhalo properties to match that of the halo
        (*sub).temperature = (*halo).temperature;
	(*sub).half_mass_virial_speed = (*halo).half_mass_virial_speed;
    }
    else if (strcmp(parameters.cooling_model,"SIS")==0) {
        mass_cooled=F_cooling_SIS((*halo).mass,(*halo).tau_dyn,(*halo).half_mass_radius,(*halo).mass_gas_hot,
				  (*halo).mass_metals_gas_hot,(*halo).temperature,(*sub).temperature,dt);
        mass_metals_cooled  = (mass_cooled/(*halo).mass_gas_hot) * (*halo).mass_metals_gas_hot;
	(*halo).mass_gas_hot -= mass_cooled;
        (*halo).mass_metals_gas_hot -= mass_metals_cooled;
        (*sub).mass_gas_hot += mass_cooled;
        (*sub).mass_metals_gas_hot += mass_metals_cooled;
        (*sub).mass_baryon += mass_cooled;
    }
    else {
        printf("cooling.F_cooling_halo: cooling model %s not implemented.\n",parameters.cooling_model);
	exit(1);
    }
    return;
}

//------------------------------------------------------------------------------------------------------------

void F_cooling_sub(struct struct_sub *sub, struct struct_gal *gal, double dt) {
    /*
    Cooling of subhalo onto galaxy.
    Also sets the radius of the disc.

    Arguments
    ---------
    gal : obj : D_gal
       The central galaxy of the subhalo currently being processed.
    sub : obj: D_sub
       Properties of the subhalo. 
    dt : double
       The timestep over which cooling takes place.
    */

    double ang_mom_gas_cold, dm_BH, dm_BH_heat_max, dm_BH_max, dm_Edd, dm_metals_BH, efac;
    double gal_temperature, mass_cooled, mass_heated, mass_metals_cooled;
    
    // Angular momentum assuming exponential disc is 2vR_dM where R_d is the exponential disk radius
    ang_mom_gas_cold = 2. * (*gal).mass_gas_cold * (*gal).v_vir * (*gal).radius_gas_cold;

    // Mass cooling in the absence of AGN heating
    if ((*sub).mass_gas_hot <= parameters.mass_minimum_internal)
	return;
				
    if (strcmp(parameters.cooling_model,"SIS")==0) {
        gal_temperature=parameters.temperature_1e4K_internal; // Cool down to 1e4 K
        mass_cooled=F_cooling_SIS((*sub).mass,(*sub).tau_dyn,(*sub).half_mass_radius,(*sub).mass_gas_hot,
				  (*sub).mass_metals_gas_hot,(*sub).temperature,gal_temperature,dt);
    }
    else {
        printf("cooling.F_cooling_sub: cooling model %s not implemented.\n",parameters.cooling_model);
	exit(1);
    }

    // Radio mode growth of BH
    // Formula from Hen15 S24
    dm_BH_max = F_BH_growth_rate_radio((*sub).mass_gas_hot,(*gal).mass_BH,parameters.c_BH_r)*dt;
    // Eddington limit
    dm_Edd = parameters.c_BH_Edd * (*gal).mass_BH * dt;
    // Amount of growth needed to fully offset cooling
    efac = parameters.c_BH_mheat_r/pow((*sub).half_mass_virial_speed,2);
    dm_BH_heat_max = mass_cooled / (1+efac);
    // Hence actual BH accretion rate
    dm_BH = fmin(fmin(dm_BH_max,dm_Edd),dm_BH_heat_max);
    // Modify the amount of gas cooled
    mass_heated = dm_BH * efac;
    mass_cooled = mass_cooled-mass_heated;
    if (mass_cooled < 0.) {
        if (mass_cooled > -1e-10)
            mass_cooled = 0.;
        else {
            printf("F_cooling_sub: negative amount of gas cooled\n");
	    exit(1);
	}
    }
    mass_metals_cooled  = (mass_cooled/(*sub).mass_gas_hot) * (*sub).mass_metals_gas_hot;
    (*sub).mass_gas_hot -= mass_cooled;
    (*sub).mass_metals_gas_hot -= mass_metals_cooled;
    (*gal).mass_gas_cold += mass_cooled;
    (*gal).mass_metals_gas_cold += mass_metals_cooled;
    // Cooled gas will add to the baryon content of galaxies
    (*gal).mass_baryon += mass_cooled;
        
    // Disc radius (assuming exponential disc for cold gas and SIS for halo gas)
    // Accreted angular momentum for SIS is (1/2)RVM*lambda where lambda=parameters.halo_angular_momentum
    if ((*gal).mass_gas_cold > parameters.mass_minimum_internal) {
	ang_mom_gas_cold += mass_cooled * (*sub).half_mass_virial_speed * (*sub).half_mass_radius * parameters.halo_angular_momentum;
	    (*gal).radius_gas_cold = ang_mom_gas_cold / (2 * (*gal).mass_gas_cold * (*gal).v_vir);
	}
    else
        (*gal).radius_gas_cold = 0.; // Set to arbitrary value
        
    // We will transfer gas from the hot gas to the BH, even if it releases more energy than required to prevent cooling
    dm_metals_BH  = (dm_BH/(*sub).mass_gas_hot) * (*sub).mass_metals_gas_hot;
    (*sub).mass_gas_hot -= dm_BH;
    (*sub).mass_metals_gas_hot -= dm_metals_BH;
    (*gal).mass_BH += dm_BH;
    (*gal).mass_metals_BH += dm_metals_BH;
    (*gal).mass_baryon += dm_BH;
    return;
}
