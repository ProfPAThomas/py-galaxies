/* cooling.c

   Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy.
   The cooling table itself is generated in the python module cooling.py and stored in the header file cooling.h

*/

#include "cooling.h"
#include "parameters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double F_get_metaldependent_cooling_rate(double log10_T, double log10_Z) {
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
    Lambda = F_get_metaldependent_cooling_rate(log10(temp_start),log10_Z);
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
