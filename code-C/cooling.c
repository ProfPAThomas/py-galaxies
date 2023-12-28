/* cooling.c

   Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy.
   The cooling table itself is generated in the python module cooling.py and stored in the header file cooling.h

*/

#include "cooling.h"
#include "math.h"
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
