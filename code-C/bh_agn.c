/*
Functions related to black holes and agn activity
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//------------------------------------------------------------------------------------------------------

double F_BH_growth_quasar(double mass_ratio, double mass_gas_cold, double v_vir, double f_BH, double v_BH) {
    /*
   Calculates the mass driven onto the bh in a merger (quasar mode).

   Mergers of galaxies trigger merger of black holes and also a transfer of cold gas onto the black hole.  
   The mergers are handled in the calling routine; this function returns the amount of extra mass 
   transferred onto the BH.

   Arguments
   ---------
   mass_ratio : double
      The ratio of the baryonic mass of the smaller to the larger galaxy in the merger.
   mass_gas_cold : double
      The mass of cold gas in the merged galaxy
   v_vir : double
      The circular speed of the host (sub)halo of the galaxy.
   f_BH : double
      The normalisation parameter from Hen15 S23.
   v_BH : double
      The transition speed from Hen15 S23.

   Returns
   -------
   double
      The mass accreted onto the BH.
    */

    double BH_growth;

    if (mass_ratio >1.) {
	printf("F_BH_growth_quasar: mass_ratio >1.\n");
	exit(1);
    }
    BH_growth = f_BH * mass_ratio * mass_gas_cold / (1 + pow(v_BH/v_vir,2));

    return BH_growth;
}

//------------------------------------------------------------------------------------------------------

double F_BH_growth_rate_radio(double mass_gas_hot,double mass_BH, double f_BH) {
    /*
      Accretion from cooling gas in the hot atmosphere as in Hen15, S24.

      Arguments
      ---------
      mass_gas_hot : double
         The mass of hot gas in the subhalo.
      mass_BH : double
         The mass of the black hole in the central galaxy of the subhalo.
      f_BH : double
         The normalisation constant k_BH from Hen15 S24, including the appropriate mass normalisation factors.

      Returns
      -------
      double
         The rate of accretion onto the BH (later adjusted for maximum required to offset cooling).
    */

    double BH_growth_rate;

    BH_growth_rate = f_BH * mass_gas_hot * mass_BH;

    return BH_growth_rate;
}
