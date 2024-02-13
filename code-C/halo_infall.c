/* halo_infall.c

   Functions for accretion of primordial gas, and of reincorporation of ejected gas onto halos.

*/

#include "all_headers.h"

//------------------------------------------------------------------------------------------------------------

void F_halo_reincorporation(struct struct_halo *halo, struct struct_var variables) {
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
