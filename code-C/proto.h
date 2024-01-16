/* Contains function prototypes.*/

// bh_agn.c
double F_BH_growth_quasar(double mass_ratio, double mass_gas_cold, double v_vir, double f_BH, double v_BH);
double F_BH_growth_rate_radio(double mass_gas_hot,double mass_BH, double f_BH);

// cooling.c
double F_cooling_get_metaldependent_cooling_rate(double log10_T, double log10_Z);
double F_cooling_SIS(double mass, double tau_dyn, double half_mass_radius, double mass_gas, double mass_metals_gas, double temp_start, double temp_end, double dt);
void F_cooling_halo(struct struct_halo *halo,struct struct_sub *sub, double dt);
void F_cooling_sub(struct struct_gal *gal, struct struct_sub *sub, double dt);

// star_formation_and_feedback.c
#ifdef SFH
double F_SFF_gal_form_stars(struct struct_gal *gal, double dt, double dt_snap, int i_bin_sfh);
#else
double F_SFF_gal_form_stars(struct struct_gal *gal, double dt, double dt_snap);
#endif
void F_SFF_gal_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_sub *sub, struct struct_halo *halo);
void F_SFF_orphan_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_halo *halo);
struct struct_SN_feedback F_SFF_SN_feedback_Hen15(double mass_stars, double v_vir, double mass_gas_cold);
double F_SFF_star_formation_unresolved(double mass_gas, double R_d, double v_vir, double dt);
