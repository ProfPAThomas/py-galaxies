/* Contains function prototypes.*/

// bh_agn.c
double F_BH_growth_quasar(double mass_ratio, double mass_gas_cold, double v_vir, double f_BH, double v_BH);
double F_BH_growth_rate_radio(double mass_gas_hot,double mass_BH, double f_BH);

// cooling.c
double F_cooling_get_metaldependent_cooling_rate(double log10_T, double log10_Z);
double F_cooling_SIS(double mass, double tau_dyn, double half_mass_radius, double mass_gas, double mass_metals_gas, double temp_start, double temp_end, double dt);
void F_cooling_halo(struct struct_halo *halo, struct struct_sub *sub, double dt);
void F_cooling_sub(struct struct_sub *sub, struct struct_gal *gal, double dt);

// halos.c
void F_halos_accrete_primordial_gas(struct struct_halo *halo);
void F_halos_central_subhalo(struct struct_halo halos[], int n_halo, struct struct_sub subs[]);
void F_halos_reincorporation(struct struct_halo *halo, struct struct_var variables);
void F_halos_set_mass_baryon(struct struct_halo halos[], int n_halo, struct struct_sub subs[], struct struct_gal gals[], int n_dt_halo);

// mergers.c
// Use first form of call if you need to change entries in the struct variables, and the second if you do not.
// Will be interesting to see if the two differ in timing test.
// int F_mergers_merge_gals(struct struct_halo *halo,struct struct_sub *sub, struct struct_gal gals[],
// 			 int n_gal, struct struct_var *variables);
int F_mergers_merge_gals(struct struct_halo *halo,struct struct_sub *sub, struct struct_gal gals[],
			 int n_gal, struct struct_var variables);
// double F_mergers_starburst(double mass_ratio, struct struct_gal *gal_main, struct struct_var *variables);
double F_mergers_starburst(double mass_ratio, struct struct_gal *gal_main, struct struct_var variables);

// process_snap.c
void F_process_snap(struct struct_halo halos[], int n_halo, struct struct_sub subs[], int n_sub,
		    struct struct_gal gals[], int n_gal, struct struct_var variables);

// star_formation_and_feedback.c
// Use first form of call if you need to change entries in the struct variables, and the second if you do not.
// Will be interesting to see if the two differ in timing test.
//double F_SFF_gal_form_stars(struct struct_gal *gal, struct struct_var *variables);
double F_SFF_gal_form_stars(struct struct_gal *gal, struct struct_var variables);
void F_SFF_gal_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_sub *sub, struct struct_halo *halo);
void F_SFF_orphan_SN_feedback(double mass_stars, struct struct_gal *gal, struct struct_halo *halo);
struct struct_SN_feedback F_SFF_SN_feedback_Hen15(double mass_stars, double v_vir, double mass_gas_cold);
double F_SFF_star_formation_unresolved(double mass_gas, double R_d, double v_vir, double dt);

#ifdef SFH
void F_sfh_update_bins(struct struct_gal gals[], int n_gal, int i_dt);
#endif
