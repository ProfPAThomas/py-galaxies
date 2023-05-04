"""
Functions related to galaxy mergers.
"""

import numpy as np
import astropy.units as u
from codetiming import Timer
from profiling import conditional_decorator
from bh_agn import F_BH_growth_quasar
from star_formation_and_feedback import F_gal_SNR_feedback

# Merge gals in subhalos

@conditional_decorator(Timer(name='F_merge_gals',logger=None),True)
def F_merge_gals(halo,sub,gals,parameters):
    """
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
    parameters : obj : C_parameters
       Instance of class containing global parameters.

    Returns
    -------
    None
    """
    # Find the most massive galaxy: we will take this to be the one onto which everything accretes
    i_main=np.argmax(gals['mass_baryon'])
    gal_main=gals[i_main]
    # The central galaxy mass may change as a result of mergers, so use the initial one to separate minor/major
    gal_main_mass_baryon=gal_main['mass_baryon']
    # At this point rename the central galaxy to be the most massive one
    sub.gal_central_sid = gal_main['gal_sid']
    # Loop over each of the other galaxies in turn, triggering the mergers
    # (I originally did this using np.sum, but better to check for starbursts one at a time)
    for i_gal in range(sub.n_gal):
        gal_sat = gals[i_gal]
        if i_gal == i_main or gal_sat['b_exists']==False: continue  # Don't merge with yourself!

        # If no galaxies have any mass then can quit early
        if gal_main_mass_baryon > parameters.mass_minimum_internal:
           mass_ratio = gal_sat['mass_baryon']/gal_main_mass_baryon
        else:
           # Newly created galaxies sharing a subhalo: nothing to do except set b_exists status
           gals[:]['b_exists']=False
           gals[i_main]['b_exists']=True
           return None

        # If satellite galaxy has no mass then can skip
        if gal_sat['mass_baryon'] < parameters.mass_minimum_internal:
           gal_sat['b_exists']=False
           continue
        
        # Store quantities needed to determine the final gas disc and stellar bulge sizes
        ang_mom_gas_cold = 2. * (gal_main['mass_gas_cold']*gal_main['v_vir']*gal_main['radius_gas_cold'] + \
                                 gal_sat['mass_gas_cold'] * gal_sat['v_vir'] * gal_sat['radius_gas_cold'])
        # Mass-weighted half mass radii (for exp disc 1.68*R_d; for Jaffe R)
        mass_bulge_sat = gal_sat['mass_stars_disc']+gal_sat['mass_stars_bulge']
        # Not sure how there can be no stars in the satellite, but there can be
        if mass_bulge_sat > parameters.mass_minimum_internal:
           radius_half_sat = (1.68*gal_sat['mass_stars_disc']*gal_sat['radius_stars_disc'] + \
                                gal_sat['mass_stars_bulge']*gal_sat['radius_stars_bulge']) / \
                          mass_bulge_sat
        else:
           radius_half_sat = 0.
        # Minor/major mergers differ in contributions to size calculation.  Equation S34 of Hen15.
        if mass_ratio < parameters.major_merger_fraction:
            # Minor merger.  Only bulge of main ends up in final bulge (excepting starburst contribution, added later)
            mass_bulge_main = gal_main['mass_stars_bulge']
            radius_half_main = gal_main['radius_stars_bulge']
        else:
            # Major merger.  All stars in main galaxy end up in bulge
            mass_bulge_main = gal_main['mass_stars_disc'] +gal_main['mass_stars_bulge']
            radius_half_main = (1.68*gal_main['mass_stars_disc']*gal_main['radius_stars_disc'] + \
                                gal_main['mass_stars_bulge']*gal_main['radius_stars_bulge']) / \
                               mass_bulge_main
        PE_bulge = 0.
        if radius_half_main > parameters.length_minimum_internal:
            PE_bulge += mass_bulge_main**2/radius_half_main
        if radius_half_sat > parameters.length_minimum_internal:
            PE_bulge += mass_bulge_sat**2/radius_half_sat
        # The following is the orbital contribution
        radius_sum = radius_half_main +radius_half_sat
        if radius_sum > parameters.length_minimum_internal:
            PE_bulge += (mass_bulge_main+mass_bulge_sat)/radius_sum
                
        # Transfer mass onto central gal.  
        # Note that only disc of satellite goes into bulge of central galaxy: if there is a major merger
        # the disc of the main galaxy will get transferred to the bulge also.
        gal_main['mass_gas_cold'] += gal_sat['mass_gas_cold']
        gal_main['mass_metals_gas_cold'] += gal_sat['mass_metals_gas_cold']
        gal_main['mass_stars_bulge'] += gal_sat['mass_stars_bulge'] + gal_sat['mass_stars_disc']
        gal_main['mass_metals_stars_bulge'] += gal_sat['mass_metals_stars_bulge'] + gal_sat['mass_metals_stars_disc']
        gal_main['mass_BH'] += gal_sat['mass_BH']
        gal_main['mass_metals_BH'] += gal_sat['mass_metals_BH']
        gal_main['mass_baryon'] += gal_sat['mass_baryon']
        # Not strictly necessary as will nullify galaxy, but tidier/safer this way
        gal_sat['mass_gas_cold']=0.
        gal_sat['mass_metals_gas_cold']=0.
        gal_sat['mass_stars_bulge']=0.
        gal_sat['mass_metals_stars_bulge']=0.
        gal_sat['mass_stars_disc']=0.
        gal_sat['mass_metals_stars_disc']=0.
        gal_sat['mass_BH']=0.
        gal_sat['mass_metals_BH']=0.
        gal_sat['mass_baryon']=0.
        gal_sat['b_exists']=False

        # Perform starburst
        if mass_ratio >= parameters.major_merger_fraction: 
            # Performs starburst and transfers stellar disc to bulge
            mass_starburst = F_starburst(mass_ratio,gal_main,parameters)
            # Note that the starburst is assumed compact and transfers ZERO angular momentum
            # Feedback associated with starburst
            F_gal_SNR_feedback(mass_starburst,gal_main,sub,halo,parameters)
        
        # BH growth: could be absorbed into above but cleaner to do it here
        # There is no BH feedback associated with mergers
        if gal_main['mass_gas_cold']>=parameters.mass_minimum_internal:
            dm_BH = F_BH_growth_quasar(mass_ratio,gal_main['mass_gas_cold'],gal_main['v_vir'],
                                parameters.BH_f_q,parameters.BH_v_q_internal)
            dm_metals_BH = dm_BH * gal_main['mass_metals_gas_cold']/gal_main['mass_gas_cold']
            gal_main['mass_BH'] += dm_BH
            gal_main['mass_metals_BH'] += dm_metals_BH
            gal_main['mass_gas_cold'] -= dm_BH
            gal_main['mass_metals_gas_cold'] += dm_metals_BH
            
        # Set new sizes for gas disc and stellar bulge
        # Note that the starburst is assumed compact and transfers ZERO angular momentum
        if gal_main['mass_gas_cold']>parameters.mass_minimum_internal:
            gal_main['radius_gas_cold']=ang_mom_gas_cold/(2. * gal_main['mass_gas_cold'] * gal_main['v_vir'])
        if gal_main['mass_stars_bulge']>parameters.mass_minimum_internal:
            gal_main['radius_stars_bulge']=gal_main['mass_stars_bulge']**2/PE_bulge
            
    return None

@conditional_decorator(Timer(name='F_BH_starburst',logger=None),True)
def F_starburst(mass_ratio,gal,parameters):
   """
   Major mergers of galaxies trigger a burst of star formation that goes into the bulge of the remnant.

   The starburst is assumed compact and transfers zero angular momentum.

   Arguments
   ---------
   mass_ratio : float
      The ratio of the baryonic mass of the smaller galaxy to the larger one.
   gal : obj : D_gal
      The galaxy that represents the merger product.

   Returns
   -------
   float
      The mass of stars formed in the starburst (before recycling).
   """

   assert parameters.major_merger_fraction <= mass_ratio <=1, str(parameters.major_merger_fraction)+','+str(mass_ratio)
   mass_stars_imf = parameters.merger_f_burst * mass_ratio**parameters.merger_beta_burst * gal['mass_gas_cold']

   if mass_stars_imf < parameters.mass_minimum_internal: return 0.

   # Record star formation rates
   gal['SFR_dt'] += mass_stars_imf/parameters.dt_gal        # zeroed at start of timestep
   gal['SFR_snap'] += mass_stars_imf/parameters.dt_snap     # This one is cumulative over the snapshot
   
   # For now assume instantaneous recycling back into the cold gas
   # Then the mass stored in stars is that AFTER recycling, not the initial mass
   mass_stars=(1.-parameters.sfr_recycle_fraction)*mass_stars_imf

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   mass_metals_stars = mass_stars * (gal['mass_metals_gas_cold'] / gal['mass_gas_cold'])
   gal['mass_gas_cold'] -= mass_stars
   gal['mass_metals_gas_cold'] -= mass_metals_stars
   gal['mass_stars_bulge'] += mass_stars
   gal['mass_metals_stars_bulge'] += mass_metals_stars

   # Now we enrich our surroundings with prompt metals returned from the newly formed stars
   # To begin with we will assume that all metal enrichment is prompt and everything goes to cold gas.
   # (The return of gas mass to the cold gas has been handled by the recycling fraction above.)
   # Subsequent feedback will then transfer some of those metals to the hot gas.
   # Later versions of L-Galaxies allowed some fraction of the metals to go straight to the hot gas phase
   gal['mass_metals_gas_cold'] += parameters.sfr_yield * mass_stars_imf
   if parameters.b_debug and gal['mass_metals_gas_cold']>gal['mass_gas_cold']*0.2:
        print('Warning, high Z for cold gas: mass,Z',
              gal['mass_gas_cold'],gal['mass_metals_gas_cold']/gal['mass_gas_cold'])

   return mass_stars_imf
