"""
Functions to make stars from galaxies and to provide feedback of cold gas from SNR.
"""

import numpy as np
import astropy.units as u
from codetiming import Timer
from profiling import conditional_decorator

import commons
b_profile_cpu=commons.load('b_profile_cpu')

#-------------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_gal_form_stars',logger=None),b_profile_cpu)
def F_gal_form_stars(gal,parameters):
   """
   Creates stars from gas in the cold gas disc.

   Arguments
   ---------
   gal : obj : D_gal
      The galaxy that represents the merger product.
   parameters : obj : C_parameters
      Instance of class containing global parameters.

   Returns
   -------
   float
      The mass of stars formed (before recycling).
   """

   b_SFH=parameters.b_SFH
   
   # The timestep
   dt_gal=commons.load('dt_gal')
   if b_SFH: i_bin_sfh=commons.load('i_bin_sfh')
   
   # We will need this later to set the new disc scale length
   ang_mom_stars_disc = 2. * gal['mass_stars_disc'] * gal['v_vir'] * gal['radius_stars_disc']

   # Determine the mass of stars formed before (mass_stars_imf) and after (mass_stars) prompt recycling
   sfr_model=parameters.sfr_model
   if sfr_model == "Unresolved":
      mass_stars_imf=F_star_formation_unresolved(gal['mass_gas_cold'],gal['radius_gas_cold'],gal['v_vir'], \
                                             dt_gal,parameters.sfr_efficiency,parameters.c_sfr_Mcrit)
   else:
      raise valueError('sfr model '+sfr_model+' not yet implemented')
   if mass_stars_imf < parameters.mass_minimum_internal: return 0.

   # Record star formation rates
   gal['SFR_dt'] += mass_stars_imf/dt_gal                    # zeroed at start of timestep
   gal['SFR_snap'] += mass_stars_imf/commons.load('dt_snap') # This one is cumulative over the snapshot
   
   # For now assume instantaneous recycling back into the cold gas
   # Then the mass stored in stars is that AFTER recycling, not the initial mass
   mass_stars=(1.-parameters.sfr_recycle_fraction)*mass_stars_imf

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   mass_metals_stars = mass_stars * (gal['mass_metals_gas_cold'] / gal['mass_gas_cold'])
   gal['mass_gas_cold'] -= mass_stars
   gal['mass_metals_gas_cold'] -= mass_metals_stars
   gal['mass_stars_disc'] += mass_stars
   gal['mass_metals_stars_disc'] += mass_metals_stars
   if b_SFH:
         gal['mass_stars_disc_sfh'][i_bin_sfh-1] += mass_stars
         gal['mass_metals_stars_disc_sfh'][i_bin_sfh-1] += mass_metals_stars

   # Now we enrich our surroundings with prompt metals returned from the newly formed stars
   # To begin with we will assume that all metal enrichment is prompt and everything goes to cold gas.
   # (The return of gas mass to the cold gas has been handled by the recycling fraction above.)
   gal['mass_metals_gas_cold'] += parameters.sfr_yield * mass_stars_imf
   if gal['mass_metals_gas_cold']>gal['mass_gas_cold']*0.2:
      print('Warning, high Z for cold gas: mass, Z =',gal['mass_gas_cold'],gal['mass_metals_gas_cold']/gal['mass_gas_cold'])

   # We want to model the size of the stellar disc.  We do this using angular momentum (assuming exponential in shape)
   ang_mom_stars_disc += mass_stars * gal['v_vir'] * gal['radius_gas_cold']
   gal['radius_stars_disc'] = ang_mom_stars_disc / (2 * gal['mass_stars_disc'] * gal['v_vir'])      

   return mass_stars_imf

#-------------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_star_formation_unresolved',logger=None),b_profile_cpu)
def F_star_formation_unresolved(mass_gas,R_d,v_vir,dt,sfr_efficiency,c_sfr_Mcrit):
   """ 
   Implements star formation assuming gas disk unresolved, using the model from Hen15 (arXiv:1410.0365) S1.6.

   Arguments
   ---------
   mass_gas : float
      The mass of the cold gas disc.
   R_d : float
      The cold gas disc exponential scale length.
   v_vir : float
      The circular speed of the halo that is associated with the galaxy
   dt : float
      The (galaxy) timestep.
   sfr_efficiency : float
      The SFR efficiency from Hen15 S14.
   c_sfr_Mcrit : float
      The mass threshold for star formation from Hen15 S14/S15.

   Returns
   -------
   float
      The mass of stars formed (before recycling).
   """

   t_dyn = R_d / v_vir
   sfr_Mcrit = c_sfr_Mcrit * v_vir * R_d
   mass_excess = mass_gas-sfr_Mcrit
   # L-Galaxies assumes that dt/t_dyn <= 1 here.
   # It seems more correct to take an exponential decline (requires an exp but saves a check to see if mass goes negative)
#    if mass_gas>10.:
#        print('mass_excess/mass_gas = {:.3g}'.format(mass_excess/mass_gas))
#        print('dt, t_dyn, dt/t_dyn = {:.3g}, {:.3g}, {:.3g}'.format(dt,t_dyn,dt/t_dyn))
#        print('R_d, v_vir, t_dyn, SFE = {:.3g}, {:.3g}, {:.3g}, {:.3g}'.format(R_d,v_vir,t_dyn,(1.-np.exp(-sfr_efficiency*dt/t_dyn))))
   if mass_excess>0.:
      mass_stars = mass_excess * (1.-np.exp(-sfr_efficiency*dt/t_dyn))
   else:
      mass_stars = 0.
   return mass_stars

#-------------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_gal_SNR_feedback',logger=None),b_profile_cpu)
def F_gal_SNR_feedback(mass_stars,gal,sub,halo,parameters):
   """
   Implements feedback from SNR after star formation.

   For now, assume prompt feedback.

   Arguments
   ---------
   mass_stars : float
      The initial mass of stars formed, before recycling.
   gal : obj : D_gal
      The galaxy currently being processed.
   sub : obj : C_sub
      The host subhalo of this galaxy (or host halo if no host subhalo exists).
   halo : obj : C_halo
      The host halo of this galaxy.
   parameters : obj : C_parameters
      Instance of class containing global parameters.
   
   Returns
   -------
   None
   """
   # Note: input sub may be a halo if no subhalo exists.
    
   snr_model=parameters.snr_model
   if snr_model == "Hen15":
      # Is it better to pass the reheating parameters that we need individually?
      mass_halo, mass_eject = F_SNR_feedback_Hen15(mass_stars,sub.half_mass_virial_speed,gal['mass_gas_cold'],parameters)
   else:
      raise valueError('snr model '+snr_model+' not yet implemented')
   mass_reheat = mass_halo + mass_eject
   #print('mass_halo, mass_eject =',mass_halo,mass_eject)

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   Zmet = gal['mass_metals_gas_cold'] / gal['mass_gas_cold']
   mass_metals_halo = mass_halo * Zmet
   mass_metals_eject = mass_eject * Zmet
   mass_metals_reheat = mass_reheat * Zmet
   gal['mass_gas_cold'] -= mass_reheat
   gal['mass_metals_gas_cold'] -= mass_metals_reheat
   sub.mass_gas_hot += mass_halo
   sub.mass_metals_gas_hot += mass_metals_halo
   gal['mass_baryon'] -= mass_reheat
   # Ejected gas leaves the subhalo (but not the halo, if the subhalo is a halo)
   try:
        sub.sub_sid
        sub.mass_baryon -= mass_eject
   except:
        pass
   halo.mass_gas_eject += mass_eject
   halo.mass_metals_gas_eject += mass_metals_eject
   
   return None

#-------------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_SNR_feedback_Hen15',logger=None),b_profile_cpu)
def F_SNR_feedback_Hen15(mass_stars,v_vir,mass_gas_cold,parameters):
   """
   Evaluates the mass returned to the corona and ejected from the halo in the Hen15 model.

   Arguments
   ---------
   mass_stars : float
      The mass of stars formed, prior to recycling.
   v_vir : float
      The circular speed of the host (sub)halo of the galaxy.
   mass_gas_cold : float
      The amount of cold gas in the galaxy prior to feedback.
   parameters : obj : C_parameters
      Instance of class containing global parameters.

   Returns
   -------
   float
      The mass reheated into the hot gas (coronal) phase.
   float
      The mass ejected from the host halo.
   """
   
   # Eq S16 & S17 SNR feedback energy (note: have absorbed minus sign into ratio by swapping num. & denom.)
   epsilon_halo = parameters.Hen15_eta*(0.5+(parameters.Hen15_v_eject_internal/v_vir)**parameters.Hen15_beta2)
   DE_snr = parameters.c_Hen15_S16*epsilon_halo*mass_stars
   # Amount of mass that this energy can raise to the virial temperature (well, actually, 0.5*v_vir^2)
   DM_snr = DE_snr / (0.5 * v_vir**2)
   # Eq S18 & S19 SNR reheated mass (note: have absorbed minus sign into ratio by swapping num. & denom.)
   epsilon_disk = parameters.Hen15_epsilon*(0.5+(parameters.Hen15_v_reheat_internal/v_vir)**parameters.Hen15_beta1)
   DM_reheat_max = epsilon_disk*mass_stars
   # Can't reheat more gas than exists, or that we have energy for
   DM_reheat = min(DM_snr,DM_reheat_max,mass_gas_cold)
   # Use excess energy to eject gas from the (sub)halo
   DM_eject = min(DM_snr-DM_reheat,DM_reheat)
   DM_halo = DM_reheat - DM_eject
    
   return DM_halo, DM_eject
 
