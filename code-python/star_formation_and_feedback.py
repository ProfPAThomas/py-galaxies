import numpy as np
import astropy.units as u
"""
Functions to make stars from galaxies and to provide feedback of cold gas from SNR
"""

def F_gal_form_stars(gal,parameters):
   """
   Creates stars from gas in the cold gas disc.
   """
   sfr_model=parameters.sfr_model
   if sfr_model == "Unresolved":
      mass_stars_imf=F_star_formation_unresolved(gal['mass_gas_cold'],gal['radius_gas_cold'],gal['v_vir'], \
                                             parameters.dt,parameters.sfr_efficiency,parameters.c_sfr_Mcrit)
   else:
      raise valueError('sfr model '+sfr_model+' not yet implemented')
   # For now assume instantaneous recycling back into the cold gas
   # Then the mass stored in stars is that AFTER recycling, not the initial mass
   mass_stars=(1.-parameters.sfr_recycle_fraction)*mass_stars_imf

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   mass_metals_stars = mass_stars * (gal['mass_metals_gas_cold'] / gal['mass_gas_cold'])
   gal['mass_gas_cold'] -= mass_stars
   gal['mass_metals_gas_cold'] -= mass_metals_stars
   gal['mass_stars_disc'] += mass_stars
   gal['mass_metals_stars_disc'] += mass_metals_stars
   
   return mass_stars

def F_star_formation_unresolved(mass_gas,R_d,v_vir,dt,sfr_efficiency,c_sfr_Mcrit):
   """ 
   Implements star formation assuming gas disk unresolved, using the model from Hen15 (arXiv:1410.0365) S1.6.
   """

   t_dyn = R_d / v_vir
   sfr_Mcrit = c_sfr_Mcrit * v_vir * R_d
   mass_excess = mass_gas-sfr_Mcrit
   # L-Galaxies assumes that dt/t_dyn <= 1 here.
   # It seems more correct to take an exponential decline (requires an exp but saves a check to see if mass goes negative)
   if mass_gas>1.:
       print('mass_excess/mass_gas = {:.3g}'.format(mass_excess/mass_gas))
       print('dt, t_dyn, dt/t_dyn = {:.3g}, {:.3g}, {:.3g}'.format(dt,t_dyn,dt/t_dyn))
       print('R_d, v_vir, t_dyn = {:.3g}, {:.3g}, {:.3g}'.format(R_d,v_vir,t_dyn))
   if mass_excess>0.:
      mass_stars = mass_excess * (1.-np.exp(-sfr_efficiency*dt/t_dyn))
   else:
      mass_stars = 0.
   return mass_stars

def F_gal_SNR_feedback(mass_stars,gal,sub,halo,parameters):
   """
   Implements feedback from SNR after star formation.
   For now, assume prompt feedback.
   """
   # Note: input sub maybe a halo if no subhalo exists.
    
   snr_model=parameters.snr_model
   if snr_model == "Hen15":
      # Is it better to pass the reheating parameters that we need individually?
      mass_halo, mass_eject = F_SNR_feedback_Hen15(mass_stars,sub.half_mass_virial_speed,gal['mass_gas_cold'],parameters)
   else:
      raise valueError('snr model '+snr_model+' not yet implemented')
   mass_reheat = mass_halo + mass_eject

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   Zmet = gal['mass_metals_gas_cold'] / gal['mass_gas_cold']
   mass_metals_halo = mass_halo * Zmet
   mass_metals_eject = mass_eject * Zmet
   mass_metals_reheat = mass_reheat * Zmet
   gal['mass_gas_cold'] -= mass_reheat
   gal['mass_metals_gas_cold'] -= mass_metals_reheat
   sub.mass_gas_hot += mass_halo
   sub.mass_metals_gas_hot += mass_metals_halo
   # Ejected gas leaves the subhalo (but not the halo, if the subhalo is a halo)
   try:
        sub.sub_sid
        sub.mass_baryon -= mass_eject
   except:
        pass
   halo.mass_gas_eject += mass_eject
   halo.mass_metals_gas_eject += mass_metals_eject
   
   return None

def F_SNR_feedback_Hen15(mass_stars,v_vir,mass_gas_cold,parameters):
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
    