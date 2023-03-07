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
      mass_stars=F_star_formation_unresolved(gal['mass_gas_cold'],gal['radius_gas_cold'],gal['v_vir'], \
                                             parameters.dt,parameters.sfr_efficiency,parameters.c_sfr_Mcrit)
   else:
      raise valueError('sfr model '+sfr_model+' not yet implemented')

   # Need generic routine (as in L-Galaxies) to do these mass transfers, especially when get many components
   mass_metals_stars = mass_stars * (gal['mass_metals_gas_cold'] / gal['mass_gas_cold'])
   gal['mass_gas_cold'] -= mass_stars
   gal['mass_metals_gas_cold'] -= mass_metals_stars
   gal['mass_stars_disc'] += mass_stars
   gal['mass_metals_stars_disc'] + mass_metals_stars
   
   return None

def F_star_formation_unresolved(mass_gas,R_d,v_vir,dt,sfr_efficiency,c_sfr_Mcrit):
   """ 
   Implements star formation assuming gas disk unresolved, using the model from Hen15 (arXiv:1410.0365) S1.6.
   """

   t_dyn = R_d / v_vir
   if t_dyn>1e10: print('mass_gas, R_d, v_vir =',mass_gas, R_d, v_vir)
   sfr_Mcrit = c_sfr_Mcrit * v_vir * R_d
   mass_excess = mass_gas-sfr_Mcrit
   # L-Galaxies assumes that dt/t_dyn <= 1 here.
   # It seems more correct to take an exponential decline (requires an exp but saves a check to see if mass goes negative)
   #print('mass_excess =',mass_excess)
   #print('dt, t_dyn =',dt,t_dyn)
   if mass_excess>0.:
      mass_stars = mass_excess * (1.-np.exp(-sfr_efficiency*dt/t_dyn))
   else:
      mass_stars = 0.
   return mass_stars

