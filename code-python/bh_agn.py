"""
Functions related to black holes and agn activity
"""
import numpy as np
from codetiming import Timer
from profiling import conditional_decorator

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_BH_growth_quasar',logger=None),True)
def F_BH_growth_quasar(mass_ratio,mass_gas_cold,v_vir,f_BH,v_BH):
   """
   Calculates the mass driven onto the bh in a merger (quasar mode).

   Mergers of galaxies trigger merger of black holes and also a transfer of cold gas onto the black hole.  
   The mergers are handled in the calling routine; this function returns the amount of extra mass 
   transferred onto the BH.

   Arguments
   ---------
   mass_ratio : float
      The ratio of the baryonic mass of the smaller to the larger galaxy in the merger.
   mass_gas_cold : float
      The mass of cold gas in the merged galaxy
   v_vir : float
      The circular speed of the host (sub)halo of the galaxy.
   f_BH : float
      The normalisation parameter from Hen15 S23.
   v_BH : float
      The transition speed from Hen15 S23.

   Returns
   -------
   float
      The mass accreted onto the BH.
   """

   assert mass_ratio <=1
   BH_growth = f_BH * mass_ratio * mass_gas_cold / (1 + (v_BH/v_vir)**2)

   return BH_growth

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_BH_growth_rate_radio',logger=None),True)
def F_BH_growth_rate_radio(mass_gas_hot,mass_BH,f_BH):
   """
   Accretion from cooling gas in the hot atmosphere as in Hen15, S24.

   Arguments
   ---------
   mass_gas_hot : float
      The mass of hot gas in the subhalo.
   mass_BH : float
      The mass of the black hole in the central galaxy of the subhalo.
   f_BH : float
      The normalisation constant k_BH from Hen15 S24, including the appropriate mass normalisation factors.

   Returns
   -------
   float
      The rate of accretion onto the BH (later adjusted for maximum required to offset cooling).
   """

   BH_growth_rate = f_BH * mass_gas_hot * mass_BH

   return BH_growth_rate
