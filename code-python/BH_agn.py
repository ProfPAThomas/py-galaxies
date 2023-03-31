"""
Functions related to black holes and agn activity
"""
import numpy as np


def F_BH_growth_quasar(mass_ratio,mass_gas_cold,v_vir,f_BH,v_BH):
   """
   Mergers of galaxies trigger merger of black holes and also a transfer of cold gas onto the black hole.  
   The mergers are handled in the calling routine; this function returns the amount of extra mass 
   transferred onto the BH.
   gal1 is the larger, 'central' galaxy; gal2 is the satellite.
   """

   assert mass_ratio <=1
   BH_growth = f_BH * mass_ratio * mass_gas_cold / (1 + (v_BH/v_vir)**2)

   return BH_growth

def F_BH_growth_rate_radio(mass_gas_hot,mass_BH,f_BH):
   """
   Accretion from cooling gas in the hot atmosphere.
   Note that the normalisation constant here accounts for the mass normalisationin eq S24 of Hen15
   """

   BH_growth_rate = f_BH * mass_gas_hot * mass_BH

   return BH_growth_rate
