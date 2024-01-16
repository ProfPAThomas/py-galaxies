"""
A class to generate and store the cooling tables.
Other cooling routines are now ported to cooling.c
"""

import numpy as np
import astropy.units as u
   
class C_cooling:
    """
    Class to generate and store the cooling tables.

    Attributes
    ----------
    log10_T_table : obj : np.array[]
       Temperatures associated with cooling table entries.
    log10_Z_table : obj : np.array[]
       Metallicities associated with cooling table entries.
    log10_Lambda_table : np.array[,]
       Cooling table in units such as to simplify cooling formula in code.
    """
    def __init__(self,parameters):
        """
        Reads in the cooling function and tabulates it.
        Normalises it such that tau_cool=T/Lambda in code units.
        Also, the temperature lookup is adjusted to code units.
        Metallicities are always assumed to be absolute (ie not relative to solar).
        """
        import astropy.constants as c
        import astropy.units as u
        data=np.load(parameters.cooling_function_file)
        self.log10_T_table=(data['log10_T']+np.log10(u.K/parameters.units_temperature_internal).si.value)  # Converted to code units
        self.log10_Z_table=data['log10_Z']
        # The first correction term converts to code units, and the second provides a dimensionless factor to reduce the cooling time
        # formula to dt=T/Lambda
        self.log10_Lambda_table=data['log10_Lambda'] + np.log10(
            (u.erg*u.cm**3 / (u.s*parameters.units_lambda_internal)).si.value) - \
            np.log10(parameters.c_cooling)
        # Fudge to test cooling
        # self.log10_Lambda_table -= 0.3
        # print('log10_lambda_table.shape =',self.log10_Lambda_table.shape)
        # print('log10_Lambda_table =',self.log10_Lambda_table)

