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

    def F_create_header_file(self):
        """
        Writes out the cooling tables as const C arrays.
     
        A more obvious way to do this would be to write out the tables as binary files, then read
        them back in as const arrays in the C cooling routine the first time that it is called:
        that would be quicker and more accurate.  I am doing it this way instead because:
        * It makes the C-code cleaner (but this routine less so).
        * We don't need high accuracy.
        * I don't think that the speed difference will be very great.
     
        Note that a trailing comma seems to be permitted in C array initialisation, which simplifies things.
     
        Attributes
        ----------
        """
        f=open('code-C/cooling.h','w')
        f.write('/* Cooling tables (fixed throughout run). */\n\n')
     
        log10_T_table=self.log10_T_table
        n_T=len(log10_T_table)
        f.write('#define n_T '+str(n_T)+'\n')
        f.write('const double log10_T_table['+str(n_T)+'] = {')
        for i_T in range(n_T):
           f.write(str(log10_T_table[i_T])+', ')
        f.write('};\n\n')
     
        log10_Z_table=self.log10_Z_table
        n_Z=len(log10_Z_table)
        f.write('#define n_Z '+str(n_Z)+'\n')
        f.write('const double log10_Z_table['+str(n_Z)+'] = {')
        # Negative infinity will cause issues, so set first value to a very small number
        f.write('-100., ')
        for i_Z in range(1,n_Z):
           f.write(str(log10_Z_table[i_Z])+', ')
        f.write('};\n\n')
     
        log10_Lambda_table=self.log10_Lambda_table
        f.write('const double log10_Lambda_table['+str(n_Z)+']['+str(n_T)+'] = {\n')
        for i_Z in range(n_Z):
           f.write('{')
           for i_T in range(n_T):
              f.write(str(log10_Lambda_table[i_Z,i_T])+', ')
           f.write('},\n')
        f.write('};\n\n')   
     
        f.close()
        return None
   
