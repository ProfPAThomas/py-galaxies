"""
Miscellaneous helper routines that do not find a place elsewhere.
"""

import astropy.units as u
import numpy as np

import commons

def F_set_dt(parameters):
   """
   Sets the number and size of the timesteps.

   Attributes
   ----------
   parameters : obj : C_parameters
       Instance of class containing global parameters
   """
   snap_table=parameters.snap_table
   n_snap=parameters.n_snap
   # Convert time to code units
   t_snap=(snap_table['time_in_years'] * u.yr / parameters.units_time_internal).value
   parameters.t_snap=t_snap
   
   # The snapshot intervals in code units
   # Note that we are stepping between this snapshot and the previous one.
   dt_snap=np.full(n_snap,np.NAN,float)
   dt_snap[0]=t_snap[0]  # Assume that bin starts at time = 0.
   dt_snap[1:]=t_snap[1:]-t_snap[:-1]
   parameters.dt_snap=dt_snap
   
   # The number and size of the halo timesteps per snapshot.
   n_dt_halo=np.full(n_snap,parameters.NO_DATA_INT,dtype=np.int32)
   dt_halo=np.full(n_snap,np.NAN,float)
   # The 0.999999 factor prevents creating an extra step when the ratio is integral.
   n_dt_halo=(0.999999*dt_snap/parameters.timestep_halo_internal+1).astype(int)
   dt_halo=dt_snap/n_dt_halo
   parameters.n_dt_halo=n_dt_halo
   parameters.dt_halo=dt_halo

   # Number of galaxy timesteps per halo timestep.
   # This could vary between snapshots as the size of the halo steps adjusts to the snapshot interval.
   n_dt_gal=np.full(n_snap,parameters.NO_DATA_INT,dtype=np.int32)
   dt_gal=np.full(n_snap,np.NAN,float)   
   # The 0.999999 factor prevents creating an extra step when the ratio is integral.
   n_dt_gal=(0.999999*dt_halo/parameters.timestep_gal_internal+1).astype(int)
   dt_gal=dt_halo/n_dt_gal
   parameters.n_dt_gal=n_dt_gal
   parameters.dt_gal=dt_gal

   return None

def F_create_cooling_header_file(cooling_table):
   """
   Writes out the cooling tables as static C arrays.

   A more obvious way to do this would be to write out the tables as binary files, then read
   them back in as static arrays in the C cooling routine the first time that it is called:
   that would be quicker and more accurate.  I am doing it this way instead because:
   * It makes the C-code cleaner (but this routine less so).
   * We don't need high accuracy.
   * I don't think that the speed difference will be very great.

   Note that a trailing comma seems to be permitted in C array initialisation, which simplifies things.

   Attributes
   ----------
   parameters : obj : C_cooling
       Instance of class containing cooling tables
   """
   f=open('code-C/cooling.h','w')
   f.write('/* Cooling tables (fixed throughout run). */\n\n')

   log10_T_table=cooling_table.log10_T_table
   n_T=len(log10_T_table)
   f.write('#define n_T '+str(n_T)+'\n')
   f.write('static double log10_T_table['+str(n_T)+'] = {')
   for i_T in range(n_T):
      f.write(str(log10_T_table[i_T])+', ')
   f.write('};\n\n')

   log10_Z_table=cooling_table.log10_Z_table
   n_Z=len(log10_Z_table)
   f.write('#define n_Z '+str(n_Z)+'\n')
   f.write('static double log10_Z_table['+str(n_Z)+'] = {')
   # Negative infinity will cause issues, so set first value to a very small number
   f.write('-100., ')
   for i_Z in range(1,n_Z):
      f.write(str(log10_Z_table[i_Z])+', ')
   f.write('};\n\n')

   log10_Lambda_table=cooling_table.log10_Lambda_table
   f.write('static double log10_Lambda_table['+str(n_Z)+']['+str(n_T)+'] = {\n')
   for i_Z in range(n_Z):
      f.write('{')
      for i_T in range(n_T):
         f.write(str(log10_Lambda_table[i_Z,i_T])+', ')
      f.write('},\n')
   f.write('};\n\n')   

   f.close()
   return None

def F_create_galaxy_struct_header_file(D_gal):
   """
   Creates a C struct definition that matches the galaxy dtype.
   Writes out to code/gals.h

   Attributes
   ----------
   D_gal : obj : numpy.dtype
       Numpy dtype used in the galaxy structured array.
   """
   f=open('code-C/gals.h','w')
   f.write('/* Contains struct definition for galaxies. */\n\n#include <stdbool.h>\n\nstruct struct_gal {\n')
   for key in D_gal.fields.keys():
      var_type=str(D_gal[key])
      # Do the awkward arrays first
      if ',))' in var_type:
         var_type=var_type.split(',')[1]
         var_type=var_type.strip(' (')
         f.write('    double '+key+'['+var_type+'];\n')
      # Now the simple variables
      elif 'bool' in var_type:
         f.write('    bool '+key+';\n')         
      elif 'int' in var_type:
         f.write('    int '+key+';\n')
      elif 'float' in var_type:
         f.write('    double '+key+';\n')
      else:
         f.write('    '+var_type+' '+key+';\n')
   f.write('}; \n')
   f.close()
   return None
   

def F_create_halo_struct_header_file(D_halo):
   """
   Creates a C struct definition that matches the subhalo dtype.
   Writes out to code/subs.h

   Attributes
   ----------
   D_halo : obj : numpy.dtype
       Numpy dtype used in the halo structured array.
   """
   f=open('code-C/halos.h','w')
   f.write('/* Contains struct definition for properties of halos. */\n\n#include <stdbool.h>\n\nstruct struct_halo {\n')
   for key in D_halo.fields.keys():
      var_type=str(D_halo[key])
      # First deal with the awkward arrays
      if '(3,)' in var_type:
         f.write('    double '+key+'[3];\n')
      # then process the simple types
      elif 'bool' in var_type:
         f.write('    bool '+key+';\n')         
      elif 'int' in var_type:
         f.write('    int '+key+';\n')
      elif 'float' in var_type:
         f.write('    double '+key+';\n')
      else:
         f.write('    '+var_type+' '+key+';\n')
   f.write('}; \n')
   f.close()
   return None
   

def F_create_parameters_header_file(parameters):
   """
   Writes out all the attributes of parameters to code/parameters.h.
   Not clear to me whether I should declare this as a static struct (memory preserved over runtime of program).

   Attributes
   ----------
   parameters : obj : C_parameters
       Instance of class containing global parameters
   """
   attributes=[a for a in dir(parameters) if not a.startswith('__') and not callable(getattr(parameters, a))]
   f=open('code-C/parameters.h','w')
   f.write('/* Runtime parameters (fixed throughout run). */\n\
\n\
#include <stdbool.h>\n\
\n\
struct struct_param {\n\
')
   for a in attributes:
      value=eval('parameters.'+a)
      a_type=str(type(value)).split('\'')[1]
      print(a,a_type,value)
      # First the parts that we wish to ignore
      if 'astropy' in a_type:
         continue
      elif 'ndarray' in a_type:
         continue
      elif 'C_sfh' in a_type:
         continue
      # Now the bits that we want to extract
      elif a_type == 'bool':
         if value==True:
            f.write('    bool '+a+';\n')
         else:
            f.write('    bool '+a+';\n')
      elif 'dict' in a_type:
         continue
      elif 'int' in a_type:
         f.write('    int '+a+';\n')
      elif 'float' in a_type:
         f.write('    float '+a+';\n')
      elif a_type == 'str':
         f.write('    char* '+a+';\n')
      else:
         f.write('    '+a_type+' '+a+';\n')
   f.write('};\n')
   f.write('static struct struct_param parameters = {\n')
   for a in attributes:
      value=eval('parameters.'+a)
      a_type=str(type(value)).split('\'')[1]
      # First the parts that we wish to ignore
      if 'astropy' in a_type:
         continue
      elif 'ndarray' in a_type:
         continue
      elif 'C_sfh' in a_type:
         continue
      # Now the bits that we want to extract
      elif a_type == 'bool':
         if value==True:
            f.write('    .'+a+'=true,\n')
         else:
            f.write('    .'+a+'=false,\n')
      elif 'dict' in a_type:
         continue
      elif a_type == 'str':
         f.write('    .'+a+'="'+value+'",\n')
      else:
         f.write('    .'+a+'='+str(value)+',\n')
   f.write('};\n')
   
   f.close()
   return None


def F_create_sub_struct_header_file(D_sub):
   """
   Creates a C struct definition that matches the subhalo dtype.
   Writes out to code/subs.h

   Attributes
   ----------
   D_sub : obj : numpy.dtype
       Numpy dtype used in the subhalo structured array.
   """
   f=open('code-C/subs.h','w')
   f.write('/* Contains struct definition for properties of subhalos. */\n\n#include <stdbool.h>\n\nstruct struct_sub {\n')
   for key in D_sub.fields.keys():
      var_type=str(D_sub[key])
      # First deal with the awkward arrays
      if '(3,)' in var_type:
         f.write('    double '+key+'[3];\n')
      # then process the simple types
      elif 'bool' in var_type:
         f.write('    bool '+key+';\n')         
      elif 'int' in var_type:
         f.write('    int '+key+';\n')
      elif 'float' in var_type:
         f.write('    double '+key+';\n')
      else:
         f.write('    '+var_type+' '+key+';\n')
   f.write('}; \n')
   f.close()
   return None
   
