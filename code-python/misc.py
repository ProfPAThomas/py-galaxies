"""
Miscellaneous helper routines that do not find a place elsewhere.
"""

import astropy.units as u
import ctypes
import numpy as np

#----------------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------------

def F_create_variables_structure_and_header_file(variables_dict):
   """
   Creates a variables structure containing all the variables that we want to pass to C;
   Populates with current values;
   Creates the equivalent variables.h file for use with C.

   Attributes
   ----------
   variables_dict : obj : python dictionary
      Dictionary of variables and values to be used in creating the struct.
   """

   # First create a list of tuples describing the fields of the structure
   fields=[]
   for key, value in variables_dict.items():
      if type(value)==bool:
         fields.append((key,ctypes.c_bool))
      elif type(value)==int:
         fields.append((key,ctypes.c_int))
      elif type(value)==float:
         fields.append((key,ctypes.c_double))
      # Not sure how to handle strings in C so test this before uncommenting
      # elif type(value)==str:
      #    fields.append((key,ctypes.c_wchar*len(value)))
      else:
         raise ValueError('Unsupported type')
   
   # Next create a python Structure to hold those variables
   class C_variables(ctypes.Structure):
      _fields_=fields
   variables=C_variables()
   # And populate with the existing values; not sure if there is a better way.
   for key, value in variables_dict.items():
      if type(value)==bool or type(value)==int or type(value)==float:
         exec('variables.'+key+'='+str(value))

   # Write out a C header file with the structure definition
   f=open('code-C/variables.h','w')
   f.write('/* Runtime variables. */\n\n')
   f.write('#include <stdbool.h>\n\n')
   f.write('struct struct_var {\n')
   for key, value in variables_dict.items():
      if type(value)==bool:
         f.write('    bool '+key+';\n')
      elif type(value)==int:
         f.write('    int '+key+';\n')
      elif type(value)==float:
         f.write('    double '+key+';\n')
      # Not sure how to handle strings in C so test this before uncommenting
      # elif type(value)==str:
      #    f.write('    str '+key+'['+str(len(value))+'];\n')
   f.write('};\n')
   f.close()

   return variables
   
#----------------------------------------------------------------------------------------------------

def F_create_cooling_header_file(cooling_table):
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
   parameters : obj : C_cooling
       Instance of class containing cooling tables
   """
   f=open('code-C/cooling.h','w')
   f.write('/* Cooling tables (fixed throughout run). */\n\n')

   log10_T_table=cooling_table.log10_T_table
   n_T=len(log10_T_table)
   f.write('#define n_T '+str(n_T)+'\n')
   f.write('const double log10_T_table['+str(n_T)+'] = {')
   for i_T in range(n_T):
      f.write(str(log10_T_table[i_T])+', ')
   f.write('};\n\n')

   log10_Z_table=cooling_table.log10_Z_table
   n_Z=len(log10_Z_table)
   f.write('#define n_Z '+str(n_Z)+'\n')
   f.write('const double log10_Z_table['+str(n_Z)+'] = {')
   # Negative infinity will cause issues, so set first value to a very small number
   f.write('-100., ')
   for i_Z in range(1,n_Z):
      f.write(str(log10_Z_table[i_Z])+', ')
   f.write('};\n\n')

   log10_Lambda_table=cooling_table.log10_Lambda_table
   f.write('const double log10_Lambda_table['+str(n_Z)+']['+str(n_T)+'] = {\n')
   for i_Z in range(n_Z):
      f.write('{')
      for i_T in range(n_T):
         f.write(str(log10_Lambda_table[i_Z,i_T])+', ')
      f.write('},\n')
   f.write('};\n\n')   

   f.close()
   return None

#----------------------------------------------------------------------------------------------------

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
   
#----------------------------------------------------------------------------------------------------

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
   
#----------------------------------------------------------------------------------------------------

def F_create_parameters_header_file(parameters):
   """
   Writes out all the attributes of parameters to code/parameters.h.

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
   # Would like to use const here rather than static, but that seems to break the linker
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

#----------------------------------------------------------------------------------------------------

def F_create_sfh_header_file(sfh,parameters):
   """
   Creates a header file for the SFH containing parameters and fixed arrays describing the time bins.

   Attributes
   ----------
   sfh : obj : C_sfh instance
       Instance of class containing SFH bin parameters and arrays describing the time bins.
   parameters : obj : C_parameters
       Instance of class containing global parameters
   """
   f=open('code-C/sfh.h','w')
   f.write('/* SFH tables (fixed throughout run). */\n\n')

   n_dt=sfh.n_dt
   n_bin=sfh.n_bin
   n_level=sfh.n_level
   n_merge=sfh.n_merge
   f.write('#define n_dt '+str(n_dt)+'\n')
   f.write('#define n_bin '+str(n_bin)+'\n')
   f.write('#define n_level '+str(n_level)+'\n')
   f.write('#define n_merge '+str(n_merge)+'\n\n')

   n_snap=parameters.n_snap
   f.write('// Index of first galaxy timestep in each snapshot.\n')
   f.write('const int i_dt_snap['+str(n_snap)+'] = {')
   for i_snap in range(n_snap):
      f.write(str(sfh.i_dt_snap[i_snap])+', ')
   f.write('};\n\n')

   f.write('// Number of SFH bins used at each galaxy timestep.\n')   
   f.write('const int i_bin_all['+str(n_dt)+'] = {')
   for i_dt in range(n_dt):
      f.write(str(sfh.i_bin[i_dt])+', ')
   f.write('};\n\n')

   f.write('// Number of SFH bins in each level of the merging hierarchy.\n')   
   f.write('const int n_bin_in_level_all['+str(n_dt)+']['+str(n_level)+'] = {\n')
   for i_dt in range(n_dt):
      f.write('{')
      for i_level in range(n_level):
         f.write(str(sfh.n_bin_in_level[i_dt,i_level])+', ')
      f.write('},\n')
   f.write('};\n\n')   

   f.write('// Level in merging hierarchy of each SFH bin.\n')   
   f.write('const int level_all['+str(n_dt)+']['+str(n_bin)+'] = {\n')
   for i_dt in range(n_dt):
      f.write('{')
      for i_bin in range(n_bin):
         f.write(str(sfh.level[i_dt,i_bin])+', ')
      f.write('},\n')
   f.write('};\n\n')   

   # Not currently needed in C routines
   # f.write('// Age of the Universe at the low-z edge of the SFH bin (code units).\n')   
   # f.write('const double t['+str(n_dt)+']['+str(n_bin)+'] = {\n')
   # for i_dt in range(n_dt):
   #    f.write('{')
   #    for i_bin in range(n_bin):
   #       f.write(str(sfh.t[i_dt,i_bin])+', ')
   #    f.write('},\n')
   # f.write('};\n\n')   

   # Not currently needed in C routines
   # f.write('// Time width of the SFH bin (code units).\n')   
   # f.write('const double dt['+str(n_dt)+']['+str(n_bin)+'] = {\n')
   # for i_dt in range(n_dt):
   #    f.write('{')
   #    for i_bin in range(n_bin):
   #       f.write(str(sfh.dt[i_dt,i_bin])+', ')
   #    f.write('},\n')
   # f.write('};\n\n')   

   f.close()
   return None

#----------------------------------------------------------------------------------------------------

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
   
#----------------------------------------------------------------------------------------------------

def F_create_Makefile(parameters):
   """
   Writes out the Makefile to be used to compile the C library

   Attributes
   ----------
   parameters : obj : C_cooling
       Instance of class containing cooling tables
   """
   f=open('code-C/Makefile','w')
   f.write('# This Makefile is dynamically generated at runtime - do not edit!\n')
   f.write('# Any options should be specified in the input.yml file and written out here in F_create_Makefile within misc.py.\n\n')
   
   f.write('all: build_library\n\n')

   f.write('build_library:\n')
   f.write('\tgcc -Wall -O -shared -o lib_c.so \\\n')
   if parameters.b_SFH:
      f.write('\t\t-D SFH \\\n')
   f.write('\t\tbh_agn.c \\\n')
   f.write('\t\tcooling.c \\\n')
   f.write('\t\tmergers.c \\\n')
   if parameters.b_SFH:
      f.write('\t\tsfh.c \\\n')
   f.write('\t\tstar_formation_and_feedback.c\n')

   f.close()
   return None

