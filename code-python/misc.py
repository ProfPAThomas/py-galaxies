"""
Miscellaneous helper routines that do not find a place elsewhere.
"""

import astropy.units as u
import ctypes
import numpy as np

#----------------------------------------------------------------------------------------------------

def F_misc_set_dt(parameters):
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

def F_misc_create_all_headers_header_file():
   """
   Writes out code/all_headers.h which contains a list of all header files to be included in the C code.
   Excludes those specialised headers that are only used in a single routine: cooling.h and sfh.h

   Attributes
   ----------
   parameters : obj : C_parameters
       Instance of class containing global parameters
   """
   f=open('code-C/all_headers.h','w')
   f.write('/* Contains a list of header files to be included in the C routines. */\n\n')
   f.write('#include <math.h>\n')
   f.write('#include <stdbool.h>\n')
   f.write('#include <stdio.h>\n')
   f.write('#include <stdlib.h>\n')
   f.write('#include <string.h>\n')
   f.write('#include "gals.h"\n')
   f.write('#include "halos.h"\n')
   f.write('#include "parameters.h"\n')
   f.write('#include "subs.h"\n')
   f.write('#include "variables.h"\n')
   # proto.h has to come last in order not to generate warnings about multiple struct definitions
   f.write('#include "proto.h"\n')
   f.close()
   return None

#----------------------------------------------------------------------------------------------------

def F_misc_create_variables_structure_and_header_file(variables_dict):
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

      def __init__(self,variables_dict):
         for key, value in variables_dict.items():
            if type(value)==bool or type(value)==int or type(value)==float:
               exec('self.'+key+'='+str(value))

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

   return C_variables
   
#----------------------------------------------------------------------------------------------------

def F_misc_create_Makefile(parameters):
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
   f.write('\tgcc -Wall ')
   if parameters.b_debug:
      f.write('-fsanitize=undefined -g ')
   else:
      f.write('-O ')
   f.write('-shared -o lib_C.so \\\n')
   if parameters.b_SFH:
      f.write('\t\t-D SFH \\\n')
   f.write('\t\tbh_agn.c \\\n')
   f.write('\t\tcooling.c \\\n')
   #f.write('\t\tdynamics.c \\\n')
   f.write('\t\thalos.c \\\n')
   f.write('\t\tmergers.c \\\n')
   f.write('\t\tprocess_snap.c \\\n')
   if parameters.b_SFH:
      f.write('\t\tsfh.c \\\n')
   f.write('\t\tstar_formation_and_feedback.c\n')

   f.close()
   return None

