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

def F_create_parameters_header_file(parameters):
   """
   Writes out all the attributes of parameters to code/parameters.h.

   Attributes
   ----------
   parameters : obj : C_parameters
       Instance of class containing global parameters
   """
   attributes=[a for a in dir(parameters) if not a.startswith('__') and not callable(getattr(parameters, a))]
   f=open('parameters.h','w')
   f.write('/* Runtime parameters (fixed throughout run). */\n\
\n\
#include <stdbool.h>\n\
\n\
struct {\n\
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
      # Now the bits that we want to extract
      elif a_type == 'bool':
         if value==True:
            f.write('    '+a_type+' '+a+'=true;\n')
         else:
            f.write('    '+a_type+' '+a+'=false;\n')
      elif 'dict' in a_type:
         continue
      elif 'int' in a_type:
         f.write('    int '+a+'='+str(value)+';\n')
      elif 'float' in a_type:
         f.write('    float '+a+'='+str(value)+';\n')
      elif a_type == 'str':
         f.write('    char* '+a+'="'+str(value)+'";\n')
      else:
         f.write('    '+a_type+' '+a+'='+str(value)+';\n')
   f.write('}; parameters;\n')
   f.close()
   return None
