"""
Galaxy numpy structured array d_type definition and template for initialisation, and python class for output.

All physical quantities are in internal code units.

Attributes
----------
graph_ID : int
   The graph_ID (from HDF5 group).
snap_ID : int
   The snapshot ID currently being processed.
halo_gid : int
   The halo ID relative to the graph.
halo_sid : int
   The halo ID relative to the snapshot.
sub_gid : int
   The subhalo ID relative to the graph
sub_sid : int
   The subhalo ID relative to the snapshot.
gal_gid : int
   The galaxy location relative to the start of this graph in the output file.
b_exists : bool
   Whether or not the galaxy exists (it may have merged with another galaxy).
desc_gid : int
   The index relative to the graph of this galaxy's descendant.
first_prog_gid: int
   The index relative to the graph of this galaxy's first progenitor.
mass_gas_cold : float
   The mass of cold gas, inclusive of metals.
mass_metals_gas_cold : float
   The mass of metals in the cold gas.
mass_metals_stars_bulge : float
   The mass of metals in bulge stars.
mass_metals_stars_disc : float
   The mass of metals in disc stars.
mass_stars_bulge : float
   The mass of bulge stars, inclusive of metals.
mass_stars_disc : float
   The mass of disc stars, inclusive of metals.
next_prog_gid : int
   The index relative to the graph of the next progenitor of the descendant galaxy.
radius_gas_cold : float
   The disc scale radius for the cold gas.
radius_stars_bulge : float
   The half mass radius for the bulge stars (Jaffe profile).
radius_stars_disc : float
   The disc scale radius for the stellar disc.
SFR_dt : float
   The star formation rate in the last galaxy timestep
SFR_snap : float
   The star formation rate averaged over the last snapshot
v_vir : float
   The half-mass circular speed of the host subhalo (or halo, if no subhalo).
"""

import ctypes
import h5py
import numpy as np

import commons
b_SFH=commons.load('b_SFH')
if b_SFH: sfh_n_bin=commons.load('sfh_n_bin')

# Create the dtype that we will need to store galaxy properties.
_quantities=[
   ('graph_ID',ctypes.c_int),
   ('snap_ID',ctypes.c_int),
   ('halo_gid',ctypes.c_int),
   ('halo_sid',ctypes.c_int),
   ('sub_gid',ctypes.c_int),
   ('sub_sid',ctypes.c_int),
   ('gal_gid',ctypes.c_int),      # The unique identifier for this galaxy within this graph; should match location in output file
   ('gal_sid',ctypes.c_int),
   ('desc_gid',ctypes.c_int),
   ('first_prog_gid',ctypes.c_int),
   ('next_prog_gid',ctypes.c_int),
   ('b_exists',ctypes.c_bool),
   ('v_vir',ctypes.c_double),      # Virial speed of host halo
   ('mass_stars_bulge',ctypes.c_double),
   ('mass_metals_stars_bulge',ctypes.c_double),
   ('mass_stars_disc',ctypes.c_double),
   ('mass_metals_stars_disc',ctypes.c_double),
   ('mass_gas_cold',ctypes.c_double),
   ('mass_metals_gas_cold',ctypes.c_double),
   ('mass_BH',ctypes.c_double),
   ('mass_metals_BH',ctypes.c_double),   # Metals have no meaning in a BH but useful for tracking
   ('mass_baryon',ctypes.c_double),      # Includes BHs.  Effectively equivalent to total mass of galaxy (assuming no DM).
   ('radius_gas_cold',ctypes.c_double),  # Exponential disc radius
   ('radius_stars_disc',ctypes.c_double), # Exponential disc radius
   ('radius_stars_bulge',ctypes.c_double), # Half mass radius
   ('SFR_dt',ctypes.c_double),
   ('SFR_dt_start',ctypes.c_double),
   ('SFR_snap',ctypes.c_double)]
# Can only append one item at a time :-(
if b_SFH:
   _quantities.append(('mass_stars_bulge_sfh',ctypes.c_double,sfh_n_bin+1))
   _quantities.append(('mass_metals_stars_bulge_sfh',ctypes.c_double,sfh_n_bin+1))
   _quantities.append(('mass_stars_disc_sfh',ctypes.c_double,sfh_n_bin+1))
   _quantities.append(('mass_metals_stars_disc_sfh',ctypes.c_double,sfh_n_bin+1))
D_gal=np.dtype(_quantities,align=True)
   
#----------------------------------------------------------------------------------------------------

def F_gals_create_header_file():
   """
   Creates a C struct definition that matches the galaxy dtype.
   Writes out to code/gals.h

   Attributes
   ----------
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
   
#------------------------------------------------------------------------------------------------------------

def F_gals_template(parameters):
   """
   Creates a template for the galaxies

   Parameters
   ----------
   parameters : obj : C_parameters
      Contains the global run parameters.

   Returns
   -------
   template : obj : D_gal
      A one row numpy structured array with dtype D_gal 
   """
   NDI=parameters.NO_DATA_INT
   template=np.empty(1,dtype=D_gal)
   template['graph_ID']=NDI
   template['snap_ID']=NDI
   template['halo_gid']=NDI
   template['halo_sid']=NDI
   template['sub_gid']=NDI
   template['sub_sid']=NDI
   template['gal_gid']=0 # Because there are no galaxies prior to the first snapshot; updated each snap.
   template['gal_sid']=NDI
   template['desc_gid']=NDI
   template['first_prog_gid']=NDI
   template['next_prog_gid']=NDI
   template['b_exists']=False
   template['v_vir']=0.
   template['mass_stars_bulge']=0.
   template['mass_metals_stars_bulge']=0.
   template['mass_stars_disc']=0.
   template['mass_metals_stars_disc']=0.
   template['mass_gas_cold']=0.
   template['mass_metals_gas_cold']=0.
   template['mass_BH']=0.
   template['mass_metals_BH']=0.
   template['mass_baryon']=0.
   template['radius_gas_cold']=0.
   template['radius_stars_disc']=0.
   template['radius_stars_bulge']=0.
   template['SFR_dt']=0.
   template['SFR_dt_start']=0.
   template['SFR_snap']=0.
   if b_SFH:
      template['mass_stars_bulge_sfh']=0.
      template['mass_metals_stars_bulge_sfh']=0.
      template['mass_stars_disc_sfh']=0.
      template['mass_metals_stars_disc_sfh']=0.
   return template

class C_gal_output:
   """
   This class contains the attributes and methods for the galaxy output file.

   Attributes
   ----------
   gal_file : obj : File
      HDF5 file for galaxy output
   i_rec : int
      Counter for how many records have been created.
   io_buffer : obj " D_gal[n_rec]
      Storage for subhalo records prior to outputting. 
   n_rec : int
      Number records to be buffered before outputting.
   dataset : obj : HDF5 dataset
      HDF5 dataset to which the data is output.
   """
   def __init__(self,parameters):
      """
      Opens the galaxy output file.
      Creates the galaxy output buffer.
      Creates the HDF5 galaxy dataset.

      Parameters
      ----------
      parameters : obj : C_parameters
         Contains the global run paramters.
      """
      # Open file for output
      self.gal_file = h5py.File(parameters.galaxy_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.n_HDF5_io_rec
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',ctypes.c_int))
      dtype.append(('snap_ID',ctypes.c_int))
      dtype.append(('halo_gid',ctypes.c_int))
      dtype.append(('sub_gid',ctypes.c_int))
      dtype.append(('gal_gid',ctypes.c_int))
      dtype.append(('desc_gid',ctypes.c_int))
      dtype.append(('first_prog_gid',ctypes.c_int))
      dtype.append(('next_prog_gid',ctypes.c_int))
      dtype.append(('b_exists',ctypes.c_bool))
      dtype.append(('mass_stars_bulge',ctypes.c_double))
      dtype.append(('mass_metals_stars_bulge',ctypes.c_double))
      dtype.append(('mass_stars_disc',ctypes.c_double))
      dtype.append(('mass_metals_stars_disc',ctypes.c_double))
      dtype.append(('mass_gas_cold',ctypes.c_double))
      dtype.append(('mass_metals_gas_cold',ctypes.c_double))
      dtype.append(('mass_BH',ctypes.c_double))
      dtype.append(('radius_gas_cold',ctypes.c_double))
      dtype.append(('radius_stars_disc',ctypes.c_double))
      dtype.append(('radius_stars_bulge',ctypes.c_double))
      dtype.append(('SFR_dt',ctypes.c_double))
      dtype.append(('SFR_dt_start',ctypes.c_double))
      dtype.append(('SFR_snap',ctypes.c_double))
      if b_SFH:
         # Note: don't include final array entry here, as that is a working value.
         dtype.append(('mass_stars_bulge_sfh',ctypes.c_double,sfh_n_bin))
         dtype.append(('mass_metals_stars_bulge_sfh',ctypes.c_double,sfh_n_bin))
         dtype.append(('mass_stars_disc_sfh',ctypes.c_double,sfh_n_bin))
         dtype.append(('mass_metals_stars_disc_sfh',ctypes.c_double,sfh_n_bin))
      # Create halo io buffer
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.gal_file.create_dataset('Galaxies', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and
      closes the halo output file
      """
      self.flush()
      # self.dataset.close() # There does not seem to be a need to close datasets.
      self.gal_file.close()
      return None

   def flush(self):
      """
      Writes io buffer to the HDF5 dataset and resets.
      """
      self.dataset.resize((self.dataset.shape[0]+self.i_rec,))
      self.dataset[-self.i_rec:]=self.io_buffer[:self.i_rec]
      self.i_rec=0
      return None

   def append(self,gals,parameters):
      """
      Extracts the quantities desired for halo_output and adds them to the io buffer,
      flushing if required.

      Parameters
      ----------
         gals : obj: D_gal[n_gal]
            Array of D_gal records to be output
         parameters : obj : C_parameters
            Contains the global run paramters.
      """
      n_gal=len(gals)
      n_offset=0
      while n_gal>0:
         n_add = min(n_gal, self.n_rec-self.i_rec)
         self.io_buffer[self.i_rec:self.i_rec+n_add]['graph_ID'] = gals[n_offset:n_offset+n_add]['graph_ID']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['snap_ID'] = gals[n_offset:n_offset+n_add]['snap_ID']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['halo_gid'] = gals[n_offset:n_offset+n_add]['halo_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['sub_gid'] = gals[n_offset:n_offset+n_add]['sub_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['gal_gid'] = gals[n_offset:n_offset+n_add]['gal_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['desc_gid'] = gals[n_offset:n_offset+n_add]['desc_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['first_prog_gid'] = gals[n_offset:n_offset+n_add]['first_prog_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['next_prog_gid'] = gals[n_offset:n_offset+n_add]['next_prog_gid']
         # Galaxy may have merged but still need to output it to avoid messing up indexing
         self.io_buffer[self.i_rec:self.i_rec+n_add]['b_exists'] = gals[n_offset:n_offset+n_add]['b_exists']  
         #self.io_buffer[self.i_rec:self.i_rec+n_add]['pos'] = gals[n_offset:n_offset+n_add]['pos']
         #self.io_buffer[self.i_rec:self.i_rec+n_add]['vel'] = gals[n_offset:n_offset+n_add]['vel']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_stars_bulge'] = gals[n_offset:n_offset+n_add]['mass_stars_bulge'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_stars_bulge'] = gals[n_offset:n_offset+n_add]['mass_metals_stars_bulge'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_stars_disc'] = gals[n_offset:n_offset+n_add]['mass_stars_disc'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_stars_disc'] = gals[n_offset:n_offset+n_add]['mass_metals_stars_disc'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_gas_cold']= gals[n_offset:n_offset+n_add]['mass_gas_cold'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_gas_cold']= gals[n_offset:n_offset+n_add]['mass_metals_gas_cold'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_BH']= gals[n_offset:n_offset+n_add]['mass_BH'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['radius_gas_cold']= gals[n_offset:n_offset+n_add]['radius_gas_cold'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['radius_stars_disc']= gals[n_offset:n_offset+n_add]['radius_stars_disc'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['radius_stars_bulge']= gals[n_offset:n_offset+n_add]['radius_stars_bulge'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['SFR_dt']=gals[n_offset:n_offset+n_add]['SFR_dt'] * parameters.mass_internal_to_output / parameters.time_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['SFR_dt_start']=gals[n_offset:n_offset+n_add]['SFR_dt_start'] * parameters.mass_internal_to_output / parameters.time_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['SFR_snap']=gals[n_offset:n_offset+n_add]['SFR_snap'] * parameters.mass_internal_to_output / parameters.time_internal_to_output
         if b_SFH:
            # if we include a range for the galaxies then the [:-1] syntax used up to v0.3 no longer seems to work.
            self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_stars_bulge_sfh'] = (gals[n_offset:n_offset+n_add]['mass_stars_bulge_sfh'])[:,:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_stars_bulge_sfh'] = (gals[n_offset:n_offset+n_add]['mass_metals_stars_bulge_sfh'])[:,:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_stars_disc_sfh'] = (gals[n_offset:n_offset+n_add]['mass_stars_disc_sfh'])[:,:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_stars_disc_sfh'] = (gals[n_offset:n_offset+n_add]['mass_metals_stars_disc_sfh'])[:,:-1] * parameters.mass_internal_to_output
         n_gal -= n_add
         n_offset += n_add
         self.i_rec += n_add
         if self.i_rec == self.n_rec: self.flush()
      return None
