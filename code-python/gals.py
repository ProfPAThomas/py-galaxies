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

import h5py
import numpy as np

import commons
b_SFH=commons.load('b_SFH')
if b_SFH: sfh_n_bin=commons.load('sfh.n_bin')

# Create the dtype that we will need to store galaxy properties.
D_gal=[
   ('graph_ID',np.int32),
   ('snap_ID',np.int32),
   ('halo_gid',np.int32),
   ('halo_sid',np.int32),
   ('sub_gid',np.int32),
   ('sub_sid',np.int32),
   ('gal_gid',np.int32),      # The unique identifier for this galaxy within this graph; should match location in output file
   ('gal_sid',np.int32),
   ('desc_gid',np.int32),
   ('first_prog_gid',np.int32),
   ('next_prog_gid',np.int32),
   ('b_exists',bool),
   ('v_vir',np.float32),      # Virial speed of host halo
   ('mass_stars_bulge',np.float32),
   ('mass_metals_stars_bulge',np.float32),
   ('mass_stars_disc',np.float32),
   ('mass_metals_stars_disc',np.float32),
   ('mass_gas_cold',np.float32),
   ('mass_metals_gas_cold',np.float32),
   ('mass_BH',np.float32),
   ('mass_metals_BH',np.float32),   # Metals have no meaning in a BH but useful for tracking
   ('mass_baryon',np.float32),      # Includes BHs.  Effectively equivalent to total mass of galaxy (assuming no DM).
   ('radius_gas_cold',np.float32),  # Exponential disc radius
   ('radius_stars_disc',np.float32), # Exponential disc radius
   ('radius_stars_bulge',np.float32), # Half mass radius
   ('SFR_dt',np.float32),
   ('SFR_snap',np.float32)]
# Can only append one item at a time :-(
if b_SFH:
   D_gal.append(('mass_stars_bulge_sfh',np.float32,sfh_n_bin+1))
   D_gal.append(('mass_metals_stars_bulge_sfh',np.float32,sfh_n_bin+1))
   D_gal.append(('mass_stars_disc_sfh',np.float32,sfh_n_bin+1))
   D_gal.append(('mass_metals_stars_disc_sfh',np.float32,sfh_n_bin+1))

def F_gal_template(parameters):
   """
   Creates a template for the galaxies

   Parameters
   ----------
   parameters : obj : C_parameters
      Contains the global run parameters.

   Returns
   -------
   : obj : D_gal
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
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_gid',np.int32))
      dtype.append(('sub_gid',np.int32))
      dtype.append(('gal_gid',np.int32))
      dtype.append(('desc_gid',np.int32))
      dtype.append(('first_prog_gid',np.int32))
      dtype.append(('next_prog_gid',np.int32))
      dtype.append(('b_exists',bool))
      dtype.append(('mass_stars_bulge',np.float32))
      dtype.append(('mass_metals_stars_bulge',np.float32))
      dtype.append(('mass_stars_disc',np.float32))
      dtype.append(('mass_metals_stars_disc',np.float32))
      dtype.append(('mass_gas_cold',np.float32))
      dtype.append(('mass_metals_gas_cold',np.float32))
      dtype.append(('mass_BH',np.float32))
      dtype.append(('radius_gas_cold',np.float32))
      dtype.append(('radius_stars_disc',np.float32))
      dtype.append(('radius_stars_bulge',np.float32))
      dtype.append(('SFR_dt',np.float32))
      dtype.append(('SFR_snap',np.float32))
      if b_SFH:
         # Note: don't include final array entry here, as that is a working value.
         dtype.append(('mass_stars_bulge_sfh',np.float32,sfh_n_bin))
         dtype.append(('mass_metals_stars_bulge_sfh',np.float32,sfh_n_bin))
         dtype.append(('mass_stars_disc_sfh',np.float32,sfh_n_bin))
         dtype.append(('mass_metals_stars_disc_sfh',np.float32,sfh_n_bin))
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
      for i_gal in range(len(gals)):
         self.io_buffer[self.i_rec]['graph_ID'] = gals[i_gal]['graph_ID']
         self.io_buffer[self.i_rec]['snap_ID'] = gals[i_gal]['snap_ID']
         self.io_buffer[self.i_rec]['halo_gid'] = gals[i_gal]['halo_gid']
         self.io_buffer[self.i_rec]['sub_gid'] = gals[i_gal]['sub_gid']
         self.io_buffer[self.i_rec]['gal_gid'] = gals[i_gal]['gal_gid']
         self.io_buffer[self.i_rec]['desc_gid'] = gals[i_gal]['desc_gid']
         self.io_buffer[self.i_rec]['first_prog_gid'] = gals[i_gal]['first_prog_gid']
         self.io_buffer[self.i_rec]['next_prog_gid'] = gals[i_gal]['next_prog_gid']
         # Galaxy may have merged but still need to output it to avoid messing up indexing
         self.io_buffer[self.i_rec]['b_exists'] = gals[i_gal]['b_exists']  
         #self.io_buffer[self.i_rec]['pos'] = gals[i_gal]['pos']
         #self.io_buffer[self.i_rec]['vel'] = gals[i_gal]['vel']
         self.io_buffer[self.i_rec]['mass_stars_bulge'] = gals[i_gal]['mass_stars_bulge'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars_bulge'] = gals[i_gal]['mass_metals_stars_bulge'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_stars_disc'] = gals[i_gal]['mass_stars_disc'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars_disc'] = gals[i_gal]['mass_metals_stars_disc'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_gas_cold']= gals[i_gal]['mass_gas_cold'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_cold']= gals[i_gal]['mass_metals_gas_cold'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_BH']= gals[i_gal]['mass_BH'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['radius_gas_cold']= gals[i_gal]['radius_gas_cold'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['radius_stars_disc']= gals[i_gal]['radius_stars_disc'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['radius_stars_bulge']= gals[i_gal]['radius_stars_bulge'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['SFR_dt']=gals[i_gal]['SFR_dt'] * parameters.mass_internal_to_output / parameters.time_internal_to_output
         self.io_buffer[self.i_rec]['SFR_snap']=gals[i_gal]['SFR_snap'] * parameters.mass_internal_to_output / parameters.time_internal_to_output
         if b_SFH:
            self.io_buffer[self.i_rec]['mass_stars_bulge_sfh'] = gals[i_gal]['mass_stars_bulge_sfh'][:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec]['mass_metals_stars_bulge_sfh'] = gals[i_gal]['mass_metals_stars_bulge_sfh'][:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec]['mass_stars_disc_sfh'] = gals[i_gal]['mass_stars_disc_sfh'][:-1] * parameters.mass_internal_to_output
            self.io_buffer[self.i_rec]['mass_metals_stars_disc_sfh'] = gals[i_gal]['mass_metals_stars_disc_sfh'][:-1] * parameters.mass_internal_to_output
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
