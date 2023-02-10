import h5py
import numpy as np

"""
Galaxies are stored in structured numpy arrays, so we will not create any instances of this class.   
But we use it to define the dtype of that strucured array.
"""
# Create the dtype that we will need to store galaxy properties.
D_gal=[
   ('graph_ID',np.int32),
   ('snap_ID',np.int32),
   ('halo_gid',np.int32),
   ('halo_sid',np.int32),
   ('sub_gid',np.int32),
   ('sub_sid',np.int32),
   ('gal_gid',np.int32),      # The unique identifier for this galaxy within this graph; should match location in output file
   ('desc_gid',np.int32),
   ('first_prog_gid',np.int32),
   ('next_prog_gid',np.int32),
   ('b_exists',np.bool),
   ('b_merger',np.bool),
   #('pos',np.float32,(3,)),  # Not sure where to set this in the main program ...
   #('vel',np.float32,(3,)),  # ... and can be found in halo/subhalo catalogues (apart maybe for orphans)
   ('mass_stars_bulge',np.float32),
   ('mass_metals_stars_bulge',np.float32),
   ('mass_stars_disc',np.float32),
   ('mass_metals_stars_disc',np.float32),
   ('mass_cold_gas',np.float32),
   ('mass_metals_cold_gas',np.float32)
]

def F_gal_template(parameters):
   """
   Creates a template for the galaxies
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
   template['b_exists']=True
   template['b_merger']=False
   #template['pos']=np.nan
   #template['vel']=np.nan
   template['mass_stars_bulge']=0.
   template['mass_metals_stars_bulge']=0.
   template['mass_stars_disc']=0.
   template['mass_metals_stars_disc']=0.
   template['mass_cold_gas']=0.
   template['mass_metals_cold_gas']=0.
   return template

class C_gal_output:
   """
   This class contains the attributes and methods for the galaxy output file.
   Attributes
   ----------

   Methods
   -------
   __init__
   append - add halos to output buffer
   close - flush io buffer then close HDF5 file
   flush - flush output buffer to HDF5 dataset
   
   """
   def __init__(self,parameters):
      """
      Opens the galaxy output file.
      Creates the galaxy output buffer.
      Creates the HDF5 galaxy dataset.

      Parameters:
      -----------
      parameters : obj : C_parameters
         Contains the gloabal run paramters.
      """
      # Open file for output
      self.gal_file = h5py.File(parameters.galaxy_output_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.D_param['performance']['n_HDF5_io_rec']['Value']
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
      dtype.append(('b_exists',np.bool))
      #dtype.append(('pos',np.float32,(3,)))
      #dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass_stars_bulge',np.float32))
      dtype.append(('mass_metals_stars_bulge',np.float32))
      dtype.append(('mass_stars_disc',np.float32))
      dtype.append(('mass_metals_stars_disc',np.float32))
      dtype.append(('mass_cold_gas',np.float32))
      dtype.append(('mass_metals_cold_gas',np.float32))
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
         halos - list of C_halo objects to be output
         parameters - C_parameters class file containing the global run parameters
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
         self.io_buffer[self.i_rec]['mass_stars_bulge'] = gals[i_gal]['mass_stars_bulge']
         self.io_buffer[self.i_rec]['mass_metals_stars_bulge'] = gals[i_gal]['mass_metals_stars_bulge']
         self.io_buffer[self.i_rec]['mass_stars_disc'] = gals[i_gal]['mass_stars_disc']
         self.io_buffer[self.i_rec]['mass_metals_stars_disc'] = gals[i_gal]['mass_metals_stars_disc']
         self.io_buffer[self.i_rec]['mass_cold_gas']= gals[i_gal]['mass_cold_gas']
         self.io_buffer[self.i_rec]['mass_metals_cold_gas']= gals[i_gal]['mass_metals_cold_gas']
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
