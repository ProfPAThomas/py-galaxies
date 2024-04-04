"""
Subhalo numpy structured array dtype definition and template for initialisation, and python class for output.

All physical quantities are in internal code units.

Methods now moved to subs.c

Attributes
----------
graph_ID : int
    The graph_ID within the input file.
snap_ID : int
     The snapshot ID currently being processed.
halo_gid : int
    The host halo location within the graph.
halo_sid : int
    The host halo location within the snapshot.
sub_gid : int
    The subhalo location within the graph.
sub_sid : int
    The subhalo location within the snapshot.
b_done : bool
    Whether or not the subhalo has been fully processed.
desc_end_gid : int
    The index relative to the graph at which this subhalo's descendents end.
desc_halo_sid : int
    The index relative to the next snapshot of the descendant of the host halo of this subhalo.
desc_start_gid : int
    The index relative to the graph at which this subhalo's descendents start.
gal_central_sid : int
    The location in the current galaxy array of the most massive galaxy in the subhalo.
gal_end_sid : int
    The location in the current galaxy array of the last galaxy in the subhalo (+1 because of python indexing)
gal_next_sid : int
    Galaxy counter used when updating galaxies.
gal_start_sid : int
    The location in the current galaxy array of the first galaxy in the subhalo
half_mass_radius : float
    The radius containing half the total mass in the DM-only sim
half_mass_virial_speed : float
    The circular speed at the half-mass radius
mass : float
    The DM-only mass of the subhalo.
mass_baryon : float
    Mass of baryons within the subhalo, inclusive of galaxies
mass_gas_hot : float
    The mass of hot gas
mass_metals_gas_hot : float
    The mass of metals in hot gas
mass_metals_stars : float
    The mass of metals in stars
mass_stars : float
    The mass of stars
n_desc : int
    The number of direct descendants.
n_dt : int
    Number of times that this halo has been processed this snapshot
n_gal : int
    The number of galaxies in the subhalo.
pos : float[3]
    The position of the subhalo.
rms_speed : float
    The rms speed of the subhalo in the DM-only sim.
tau_dyn : float
    The dynamical time at twice the half-mass radius.
temperature : float
    The virial temperature as derived from the virial speed.
vel : float[3]
    The velocity of the subhalo.

"""

import ctypes
import h5py
import numpy as np

# First define the dtype for subhalo properties in a way that is compatible with ctypes
D_sub=np.dtype([
   ('graph_ID',ctypes.c_int),
   ('snap_ID',ctypes.c_int),
   ('halo_gid',ctypes.c_int),
   ('halo_sid',ctypes.c_int),
   ('sub_gid',ctypes.c_int),
   ('sub_sid',ctypes.c_int),
   ('n_desc',ctypes.c_int),
   ('desc_start_gid',ctypes.c_int),
   ('desc_end_gid',ctypes.c_int),
   ('desc_halo_sid',ctypes.c_int),  # The descendent of the halo this subhalo sits in
   ('desc_main_sid',ctypes.c_int),  # The main descendent of this subhalo
   ('n_dt',ctypes.c_int),
   ('b_done',ctypes.c_bool),
   ('mass',ctypes.c_double),
   ('pos',ctypes.c_double*3),
   ('vel',ctypes.c_double*3),
   ('rms_speed',ctypes.c_double),
   ('half_mass_radius',ctypes.c_double),
   ('half_mass_virial_speed',ctypes.c_double),
   # Derived properties
   ('temperature',ctypes.c_double),
   ('tau_dyn',ctypes.c_double),
   # SAM properties
   ('mass_baryon',ctypes.c_double),
   ('mass_gas_hot',ctypes.c_double),
   ('mass_metals_gas_hot',ctypes.c_double),
   ('mass_stars',ctypes.c_double),
   ('mass_metals_stars',ctypes.c_double),
   ('n_gal',ctypes.c_int),
   ('gal_start_sid',ctypes.c_int),
   ('gal_end_sid',ctypes.c_int),
   ('gal_next_sid',ctypes.c_int),
   ('gal_central_sid',ctypes.c_int)],
   align=True)
    
#----------------------------------------------------------------------------------------------------

def F_subs_create_header_file():
   """
   Creates a C struct definition that matches the subhalo dtype.
   Writes out to code/subs.h

   Attributes
   ----------
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
   
#--------------------------------------------------------------------------------------------------

def F_subs_initialise(subs_this_snap,graph_ID,snap_ID,graph,parameters):
   """
   Read in the subhalo properties from the graph, including ranges for descendants;
    
   Parameters
   ----------
   subs_this_snap : obj : D_sub[n_sub]
       The properties of subhalos in this graph/snapshot
   graph_ID : int
       The graph_ID within the input data file.
   snap_ID : int
       The snapshot ID currently being processed.
   sub_gid : int
       The subhalo ID relative to the graph.
   graph : obj : C_graph
       The graph containing this halo.
   parameters : obj : C_parameters
       The global parameters for this SAM run.
   """
   halo_offset=graph.snap_first_halo_gid[snap_ID]
   sub_offset=graph.snap_first_sub_gid[snap_ID]
   n_sub=len(subs_this_snap)
   # Read in halo properties from graph instance.  These should already be in internal code units.
   subs_this_snap['graph_ID'] = graph_ID
   subs_this_snap['snap_ID'] = snap_ID
   for i_sub in range(n_sub):
      
      # Indexing
      sub_sid = i_sub
      sub_gid = sub_offset+i_sub
      subs_this_snap[i_sub]['halo_gid'] =  graph.sub_host_gid[sub_gid]
      subs_this_snap[i_sub]['halo_sid'] = subs_this_snap[i_sub]['halo_gid'] - halo_offset
      subs_this_snap[i_sub]['sub_gid'] =  sub_gid
      # Do we want a whole separate entry just for sub_sid, which is the same as i_sub?
      subs_this_snap[i_sub]['sub_sid'] = sub_sid

      # Graph properties
      subs_this_snap[i_sub]['n_desc'] = graph.sub_n_desc[sub_gid]
      subs_this_snap[i_sub]['desc_start_gid'] = graph.sub_first_desc_gid[sub_gid]
      subs_this_snap[i_sub]['desc_end_gid'] = subs_this_snap[i_sub]['desc_start_gid'] + subs_this_snap[i_sub]['n_desc']

      # Astrophysical properties
      subs_this_snap[i_sub]['mass'] = graph.sub_mass[sub_gid]
      subs_this_snap[i_sub]['pos'] = graph.sub_pos[sub_gid]
      subs_this_snap[i_sub]['vel'] = graph.sub_vel[sub_gid]
      subs_this_snap[i_sub]['rms_speed'] = graph.sub_rms_speed[sub_gid]
      subs_this_snap[i_sub]['half_mass_radius'] = graph.sub_half_mass_radius[sub_gid]
      subs_this_snap[i_sub]['half_mass_virial_speed'] = (0.5*parameters.c_G*subs_this_snap[i_sub]['mass']/subs_this_snap[i_sub]['half_mass_radius'])**(0.5)
      subs_this_snap[i_sub]['temperature'] = subs_this_snap[i_sub]['half_mass_virial_speed']**2 * parameters.c_half_mass_virial_speed_to_temperature
      subs_this_snap[i_sub]['tau_dyn'] = 2.*subs_this_snap[i_sub]['half_mass_radius']/subs_this_snap[i_sub]['half_mass_virial_speed']

#--------------------------------------------------------------------------------------------------

def F_subs_mass_baryon(sub,gals):
       """
       Calculates the total baryonic mass of the subhalo, including galaxies.
       Returns value rather than setting it because being used as a check on simpler method.

       Parameters
       ----------
       gals : obj : D_gal[]
          Array of records for the galaxies in this subhalo.

       Returns
       -------
       float
          The baryonic mass of the subhalo, inclusive of galaxies.
       """
       mass_baryon = sub['mass_gas_hot'] + sub['mass_stars'] + \
                       np.sum(gals[sub['gal_start_sid']:sub['gal_end_sid']]['mass_baryon'])
       return mass_baryon
        
#--------------------------------------------------------------------------------------------------

def F_subs_template(parameters):
   """
   Define template for new subhalos

   Parameters
   ----------
   parameters : obj : C_parameters
      Contains the global run parameters.

   Returns
   -------
   template : obj : D_sub
      A one row numpy structured array with dtype D_sub 
   """
   # Define template for new subhalo instances
   template=np.empty(1,dtype=np.dtype(D_sub,align=True))
   template['graph_ID'] = parameters.NO_DATA_INT
   template['snap_ID'] = parameters.NO_DATA_INT
   template['halo_gid'] = parameters.NO_DATA_INT
   template['halo_sid'] = parameters.NO_DATA_INT
   template['sub_gid'] = parameters.NO_DATA_INT
   template['sub_sid'] = parameters.NO_DATA_INT
   template['n_desc'] = 0
   template['desc_start_gid'] = parameters.NO_DATA_INT
   template['desc_end_gid'] = parameters.NO_DATA_INT
   template['desc_halo_sid'] = parameters.NO_DATA_INT
   template['desc_main_sid'] = parameters.NO_DATA_INT
   template['n_dt'] = 0        # Number of times that this halo has been processed
   template['b_done'] = False  # Has this halo been fully processed or not.
   template['mass'] = 0.
   template['pos'] = [0.,0.,0.]
   template['vel'] = [0.,0.,0.]
   template['rms_speed'] = 0.
   template['half_mass_radius'] = 0.
   template['half_mass_virial_speed'] = 0.
   template['temperature'] = 0.
   template['tau_dyn'] = 0.
   template['mass_baryon']=0.
   template['mass_gas_hot']=0.
   template['mass_metals_gas_hot']=0.
   template['mass_stars']=0.
   template['mass_metals_stars']=0.
   template['n_gal'] = 0
   template['gal_start_sid'] = parameters.NO_DATA_INT
   template['gal_end_sid'] = parameters.NO_DATA_INT
   template['gal_next_sid'] = parameters.NO_DATA_INT
   template['gal_central_sid'] = parameters.NO_DATA_INT
   return template

#--------------------------------------------------------------------------------------------------

class C_sub_output:
   
   """
   This class contains the attributes and methods for the subhalo output files.

   Attributes
   ----------
   sub_file : obj : File
      HDF5 file for subhalo output
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
      Opens the subhalo output file.
      Creates the subhalo output buffer.
      Creates the HDF5 halo dataset.

      Parameters
      ----------
      parameters : obj : C_parameters
         Contains the global run parameters.
      """
      # Open file for output
      self.sub_file = h5py.File(parameters.subhalo_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.n_HDF5_io_rec
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_gid',np.int32))
      dtype.append(('sub_gid',np.int32))
      dtype.append(('pos',np.float32,(3,)))
      dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass',np.float32))
      dtype.append(('rms_speed',np.float32))
      dtype.append(('half_mass_virial_speed',np.float32))
      dtype.append(('temperature',np.float32))
      dtype.append(('mass_gas_hot',np.float32))
      dtype.append(('mass_metals_gas_hot',np.float32))
      dtype.append(('mass_stars',np.float32))
      dtype.append(('mass_metals_stars',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.sub_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and closes the subhalo output file.
      """
      self.flush()
      # self.dataset.close() # There does not seem to be a need to close datasets.
      self.sub_file.close()
      return None

   def flush(self):
      """
      Writes io buffer to the HDF5 dataset and resets.
      """
      self.dataset.resize((self.dataset.shape[0]+self.i_rec,))
      self.dataset[-self.i_rec:]=self.io_buffer[:self.i_rec]
      self.i_rec=0
      return None

   def append(self,subs,parameters):
      """
      Extracts the quantities desired for halo_output and adds them to the io buffer, flushing if required.

      Parameters
      ----------
         subs : obj : C_sub[]
            list of C_sub objects to be output.
         parameters : obj : C_parameters
            The global run parameters.
      """
      for i_sub in range(len(subs)):
         self.io_buffer[self.i_rec]['graph_ID'] = subs[i_sub]['graph_ID']
         self.io_buffer[self.i_rec]['snap_ID'] = subs[i_sub]['snap_ID']
         self.io_buffer[self.i_rec]['halo_gid'] = subs[i_sub]['halo_gid']
         self.io_buffer[self.i_rec]['sub_gid'] = subs[i_sub]['sub_gid']
         self.io_buffer[self.i_rec]['pos'] = subs[i_sub]['pos'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec]['vel'] = subs[i_sub]['vel'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['mass'] = subs[i_sub]['mass'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['rms_speed'] = subs[i_sub]['rms_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['half_mass_virial_speed'] = subs[i_sub]['half_mass_virial_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec]['temperature'] = subs[i_sub]['temperature'] * parameters.temperature_internal_to_output
         self.io_buffer[self.i_rec]['mass_gas_hot']= subs[i_sub]['mass_gas_hot'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_gas_hot']= subs[i_sub]['mass_metals_gas_hot'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_stars']= subs[i_sub]['mass_stars'] * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec]['mass_metals_stars']= subs[i_sub]['mass_metals_stars'] * parameters.mass_internal_to_output
         self.i_rec+=1
         if self.i_rec == self.n_rec: self.flush()
      return None
