"""
Halo numpy structured array dtype definition and template for initialisation, and python class for output.

All physical quantities are in internal code units.

Methods now moved to halos.c

Attributes
----------
graph_ID : int
    The graph_ID (from HDF5 group).
snap_ID : int
    The snapshot ID currently being processed.
halo_gid : int
    The halo location within the graph of the halo currently being processed.
halo_sid : int
    The halo location within the snapshot of the halo currently being processed.
b_done : bool
    Whether or not the halo has been fully processed.
desc_end_gid : int
    The index at which this halo's descendents end (+1 as is usual for python).
desc_main_sid : int
    The main descendant location in this snapshot (ie in halos_this_snap).
desc_start_gid : int
    The index at which this halo's descendents start.
half_mass_radius : float
    The radius containing half the total mass of the halo in the DM-only sim.
half_mass_virial_speed : float
    The circular speed at the half-mass radius.
mass : float
    The DM-only mass of this halo.
mass_baryon : float
    Mass of baryons within the halo, inclusive of subhalos and galaxies.
mass_baryon_from_progenitors : float
    Total mass of all the baryons contained within the progenitor halos.
mass_baryon_delta_dthalo : float
    The accretion needed per halo timestep to bring the baryon content up to the universal mean.
mass_from_progenitors : float 
    Total DM-only mass of all the progenitor halos.
mass_gas_eject : float
    The mass of ejected gas
mass_gas_hot : float
    The mass of hot gas in the halo, exclusive of subhalos.
mass_stars : float
    The mass of stars in the halo, exclusive of subhalos.
mass_metals_gas_eject : float
    The mass of metals in ejected gas.
mass_metals_gas_hot : float
    The mass of metals in hot gas in the halo, exclusive of subhalos.
mass_metals_stars : float
    The mass of metals in stars in the halo, exclusive of subhalos.
n_desc : int
    The number of direct descendants of this halo
n_dt : int
    Number of times that this halo has been processed this snapshot
n_gal : int
    The number of galaxies in the halo, inclusive of subhalos
n_orphan : int
    The number of orphan galaxies (i.e. galaxies not contained in subhalos)
n_sub : int
    The number of subhalos.
orphan_start_sid : int
    The location of the first orphan galaxy within the current snapshot galaxy array
orphan_next_sid : int
    To track the processing of orphan galaxies during the halo_update phase.
pos : float[3]
    The position of the halo.
rms_radius : float
    The rms radius of the halo particles in the DM-only sim.
rms_speed : float
    The rms speed of the halo particles in the DM-only sim.
sub_central_sid : int
    The location in the snapshot of the subhalo at the centre of the halo (if any)
sub_end_gid : int
    The location in the graph of the last subhalo (+1 because of python indexing)
sub_end_sid : int
    The location in the snapshot of the last subhalo (+1 because of python indexing)
sub_mass : float[n_sub]
    The DM-only masses of the subhalos
sub_rel_pos : float[n_sub,3]
    The positions of the subhalos relative to that of the halo
sub_rel_vel : float[n_sub,3]
    The velocities of the subhalos relative to that of the halo
sub_start_gid : int
    The location in the graph of the first subhalo
sub_start_sid : int
    The location in the snapshot of the first subhalo
tau_dyn : float
    Twice the dynamical time at the half-mass radius (= dynamical time at twice the half mass radius for isothermal sphere).
temperature : float
    The temperature as derived from the virial speed.
vel : float[3]
    Velocity of the halo
"""

import ctypes
import h5py
import numpy as np

# First define the dtype for halo properties in a way that is compatible with ctypes
D_halo=np.dtype([
   ('graph_ID',ctypes.c_int),
   ('snap_ID',ctypes.c_int),
   ('halo_gid',ctypes.c_int),
   ('halo_sid',ctypes.c_int),
   ('n_desc',ctypes.c_int),
   ('desc_start_gid',ctypes.c_int),
   ('desc_end_gid',ctypes.c_int),
   ('desc_main_sid',ctypes.c_int),
   ('n_dt',ctypes.c_int),
   ('b_done',ctypes.c_bool),
   ('mass',ctypes.c_double),
   ('mass_from_progenitors',ctypes.c_double),
   ('pos',ctypes.c_double*3),
   ('vel',ctypes.c_double*3),
   ('rms_speed',ctypes.c_double),
   ('half_mass_radius',ctypes.c_double),
   ('half_mass_virial_speed',ctypes.c_double),
   ('rms_radius',ctypes.c_double),
   # Derived properties
   ('temperature',ctypes.c_double),
   ('tau_dyn',ctypes.c_double),
   # SAM properties
   ('mass_baryon',ctypes.c_double),
   ('mass_baryon_from_progenitors',ctypes.c_double),
   ('mass_baryon_delta_dthalo',ctypes.c_double),
   ('mass_gas_eject',ctypes.c_double),
   ('mass_metals_gas_eject',ctypes.c_double),
   ('mass_gas_hot',ctypes.c_double),
   ('mass_metals_gas_hot',ctypes.c_double),
   ('mass_stars',ctypes.c_double),
   ('mass_metals_stars',ctypes.c_double),
   ('n_sub',ctypes.c_int),
   ('sub_start_sid',ctypes.c_int),
   ('sub_end_sid',ctypes.c_int),
   ('sub_central_sid',ctypes.c_int),
   ('n_gal',ctypes.c_int),
   ('gal_start_sid',ctypes.c_int),
   ('gal_end_sid',ctypes.c_int),
   ('n_orphan',ctypes.c_int),
   ('orphan_start_sid',ctypes.c_int),
   ('orphan_end_sid',ctypes.c_int),
   ('orphan_next_sid',ctypes.c_int)],
   align=True)

#----------------------------------------------------------------------------------------------------

def F_halos_create_header_file():
   """
   Creates a C struct definition that matches the subhalo dtype.
   Writes out to code/subs.h

   Attributes
   ----------
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
   
#--------------------------------------------------------------------------------------------------

def F_halos_initialise(halos_this_snap,graph_ID,snap_ID,graph,parameters):
   """
   Read in the halos properties from the graph, including ranges for decendants, subhalos and galaxies.
 
   Parameters
   ----------
   halos_this_snap : obj : D_halo[n_halo]
       The properties of halos in this graph/snapshot
   graph_ID : str
       The graph_ID (from HDF5 group).
   snap_ID : int
       The snapshot ID currently being processed.
   graph : an instance of the class C_graph
       The graph containing this halo.
   parameters : an instance of the class C_parameters
       The global parameters for this SAM run.
   """
   halo_offset=graph.snap_first_halo_gid[snap_ID]
   n_halo=len(halos_this_snap)
   # Read in halo properties from graph instance.  These should already be in internal code units.
   halos_this_snap['graph_ID'] = graph_ID
   halos_this_snap['snap_ID'] = snap_ID
   for i_halo in range(n_halo):
      
      # Indexing
      halo_sid = i_halo
      halo_gid = halo_offset+i_halo
      halos_this_snap[i_halo]['halo_gid'] = halo_gid
      # Do we want a whole separate entry just for halo_sid, which is the same as i_halo?
      halos_this_snap[i_halo]['halo_sid'] = halo_sid

      # Graph properties
      halos_this_snap[i_halo]['n_desc'] = graph.halo_n_desc[halo_gid]
      halos_this_snap[i_halo]['desc_start_gid'] = graph.halo_first_desc_gid[halo_gid]
      halos_this_snap[i_halo]['desc_end_gid'] = halos_this_snap[i_halo]['desc_start_gid'] + halos_this_snap[i_halo]['n_desc']
   
      # Astrophysical properties
      halos_this_snap[i_halo]['mass'] = graph.halo_mass[halo_gid]
      halos_this_snap[i_halo]['pos'] = graph.halo_mean_pos[halo_gid]
      halos_this_snap[i_halo]['vel'] = graph.halo_mean_vel[halo_gid]
      halos_this_snap[i_halo]['half_mass_radius'] = graph.halo_half_mass_radius[halo_gid]
      halos_this_snap[i_halo]['rms_radius'] = graph.halo_rms_radius[halo_gid]
      halos_this_snap[i_halo]['rms_speed'] = graph.halo_rms_speed[halo_gid]
      # Derived quantities
      # Using v^2=GM/r but for half mass
      halos_this_snap[i_halo]['half_mass_virial_speed'] = (0.5*parameters.c_G*halos_this_snap[i_halo]['mass']/halos_this_snap[i_halo]['half_mass_radius'])**(0.5)
      halos_this_snap[i_halo]['temperature'] = halos_this_snap[i_halo]['half_mass_virial_speed']**2 * parameters.c_half_mass_virial_speed_to_temperature
      halos_this_snap[i_halo]['tau_dyn'] = 2.*halos_this_snap[i_halo]['half_mass_radius']/halos_this_snap[i_halo]['half_mass_virial_speed']
   
      # Subhalos
      n_sub = graph.halo_n_sub[halo_gid]
      halos_this_snap[i_halo]['n_sub'] = n_sub
      halos_this_snap[i_halo]['sub_start_sid'] = graph.halo_first_sub_gid[halo_gid]-graph.snap_first_sub_gid[snap_ID]
      halos_this_snap[i_halo]['sub_end_sid'] = halos_this_snap[i_halo]['sub_start_sid']+n_sub
      ################################################################################
      # Have removed 'sub_mass', 'sub_rel_pos' and 'sub_rel_vel' as variable length! #
      # Will need to look them up in subhalos_this_snap as required.                 #
      ################################################################################

   return None

#--------------------------------------------------------------------------------------------------

def F_halos_template(parameters):
   """
   Define template for new halos

   Parameters
   ----------
   parameters : obj : C_parameters
      Contains the global run parameters.

   Returns
   -------
   template : obj : D_halo
      A one row numpy structured array with dtype D_halo
   """
   template=np.empty(1,dtype=D_halo)
   template['graph_ID'] = parameters.NO_DATA_INT
   template['snap_ID'] = parameters.NO_DATA_INT
   template['halo_gid'] = parameters.NO_DATA_INT
   template['halo_sid'] = parameters.NO_DATA_INT
   template['n_desc'] = 0
   template['desc_start_gid'] = parameters.NO_DATA_INT
   template['desc_end_gid'] = parameters.NO_DATA_INT
   template['desc_main_sid'] = parameters.NO_DATA_INT
   template['n_dt'] = 0        # Number of times that this halo has been processed
   template['b_done'] = False  # Has this halo been fully processed or not.
   template['mass'] = 0.
   template['mass_from_progenitors'] = 0.
   template['pos'] = [0.,0.,0.]
   template['vel'] = [0.,0.,0.]
   template['rms_speed'] = 0.
   template['half_mass_radius'] = 0.
   template['half_mass_virial_speed'] = 0.
   template['rms_radius'] = 0.
   template['temperature'] = 0.
   template['tau_dyn'] = 0.
   template['mass_baryon']=0.
   template['mass_baryon_from_progenitors']=0.
   template['mass_baryon_delta_dthalo']=0.
   template['mass_gas_eject']=0.
   template['mass_metals_gas_eject']=0.
   template['mass_gas_hot']=0.
   template['mass_metals_gas_hot']=0.
   template['mass_stars']=0.
   template['mass_metals_stars']=0.
   template['n_sub'] = 0
   template['sub_start_sid'] = parameters.NO_DATA_INT
   template['sub_end_sid'] = parameters.NO_DATA_INT
   template['sub_central_sid'] = parameters.NO_DATA_INT
   template['n_gal'] = 0     # Total number of galaxies in halo + subhalos
   template['gal_start_sid'] = parameters.NO_DATA_INT
   template['gal_end_sid'] = parameters.NO_DATA_INT
   template['n_orphan'] = 0  # Galaxies not associated with a subhalo
   template['orphan_start_sid'] = parameters.NO_DATA_INT
   template['orphan_end_sid'] = parameters.NO_DATA_INT
   template['orphan_next_sid'] = parameters.NO_DATA_INT
   return template
   
#--------------------------------------------------------------------------------------------------

class C_halo_output:
   
   """
   This class contains the attributes and methods for the halo output files. 

   Attributes
   ----------
   halo_file : obj : File
      HDF5 file for halo output
   i_rec : int
      Counter for how many records have been created.
   io_buffer : obj : D_gal[n_rec]
      Storage for halo records prior to outputting. 
   n_rec : int
      Number records to be buffered before outputting.
   dataset : obj : HDF5 dataset
      HDF5 dataset to which the data is output.
   """
   def __init__(self,parameters):
      """
      Opens the halo output file.
      Creates the halo output buffer.
      Creates the HDF5 halo dataset.

      Parameters
      ----------
      parameters : obj : C_parameters
         Contains the global run parameters.
      Returns
      -------
      None
      """
      # Open file for output
      self.halo_file = h5py.File(parameters.halo_file,'w')
      # Counter for and max number of records in io buffer
      self.i_rec = 0
      self.n_rec = parameters.n_HDF5_io_rec
      # dtype of io buffer
      dtype=[]
      dtype.append(('graph_ID',np.int32))
      dtype.append(('snap_ID',np.int32))
      dtype.append(('halo_gid',np.int32))
      dtype.append(('pos',np.float32,(3,)))
      dtype.append(('vel',np.float32,(3,)))
      dtype.append(('mass',np.float32))
      dtype.append(('temperature',np.float32))
      dtype.append(('rms_speed',np.float32))
      dtype.append(('half_mass_virial_speed',np.float32))
      dtype.append(('mass_baryon',np.float32))
      dtype.append(('mass_gas_hot',np.float32))
      dtype.append(('mass_metals_gas_hot',np.float32))
      dtype.append(('mass_gas_eject',np.float32))
      dtype.append(('mass_metals_gas_eject',np.float32))
      dtype.append(('mass_stars',np.float32))
      dtype.append(('mass_metals_stars',np.float32))
      # Create halo io buffer
      print('self.n_rec =',self.n_rec)
      self.io_buffer=np.empty(self.n_rec,dtype=dtype)
      # Create HDF5 dataset
      self.dataset = self.halo_file.create_dataset('Halos', \
         (0,),maxshape=(None,),dtype=dtype,compression='gzip')

   def close(self):
      """
      Empties the halo io buffer, closes the halo dataset, and
      closes the halo output file.
      """
      self.flush()
      # self.dataset.close() # There does not seem to be a need to close datasets.
      self.halo_file.close()
      return None

   def flush(self):
      """
      Writes io buffer to the HDF5 dataset and resets.
      """
      self.dataset.resize((self.dataset.shape[0]+self.i_rec,))
      self.dataset[-self.i_rec:]=self.io_buffer[:self.i_rec]
      self.i_rec=0
      return None

   def append(self,halos,parameters):
      """
      Extracts the quantities desired for halo_output and adds them to the io buffer,
      flushing if required.

      Parameters
      ----------
         halos : obj : C_halo[]
            list of C_halo objects to be output.
         parameters : obj : C_parameters
            The global run parameters.
      """
      n_halo=len(halos)
      n_offset=0
      while n_halo>0:
         n_add = min(n_halo, self.n_rec-self.i_rec)
         self.io_buffer[self.i_rec:self.i_rec+n_add]['graph_ID'] = halos[n_offset:n_offset+n_add]['graph_ID']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['snap_ID'] = halos[n_offset:n_offset+n_add]['snap_ID']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['halo_gid'] = halos[n_offset:n_offset+n_add]['halo_gid']
         self.io_buffer[self.i_rec:self.i_rec+n_add]['pos'] = halos[n_offset:n_offset+n_add]['pos'] * parameters.length_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['vel'] = halos[n_offset:n_offset+n_add]['vel'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass'] = halos[n_offset:n_offset+n_add]['mass'] * parameters.mass_internal_to_output      
         self.io_buffer[self.i_rec:self.i_rec+n_add]['temperature'] = halos[n_offset:n_offset+n_add]['temperature'] * parameters.temperature_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['rms_speed'] = halos[n_offset:n_offset+n_add]['rms_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['half_mass_virial_speed'] = halos[n_offset:n_offset+n_add]['half_mass_virial_speed'] * parameters.speed_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_baryon']= halos[n_offset:n_offset+n_add]['mass_baryon']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_gas_hot'] = halos[n_offset:n_offset+n_add]['mass_gas_hot']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_gas_hot'] = halos[n_offset:n_offset+n_add]['mass_metals_gas_hot']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_gas_eject'] = halos[n_offset:n_offset+n_add]['mass_gas_eject']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_gas_eject'] = halos[n_offset:n_offset+n_add]['mass_metals_gas_eject']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_stars'] = halos[n_offset:n_offset+n_add]['mass_stars']  * parameters.mass_internal_to_output
         self.io_buffer[self.i_rec:self.i_rec+n_add]['mass_metals_stars'] = halos[n_offset:n_offset+n_add]['mass_metals_stars']  * parameters.mass_internal_to_output
         n_halo -= n_add
         n_offset += n_add
         self.i_rec += n_add
         if self.i_rec == self.n_rec: self.flush()
      return None
