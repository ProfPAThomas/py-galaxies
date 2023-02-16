""" 
Program to read in the cooling tables in raw binary format and write them out in npz and HDF5 formats.
Also saves the temperatures and metallicities corresponding to the rows and columns.
NOte that there are a lot of numbers in here specific to these particular tables.
"""

import h5py
import numpy as np
input_dir='../../CoolFunctions/'
input_files=[ 
  'stripped_mzero.cie',
  'stripped_m-30.cie',
  'stripped_m-20.cie',
  'stripped_m-15.cie',
  'stripped_m-10.cie',
  'stripped_m-05.cie',
  'stripped_m-00.cie',
  'stripped_m+05.cie'
]
input_files_dtype=np.dtype([
   ('sd_logT',np.float32),('sd_ne', np.float32),('sd_nh',np.float32),('sd_nt',np.float32),
   ('sd_logLnet',np.float32),('sd_logLnorm',np.float32),('sd_logU',np.float32),('sd_logTau',np.float32),
   ('sd_logP12',np.float32),('sd_logRho24',np.float32),('sd_ci',np.float32),('sd_mubar',np.float32) ])
n_Z=8
n_T=91
# Metallicies in input tables with repect to solar.  Convert to absolute by adding Zsun=0.02
log10_Z=np.array([-np.Inf,-3.0,-2.0,-1.5,-1.0,-0.5,0.0,0.5])+np.log10(0.02)
assert len(log10_Z)==n_Z
output_file_prefix='Lambda_table'

log10_Lambda=np.empty([n_Z,n_T],np.float32)

# Read in cooling table
for i_Z in range(n_Z):
   data=np.loadtxt(input_dir+input_files[i_Z],dtype=input_files_dtype)
   log10_T=data['sd_logT'][:]
   log10_Lambda[i_Z,:] = data['sd_logLnorm'][:]

# Now simply save the results
# Numpy compressed binary format
np.savez(output_file_prefix+'.npz',log10_T=log10_T,log10_Z=log10_Z,log10_Lambda=log10_Lambda)
# HDF5 format
output_file = h5py.File(output_file_prefix+'.hdf5','w')
dataset0 = output_file.create_dataset('log10_Z',data=log10_Z,compression='gzip')
dataset1 = output_file.create_dataset('log10_T',data=log10_T,compression='gzip')
dataset2 = output_file.create_dataset('log10_Lambda',data=log10_Lambda,compression='gzip')
output_file.close()
