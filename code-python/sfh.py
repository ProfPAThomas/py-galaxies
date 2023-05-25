"""
Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy
"""

import numpy as np
import h5py
from codetiming import Timer
from profiling import conditional_decorator

import commons

class C_sfh:
   """
   Class to generate and store the generic star formation history tables.
   
   Attributes
   ----------
   """
   def __init__(self,parameters):
      """
      Generates the reference structure for storing the star formation histories in
      logarithmic bins (for each snapshot/time step combination). In the code galaxy
      structures are adjusted with respect to this structure at each step.
      t[MAXSNAPS][STEPS][SFH_NBIN] : float
         Time to present (i.e. z=0 ?) at the low-z edge of the bin (code units)
      dt[MAXSNAPS][STEPS][SFH_NBIN] : float
         Time width of the bin (code units)
      n_bin[MAXSNAPS][STEPS][SFH_NBIN] : int
         Number of bins merged in each bin (only useful for the merging algorithm)
      i_bin[MAXSNAPS][STEPS] : int
         Last active bin
      """
      n_merge=parameters.SFH_n_merge
      self.n_merge=n_merge  # Save here for convenience
      assert n_merge > 2
      n_max=n_merge-1  # Needs a better name - the maximum number of bins retained per level

      # First let's determine the number of bins that we need.
      n_snap=parameters.n_snap
      n_dt_halo_arr=parameters.n_dt_halo  # Number of halo steps in each snapshot
      n_dt_gal_arr=parameters.n_dt_gal    # Number of galaxy steps per halo step, in each snapshot
      # For convenience, it is useful to record the starting timestep number for each snapshot
      i_dt_snap=np.full(n_snap,parameters.NO_DATA_INT)
      i_dt_snap[0]=0
      for i_snap in range(1,n_snap):
         i_dt_snap[i_snap]=i_dt_snap[i_snap-1]+n_dt_halo_arr[i_snap-1]*n_dt_gal_arr[i_snap-1]
      n_dt=i_dt_snap[-1]+n_dt_halo_arr[-1]*n_dt_gal_arr[-1]
      assert n_dt==np.dot(n_dt_halo_arr,n_dt_gal_arr)  # Total number of steps taken in the code
      self.n_dt=n_dt
      self.i_dt_snap=i_dt_snap
        
      # Now there can be at most SFH_n_merge_max bins of any given size, + temporarily 1 extra bin prior to merging
      # So this is the maximum number of SFH bins that we could need at each step
      # n_level is the number of levels of the merging hierarchy
      n_level=int(1+0.999999*np.log2(n_dt/n_max))
      # And this the maximum number of SFH bins that we could need
      n_bin_max=n_max*n_level+1

      # Assign arrays to hold SFH information; we will later truncate to actual size needed
      self.t=np.full([n_dt,n_bin_max],np.NAN)      # Time at low-redshift end of bin
      self.dt=np.full([n_dt,n_bin_max],np.NAN)     # Width of the bin
      # The following three are used to track the bin structure;
      # Strictly speaking, only the first (or second) is needed; the others are tautologous.
      self.level=np.full([n_dt,n_bin_max],parameters.NO_DATA_INT)   # Level in the merging hierarchy
      self.n_bin_in_level=np.full([n_dt,n_level],parameters.NO_DATA_INT)  # Number of bins at this level
      self.i_bin=np.full(n_dt,parameters.NO_DATA_INT) # Highest occupied bin

      # Read in times from parameter file
      t_snap_arr=parameters.t_snap
      dt_snap_arr=parameters.dt_snap
      dt_gal_arr=parameters.dt_gal
        
      # Bins will be extend from the earliest times to the present day.
      # As bins get merged, they get shunted from high to low indices, and get fatter in size.
      # We use the following array to hold how many bins exist at each level of the merging hierarchy
      n_bin_in_level=np.full(n_level,0)
      # And these arrays hold the properties of the bins at each stage
      t=np.full(n_bin_max,np.NAN)
      dt=np.full(n_bin_max,0.)
      level=np.full(n_bin_max,parameters.NO_DATA_INT)
      # Now run through the timesteps, creating bins, merging as necessary, and saving
      i_bin=0           # Number of bins used at this timestep
      i_bin_max=0       # Maximum number of bins ever used
      i_dt=0
      for i_snap in range(n_snap):
         assert self.i_dt_snap[i_snap]==i_dt
         t_snap=t_snap_arr[i_snap]
         dt_snap=dt_snap_arr[i_snap]
         dt_gal=dt_gal_arr[i_snap]
         n_dt_halo=n_dt_halo_arr[i_snap]          
         n_dt_gal=n_dt_gal_arr[i_snap]          
         for i_dt_gal in range(n_dt_gal*n_dt_halo):
            # Create a new bin at merge level 0
            # Note, because of python counting, do not update until end of loop
            assert i_bin<n_bin_max  # Should not exceed our estimated maximum
            level[i_bin]=0
            n_bin_in_level[0]+=1
            t[i_bin]=t_snap+(i_dt_gal+1-n_dt_gal*n_dt_halo)*dt_gal # +1 because records time at lower edge of time bin.
            dt[i_bin]=dt_gal
            # Now run through levels, merging as required
            # Index of final bin at this level (counting backwards from high indices)
            j_bin=i_bin+1
            for i_level in range(n_level):
               j_bin-=n_bin_in_level[i_level]
               if n_bin_in_level[i_level] == n_merge:
                  # Need to merge bins j_bin and j_bin+1
                  # This essentially here means first resetting the times and counts of bins at each level,
                  # then shuffling bins downward.
                  i_bin-=1
                  n_bin_in_level[i_level+1]+=1
                  n_bin_in_level[i_level]-=2
                  level[j_bin]+=1
                  level[j_bin+1:-1]=level[j_bin+2:]
                  level[-1]=parameters.NO_DATA_INT
                  t[j_bin]=t[j_bin+1]
                  t[j_bin+1:-1]=t[j_bin+2:]
                  t[-1]=np.NAN 
                  dt[j_bin]+=dt[j_bin+1]
                  dt[j_bin+1:-1]=dt[j_bin+2:]
                  dt[-1]=0.
                  j_bin+=1      # First bin at this level has been pushed upwards by one slot
            # Merging complete
            # Update bin count from largest index used to number of bins used
            i_bin+=1
            i_bin_max=max(i_bin,i_bin_max)
            # Save results
            self.t[i_dt,:]=t
            self.dt[i_dt,:]=dt
            self.n_bin_in_level[i_dt,:]=n_bin_in_level
            self.i_bin[i_dt]=i_bin
            self.level[i_dt,:]=level
            # For testing
            if parameters.verbosity >=2: self.__repr__(n_step=[i_dt,i_dt+1])
            i_dt +=1
      assert i_dt==n_dt       
      # Can now truncate arrays to actual size used
      self.t=self.t[:,:i_bin_max]
      self.dt=self.dt[:,:i_bin_max]
      self.level=self.level[:,:i_bin_max]
      self.n_bin=i_bin_max   # Need to add 1 for the internal arrays that may temporarily have an extra bin.

      # Save the results
      # Would like to save time units but astropy units seem incompatible with hdf5.
      # Create mask corresponding to snapshots; we want the end, not the beginning of the snaps
      mask=i_dt_snap.copy()
      mask[:-1]=mask[1:]-1  # This is the index of the last timestep in this snapshot
      mask[-1]=-1
      sfh_file = h5py.File(parameters.sfh_file,'w')
      sfh_file.attrs['Number_of_time_bins']=self.n_bin
      dset=sfh_file.create_dataset('t',data=self.t[mask]*parameters.time_internal_to_output,compression='gzip')
      dset.attrs['Description']='Time at lower redshift edge of bin'
      dset=sfh_file.create_dataset('dt',data=self.dt[mask]*parameters.time_internal_to_output,compression='gzip')
      dset.attrs['Description']='Width of time bin'
      dset=sfh_file.create_dataset('i_bin',data=self.i_bin[mask],compression='gzip')
      dset.attrs['Description']='Number of active time bins'
      sfh_file.close()
      
      return None

   def __repr__(self,n_step=[0,7]):
      """
      Prints out the structure of the first n_step entries in the SFH arrays.
      Useful mainly as a check that things have been implemented properly.
      Would be nice to produce graphics at some point.
      """
      for i_step in range(n_step[0],n_step[1]):
         i_bin=self.i_bin[i_step]
         print('Step ',i_step,':')
         print('Level',end=': ')
         for i in range(i_bin):
            print(self.level[i_step,i],end=', ')
         print('\ndt',end=': ')
         for i in range(i_bin):
            print('{:.3f}'.format(self.dt[i_step,i]),end=', ')
         print('\nt',end=': ')
         for i in range(i_bin):
            print('{:.3f}'.format(self.t[i_step,i]),end=', ')
         print('\n')
      return ''

def F_sfh_update_bins(gals,sfh,parameters):
   """
   Merges galaxy bins if required.
   Parameters
   ----------
   gals : obj : np.array(n_gal)
      Structured array of galaxies
   sfh : obj : class C_SFH
      Instance of C_SFH class containing information about the time-binning structure
   """
   i_dt=commons.load('i_dt')   # This should hold the ministep ID BEFORE updating
   #print('F_sfh_update_bins: i_dt (before updating) =',i_dt)
   n_merge=sfh.n_merge
   n_bin_in_level=sfh.n_bin_in_level[i_dt].copy()
   level=np.full([sfh.n_bin+1],parameters.NO_DATA_INT)
   level[:-1]=sfh.level[i_dt,:]
   i_bin=sfh.i_bin[i_dt].copy()
   # Create new bin
   level[i_bin]=0
   n_bin_in_level[0]+=1
   # Run through levels merging as required
   j_bin=i_bin+1
   for i_level in range(len(n_bin_in_level)):
      j_bin-=n_bin_in_level[i_level]
      if n_bin_in_level[i_level] == n_merge:
         # Need to merge bins j_bin and j_bin+1
         # This essential here means first resetting the times and counts of bins at each level,
         # then shuffling bins downward.
         i_bin-=1
         n_bin_in_level[i_level+1]+=1
         n_bin_in_level[i_level]-=2
         level[j_bin]+=1
         level[j_bin+1:-1]=level[j_bin+2:]
         level[-1]=parameters.NO_DATA_INT
         # Now combine the data.
         # Would this be faster (and make the code look simpler) if all the SFH data was a sub-array?
         gals['mass_stars_bulge_sfh'][:,j_bin]+=gals['mass_stars_bulge_sfh'][:,j_bin+1]
         gals['mass_stars_bulge_sfh'][:,j_bin+1:-1]=gals['mass_stars_bulge_sfh'][:,j_bin+2:]
         gals['mass_stars_bulge_sfh'][:,-1]=0.
         gals['mass_metals_stars_bulge_sfh'][:,j_bin]+=gals['mass_metals_stars_bulge_sfh'][:,j_bin+1]
         gals['mass_metals_stars_bulge_sfh'][:,j_bin+1:-1]=gals['mass_metals_stars_bulge_sfh'][:,j_bin+2:]
         gals['mass_metals_stars_bulge_sfh'][:,-1]=0.
         gals['mass_stars_disc_sfh'][:,j_bin]+=gals['mass_stars_disc_sfh'][:,j_bin+1]
         gals['mass_stars_disc_sfh'][:,j_bin+1:-1]=gals['mass_stars_disc_sfh'][:,j_bin+2:]
         gals['mass_stars_disc_sfh'][:,-1]=0.
         gals['mass_metals_stars_disc_sfh'][:,j_bin]+=gals['mass_metals_stars_disc_sfh'][:,j_bin+1]
         gals['mass_metals_stars_disc_sfh'][:,j_bin+1:-1]=gals['mass_metals_stars_disc_sfh'][:,j_bin+2:]
         gals['mass_metals_stars_disc_sfh'][:,-1]=0.
         j_bin+=1      # First bin at this level has been pushed upwards by one slot
      
