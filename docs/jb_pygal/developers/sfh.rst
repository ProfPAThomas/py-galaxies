Star formation histories
========================

This section describes the use of an evolving time bin structure to capture star formation histories (SFHs).  This was first implemented by `Shamshiri etal (2015) <https://arxiv.org/abs/1501.05649>`_ although previously used by `Yates etal (2013) <https://arxiv.org/abs/1305.7231>`_.  The use or otherwise of SFH bins is controlled by a run-time boolean parameter :code:`b_SFH`.

Galaxy formation in py-gal (and L-Galaxies) proceeds through a series of mini time-steps: in L-Galaxies these were of length :code:`dt_snap/STEPS`; in py-gal they are of length :code:`dt_gal` chosen to be less that an input parameter :code:`model_parameters:timestep_galaxies` such that an integral number fits within each halo timestep :code:`dt_halo` (which is likewise chosen to be less that an input parameter :code:`model_parameters:timestep_halo` such that an integral number fits within the snapshot interval :code:`dt_snap`).

Now there would typically be far too many mini-steps to permit recording of all of them as that would swamp memory, so bins are merged in the code to produce a hierarchy of levels, :math:`l=0,\ 1,\ 2,\ \ldots` where within each level :math:`2^l` bins have been merged together.  The maximum number of bins allowed at each level is controlled by the input parameter :code:`SFH_n_merge` (stored internally as :code:`sfh.n_merge`) which triggers merging of the two oldest bins at each level once this number is reached (the=us the maximum number of bins at each level is restricted to :code:`SFH_n_merge-1`.  Experience has shown that :code:`SFH_n_merge=3`, the minimum possible, is sufficient to be able to reconstruct SEDs of galaxies with high accuracy, except in the far UV where it may be necessary to keep more high-resolution bins.

Tracking time in the code
-------------------------

Being constant throughout the execution of the code, all these parameters are stored in the :code:`parameters` instance of the :code:`C_parameters` class.  There are the following quantities associated with snapshots, with all quantities in internal code units; they are set be the function :code:`F_set_dt` within :code:`misc.py`:

  * :code:`n_snap` -- the number of snapshots
  * :code:`t_snap[n_snap]` -- the age at the time of the snapshot
  * :code:`dt_snap[n_snap]` -- the difference in age between this snapshot and the previous one
  * :code:`n_dt_halo[n_snap]` -- the number of halo timesteps between this snapshot and the previous one
  * :code:`dt_halo[n_snap]` -- the corresponding halo timestep size
  * :code:`n_dt_gal[n_snap]` -- the number of galaxy timesteps during each halo timestep
  * :code:`dt_gal[n_snap]` -- the corresponding galaxy timestep size

All code relating to SFH binning is located in :code:`sfh.py`.  That contains a class :code:`C_SFH` which upon initialisation calculates and saves the SFH binning structure at each timestep:

  * :code:`n_level` -- the number of levels in the merging hierarchy
  * :code:`n_dt` -- the total number of galaxy timesteps
  * :code:`i_dt_snap[n_snap]` -- the index of the first galaxy timestep in this snapshot interval (i.e. the total number of galaxy timesteps at the time of the previous snapshot.
  * :code:`n_bin` -- the number of time bins stored in the output file (note that the working arrays need one extra element to account for bins created before merging)
  * :code:`i_bin[n_dt]` -- number of active SFH bins at each timestep
  * :code:`n_bin_in_level[n_dt,n_level]` -- number of bins in each level of the merging hierarchy at each timestep
  * :code:`level[n_dt,n_bin]` -- the level (see above) of the merging hierarchy of each time bin at each timestep
  * :code:`t[n_dt,n_bin]` -- the time at the lower redshift edge of each time bin at each timestep
  * :code:`dt[n_dt]` -- the width of the time bin at each timestep

As the galaxy results are output only on snapshots [#]_ then that is all we need to store in the SFH output file:

  * :code:`i_bin[n_snap]` -- number of active SFH bins at each snapshot
  * :code:`t[n_snap,n_bin]` -- the time at the lower redshift edge of each time bin at each snapshot
  * :code:`dt[n_dt]` -- the width of the time bin at each snapshot
    
.. [#] There is no reason why this needs to be so, but it is simplest and means that galaxies can be matched to halo properties. 

Creating the SFH binning arrays
-------------------------------

The following code block shows the construction of the SFH binning arrays mentioned above.  It is executed during the initialisation of :code:`C_sfh`.  The basic idea is that one loops through binning levels, merging cells if the number of bins in that level equals :code:`n_merge`; the rest is just messy book-keeping.  Note that later time bins are added to the end of the SFH bin arrays, so we need to count backwards when doing the merging.

.. code-block:: python3

      i_bin=0           # Number of bins used at this timestep
      i_bin_max=0       # Maximum number of bins ever used
      i_dt=0
      for i_snap in range(n_snap):
         self.i_dt_snap[i_snap]=i_dt
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

Merging SFH bins during galaxy evolution
----------------------------------------

Galaxies inherit the SFH binning structure from ancestors in the previous timestep.  Any new stars formed are added to a new SFH bin.  Then, at the end the timestep, these bins are merged, if required.  The merging algorithm mirrors that used initially to create the SFH binning arrays during initialisation.

.. code-block:: python3

   i_dt=commons.load('i_dt')   # This should hold the ministep ID BEFORE updating
   n_merge=sfh.n_merge
   n_bin_in_level=sfh.n_bin_in_level[i_dt]
   j_bin=sfh.i_bin[i_dt]+1
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
         gals['mass_stars_bulge_sfh'][j_bin]+=gals['mass_stars_bulge_sfh'][j_bin+1]
         gals['mass_stars_bulge_sfh'][j_bin+1:-1]=gals['mass_stars_bulge_sfh'][j_bin+2:]
         gals['mass_stars_bulge_sfh'][-1]=0.
         gals['mass_metals_stars_bulge_sfh'][j_bin]+=gals['mass_metals_stars_bulge_sfh'][j_bin+1]
         gals['mass_metals_stars_bulge_sfh'][j_bin+1:-1]=gals['mass_metals_stars_bulge_sfh'][j_bin+2:]
         gals['mass_metals_stars_bulge_sfh'][-1]=0.
         gals['mass_stars_disc_sfh'][j_bin]+=gals['mass_stars_disc_sfh'][j_bin+1]
         gals['mass_stars_disc_sfh'][j_bin+1:-1]=gals['mass_stars_disc_sfh'][j_bin+2:]
         gals['mass_stars_disc_sfh'][-1]=0.
         gals['mass_metals_stars_disc_sfh'][j_bin]+=gals['mass_metals_stars_disc_sfh'][j_bin+1]
         gals['mass_metals_stars_disc_sfh'][j_bin+1:-1]=gals['mass_metals_stars_disc_sfh'][j_bin+2:]
         gals['mass_metals_stars_disc_sfh'][-1]=0.
         j_bin+=1      # First bin at this level has been pushed upwards by one slot
