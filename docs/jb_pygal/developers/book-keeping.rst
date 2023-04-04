Book-keeping of halos, subhalos and galaxies
============================================

This section outlines the book-keeping associated with:

* the ordering of galaxies
* the inheritance of properties from progenitors

I have omitted lines from the code blocks that are not relevant to this.


Looping over graphs and snapshots
---------------------------------


.. code-block:: python3
   
    # Loop over graphs
    for i_graph in range(min(parameters.n_graph,n_GRAPH)):
        graph = C_graph(i_graph,graph_file,parameters)
    
    # Loop over snapshots
    halos_last_snap = None
    subs_last_snap = None
    gals_last_snap = None
    for i_snap in graph.snap_ID:
           
        # Initialise halo and subhalo properties.
        halos_this_snap = [C_halo(i_graph,i_snap,i_halo,graph,parameters) for i_halo in 
                         graph.halo_start_gid[i_snap]+range(graph.n_halo_snap[i_snap])]
        subs_this_snap = None
        if graph.n_sub > 0:
            if graph.n_sub_snap[i_snap] > 0:
                subs_this_snap = [C_sub(i_graph,i_snap,i_sub,graph,parameters) 
                                     for i_sub in graph.sub_start_gid[i_snap]+range(graph.n_sub_snap[i_snap])]
        
        # Propagate information from progenitors to this generation and initialise galaxies
        gals_this_snap=F_update_halos(halos_last_snap,halos_this_snap,subs_last_snap,
                                          subs_this_snap,gals_last_snap,graph,parameters)
        del halos_last_snap
        del subs_last_snap
        del gals_last_snap

        # Process the halos
        F_process_halos(halos_this_snap,subs_this_snap,gals_this_snap,graph,parameters)
            
        # Once all halos have been done, output results
        # This could instead be done on a halo-by-halo basis in F_process_halos
        halo_output.append(halos_this_snap,parameters)
        if subs_this_snap != None: sub_output.append(subs_this_snap,parameters)
        if gals_this_snap != None: gal_output.append(gals_this_snap,parameters)
            
        # Rename this_snap data structures to last_snap
        halos_last_snap=halos_this_snap
        subs_last_snap=subs_this_snap
        gals_last_snap=gals_this_snap

        # Delete old references (so that create new objects on next snapshot)
        del halos_this_snap
        del subs_this_snap
        del gals_this_snap

    # Tidy up
    del halos_last_snap
    del subs_last_snap
    del gals_last_snap



Initialisation of halos and subhalos
------------------------------------

Within each snapshot, halos, subhalo and galaxies are stored as python lists.  When passing forward information from one snapshot to the next, two copies of these lists are required: one for the progenitor and one for the new snapshot; for subsequent processing, only the current snapshot is required.

.. code-block:: python3
   
  halos_this_snap = [C_halo(i_graph,i_snap,i_halo,graph,parameters) for i_halo in
                       graph.halo_start_gid[i_snap]+range(graph.n_halo_snap[i_snap])]

.. code-block:: python3
   
        subs_this_snap = None
        if graph.n_sub > 0:
            if graph.n_sub_snap[i_snap] > 0:
                subs_this_snap = [C_sub(i_graph,i_snap,i_sub,graph,parameters) 
                                     for i_sub in graph.sub_start_gid[i_snap]+range(graph.n_sub_snap[i_snap])]

During the initialisation, pointers the host halo and related subhalos, that are available from the graph, are read in.

Propagation of information from progenitors
-------------------------------------------

Note that this also generates the galaxy array:

.. code-block:: python3
   
        # Have to do this even if no progenitors in order to initialise galaxy array
        gals_this_snap=F_update_halos(halos_last_snap,halos_this_snap,subs_last_snap,
                                          subs_this_snap,gals_last_snap,graph,parameters)

Galaxies are stored in halo order.  Within each halo, we first have orphan galaxies (i.e. those who have lost their subhalos) followed by the galaxies associated with each subhalo, in subhalo order).

We need to do a first pass to push halo/subhalo properties and to determine the number of galaxies.  This also sets pointers in the halo and subhalo instances of where the associated galaxy and orphan galaxy counts start.

The first code block loops over halos, giving mass and hot gas to it's descendants in proportion to their overlap, and all galaxies to the main descendant (the one with the most overlap).  I have omitted some of the lines, for clarity; this just shows the overall structure of the block.
					  
.. code-block:: python3

    # Loop over halos
    if halos_last_snap != None:
       for halo in halos_last_snap:
          # First determine what fraction to give to each descendant
          desc_start_gid=halo.desc_start_gid
          desc_end_gid=halo.desc_end_gid
          if (halo.n_desc==0): 
             # For now just skip this halo; might want in future to log these occurrences
             continue
          fractions=graph.desc_contribution[desc_start_gid:desc_end_gid]/ \
             np.sum(graph.desc_contribution[desc_start_gid:desc_end_gid])
          # The main descendant is the one that inherits the greatest contribution
          desc_main_sid=graph.desc_IDs_gid[desc_start_gid+np.argmax(fractions)]-halo_offset
          halo.desc_main_sid=desc_main_sid
          # All orphans gals go to main descendant so increase relevant orphan count
          halos_this_snap[desc_main_sid].n_orphan+=halo.n_orphan 
          # Now loop over descendants transferring properties to them:
          for i_desc in range(desc_start_gid,desc_end_gid):
             desc_halo=halos_this_snap[graph.desc_IDs_gid[i_desc]-halo_offset]
             # Distribute mass to descendants in proportion to fractional contributions
             i_frac=i_desc-desc_start_gid # fraction index corresponding to descendent index i_desc
             desc_halo.mass_from_progenitors+=fractions[i_frac]*halo.mass
	     # ... repeat for other baryonic properties ...

Next we loop over subhalos.  For now the main descendent subhalo gets everything.  If there is no descendent then the hot gas and galaxies get given to the descendent of the host halo.

.. code-block:: python3	  

    # Now loop over the subhalos
    if subs_last_snap != None:
        for sub in subs_last_snap:
            sub_desc_start_gid=sub.desc_start_gid
            sub_desc_end_gid=sub.desc_end_gid
            host_sid=sub.host-halo_offset_last
            desc_main_sid=halos_last_snap[host_sid].desc_main_sid
            if sub.n_desc==0:
                # If no descendant, subhalo components get given to the (main descendant of) the host halo
                # and gals become orphans of that halo.  So add to relevant orphan count.
                halos_this_snap[desc_main_sid].n_orphan+=sub.n_gal
	        # ... add subhalo baryons to the descendent halo ...
            else:
                # Otherwise the main subhalo descendant gets all the gals and hot gas - 
                # i.e. assume that subhalos cannot split.
                fractions=graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid]/ \
                    np.sum(graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid])
                sub_desc_main_sid=graph.sub_desc_IDs_gid[sub_desc_start_gid+np.argmax(fractions)]-sub_offset
                sub.desc_main_sid=sub_desc_main_sid
                subs_this_snap[sub_desc_main_sid].n_gal+=sub.n_gal
                # ... add subhalo baryons to the main descendent subhalo ...

Next we count the total number of galaxies and initialise the galaxy numpy array:

.. code-block:: python3

    # Now count the total number of gals and generate the gal array.
    # This is done as a loop over subhalos within halos so as to keep all gals in a halo 
    # closely associated in the array.
    n_gal=0
    for halo in halos_this_snap:
        n_gal_start=n_gal
        if halo.n_sub>0:
            for sub in subs_this_snap[halo.sub_start_sid:halo.sub_end_sid]:
                # Record the location of this subhalo's gals in the gal lookup table.  This also updates n_gal.
                n_gal=sub.gal_loc(n_gal)
        # Record the starting location of all this halos gals, and of of its orphans, in the gal lookup table, 
        # and update n_gal to include the orphans.
        n_gal=halo.gal_loc(n_gal_start,n_gal)
    if n_gal==0: return None
    # Create new gal array and initially set all entries to empty and existence to True
    gals_this_snap=np.empty(n_gal,dtype=D_gal)
    gals_this_snap[:]=gal_template
    # Set galaxy gids and update graph galaxy counter (in that order).
    gals_this_snap['gal_gid']=graph.n_gal+np.arange(n_gal)
    graph.n_gal+=n_gal

Now we do a second pass to populate galaxies with inherited properties:

.. code-block:: python3

    # Second loop to pass on gal properties.
    if gals_last_snap != None:
        if parameters.b_debug: 
            print('Pushing gals',flush=True)
        for halo in halos_last_snap:
            if halo.b_desc_exists == False: continue
            n_orphan=halo.n_orphan
            if n_orphan > 0:
                # match up orphans
                desc_halo=halos_this_snap[halo.desc_main_sid]
                # The is the location of orphan galaxies in the previous snapshot
                gal_last_start_sid=halo.orphan_start_sid
                gal_last_end_sid=gal_last_start_sid+n_orphan
                # and in the current snapshot
                gal_this_start_sid=desc_halo.orphan_count(n_orphan)
                gal_this_end_sid=gal_this_start_sid+n_orphan
                # Copy over all properties
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                # Update the tracking pointers
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=desc_halo.halo_gid
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=desc_halo.halo_gid-halo_offset
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=parameters.NO_DATA_INT
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=parameters.NO_DATA_INT
                # Inherited orphans will not have merged (I think); otherwise the following line could be overwritten
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT
        if subs_last_snap != None:
            for sub in subs_last_snap:
                n_sub_gal=sub.n_gal
                sub_desc_start_gid=sub.desc_start_gid
                sub_desc_end_gid=sub_desc_start_gid+sub.n_desc
                gal_last_start_sid=sub.gal_start_sid
                gal_last_end_sid=gal_last_start_sid+n_sub_gal
                if sub.n_desc==0:
                    # If no descendant gals become orphans of (the main descendant of) the host halo
                    desc_halo=halos_this_snap[sub.desc_host_sid]
                    gal_this_start_sid=desc_halo.orphan_count(n_sub_gal)
                    gal_this_end_sid=gal_this_start_sid+n_sub_gal
                    # Copy over all properties
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                    # Update the tracking pointers
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=desc_halo.halo_gid
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=desc_halo.halo_gid-halo_offset
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=parameters.NO_DATA_INT
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=parameters.NO_DATA_INT
                    # New orphans will not have merged (I think); otherwise the following line could be overwritten
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT
                else:
                    # Otherwise the main subhalo descendant gets all the gals
                    desc_sub=subs_this_snap[sub.desc_main_sid]
                    desc_halo=halos_this_snap[sub.desc_host_sid]
                    # Obtain current galaxy counter for this subhalo
                    gal_this_start_sid=desc_sub.gal_count(n_sub_gal)
                    gal_this_end_sid=gal_this_start_sid+n_sub_gal
                    # Copy over all properties
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                    # Update the tracking pointers
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=desc_halo.halo_gid
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=desc_halo.halo_gid-halo_offset
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=desc_sub.sub_gid
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=desc_sub.sub_gid-sub_offset
                    # This is probably wrong: we need to check if there is already an entry for
                    # first_prog_sid for these galaxies and, if so, update next_prog_sid to point to it.
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                    gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT
    gals_this_snap['graph_ID']=graph.graph_ID
    gals_this_snap['snap_ID']=halos_this_snap[0].snap_ID

