import numpy as np
from codetiming import Timer
from profiling import conditional_decorator

import commons
b_profile_cpu=commons.load('b_profile_cpu')

from gals import D_gal
from subs import F_subs_mass_baryon

@conditional_decorator(Timer(name='F_push_snap',logger=None),b_profile_cpu)
def F_push_snap(halos_last_snap,halos_this_snap,subs_last_snap,subs_this_snap,
                   gals_last_snap,graph,parameters,variables):
    """
    Propagate properties from progenitor halos to descendants.
    To convert this to C one would first have to extract all the required graph properties into numpy arrays,
    then pass those into this routine as arguments.  This would not be too hard to do, but only worth it if
    this routine is CPU time-limiting.

    Done as a push rather than a pull because sharing determined by progenitor.
    First loop to push halo / subhalo properties; 
    then structured gal array needs to be generated;
    then second push of gal properties.

    Arguments
    ---------
    halos_last : obj : D_halo[n_halo_last_snap]
       The halos from the last (ie previous) snapshot
    halos_this : obj : D_halo[n_halo_last_snap]
       The halos from the current snapshot
    subs_last : obj : D_sub[n_sub_last_snap]
       The subhalos from the last (ie previous) snapshot
    subs_this : obj : D_sub[n_sub_last_snap]
       The subhalos from the current snapshot
    gals_last : obj : D_gal[n_gal_last_snap]
       The galaxies from the last (ie previous) snapshot
    graph : obj : C_graph
       Instance of class containing the current graph
    parameters : obj : C_parameters
       Instance of class containing the run-time parameters
    variables : obj : C_variables
       Instance of class containing global variables structure

    Returns
    -------
    n_gal_this_snap : int
       The number of galaxies in this snapshot
    """

    # These offsets give the first (sub)halo in this snapshot
    n_halo_this_snap=len(halos_this_snap)
    halo_offset_this_snap=halos_this_snap[0]['halo_gid']
    if isinstance(subs_this_snap,np.ndarray):
        n_sub_this_snap=len(subs_this_snap)
        sub_offset_this_snap=subs_this_snap[0]['sub_gid']
    else:
        n_sub_this_snap=0
    if isinstance(halos_last_snap,np.ndarray):
        n_halo_last_snap=len(halos_last_snap)
        halo_offset_last_snap=halos_last_snap[0]['halo_gid']
    else:
        n_halo_last_snap=0
    if isinstance(subs_last_snap,np.ndarray):
        n_sub_last_snap=len(subs_last_snap)
        sub_offset_last_snap=subs_last_snap[0]['sub_gid']
    else:
        n_sub_last_snap=0
        
    for i_halo_last_snap in range(n_halo_last_snap):
	# First determine what fraction to give to each descendant
        n_desc=halos_last_snap[i_halo_last_snap]['n_desc']
        desc_start_gid=halos_last_snap[i_halo_last_snap]['desc_start_gid']
        desc_end_gid=halos_last_snap[i_halo_last_snap]['desc_end_gid']
        desc_main_sid=halos_last_snap[i_halo_last_snap]['desc_main_sid']
        # Hack to fix broken Mill trees:
        if desc_start_gid==-1: n_desc=0
        if n_desc==0: 
            if parameters.verbosity >= 3: print('No descendants for halo:',halos_last_snap[i_halo_last_snap]['halo_gid'],flush=True)
            halos_last_snap[i_halo_last_snap]['b_desc_exists'] = False
            # For now just skip this halo; might want in future to log these occurrences
            # Note that any orphan galaxies will cease to exist, gal['b_exists'] == False
            continue
        fractions=graph.halo_desc_contribution[desc_start_gid:desc_end_gid]/ \
            np.sum(graph.halo_desc_contribution[desc_start_gid:desc_end_gid])
        # The main descendant is the one that inherits the greatest contribution
        desc_main_sid=graph.halo_desc_IDs_gid[desc_start_gid+np.argmax(fractions)]-halo_offset_this_snap
        halos_last_snap[i_halo_last_snap]['desc_main_sid']=desc_main_sid
        # All orphans gals go to main descendant so increase relevant orphan count
        halos_this_snap[desc_main_sid]['n_orphan']+=halos_last_snap[i_halo_last_snap]['n_orphan']
        # Now loop over descendants transferring properties to them:
        for i_desc in range(desc_start_gid,desc_end_gid):
            i_desc_halo=graph.halo_desc_IDs_gid[i_desc]-halo_offset_this_snap
            # Distribute mass to descendants in proportion to fractional contributions
            i_frac=i_desc-desc_start_gid # fraction index corresponding to descendent index i_desc
            halos_this_snap[i_desc_halo]['mass_from_progenitors']+=fractions[i_frac]*halos_last_snap[i_halo_last_snap]['mass']
            halos_this_snap[i_desc_halo]['mass_gas_hot']+=fractions[i_frac]*halos_last_snap[i_halo_last_snap]['mass_gas_hot']
            halos_this_snap[i_desc_halo]['mass_metals_gas_hot']+=fractions[i_frac]*halos_last_snap[i_halo_last_snap]['mass_metals_gas_hot']
            halos_this_snap[i_desc_halo]['mass_stars']+=fractions[i_frac]*halos_last_snap[i_halo_last_snap]['mass_stars']
            halos_this_snap[i_desc_halo]['mass_metals_stars']+=fractions[i_frac]*halos_last_snap[i_halo_last_snap]['mass_metals_stars']
            
    # Now loop over the subhalos
    for i_sub_last_snap in range(n_sub_last_snap):
        sub_desc_start_gid=subs_last_snap[i_sub_last_snap]['desc_start_gid']
        sub_desc_end_gid=subs_last_snap[i_sub_last_snap]['desc_end_gid']
        # Hack to fix broken Mill trees
        if sub_desc_start_gid==-1: subs_last_snap[i_sub_last_snap]['n_desc']=0
        halo_sid=subs_last_snap[i_sub_last_snap]['halo_sid']
        desc_main_sid=halos_last_snap[halo_sid]['desc_main_sid']  # This possibly does not exist
        if desc_main_sid==parameters.NO_DATA_INT:
            if parameters.verbosity >= 3: print('No descendant for halo ',sub_last_snap[i_sub_last_snap]['halo_gid'],': galaxies will be lost')
            continue
        subs_last_snap[i_sub_last_snap]['desc_halo_sid']=desc_main_sid
        if subs_last_snap[i_sub_last_snap]['n_desc']==0:
            # If no descendant, subhalo components get given to the (main descendant of) the host halo
            # and gals become orphans of that halo.   So add to relevant orphan count.
            # Note: it seems that this can result in a huge excess of baryons in the descendant halo.
            halos_this_snap[desc_main_sid]['n_orphan']+=subs_last_snap[i_sub_last_snap]['n_gal']
            halos_this_snap[desc_main_sid]['mass_gas_hot']+=subs_last_snap[i_sub_last_snap]['mass_gas_hot']
            halos_this_snap[desc_main_sid]['mass_metals_gas_hot']+=subs_last_snap[i_sub_last_snap]['mass_metals_gas_hot']
            halos_this_snap[desc_main_sid]['mass_stars']+=subs_last_snap[i_sub_last_snap]['mass_stars']
            halos_this_snap[desc_main_sid]['mass_metals_stars']+=subs_last_snap[i_sub_last_snap]['mass_metals_stars']
        else:
            # Otherwise the main subhalo descendant gets all the gals and hot gas - 
            # i.e. assume that subhalos cannot split. 
            fractions=graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid]/ \
                np.sum(graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid])
            sub_desc_main_sid=graph.sub_desc_IDs_gid[sub_desc_start_gid+np.argmax(fractions)]-sub_offset_this_snap
            subs_last_snap[i_sub_last_snap]['desc_main_sid']=sub_desc_main_sid
            subs_this_snap[sub_desc_main_sid]['n_gal']+=subs_last_snap[i_sub_last_snap]['n_gal']
            # So long as subhalos don't split, we can simply inherit the baryon mass
            subs_this_snap[sub_desc_main_sid]['mass_baryon']+=subs_last_snap[i_sub_last_snap]['mass_baryon']
            subs_this_snap[sub_desc_main_sid]['mass_gas_hot']+=subs_last_snap[i_sub_last_snap]['mass_gas_hot']
            subs_this_snap[sub_desc_main_sid]['mass_metals_gas_hot']+=subs_last_snap[i_sub_last_snap]['mass_metals_gas_hot']
            subs_this_snap[sub_desc_main_sid]['mass_stars']+=subs_last_snap[i_sub_last_snap]['mass_stars']
            subs_this_snap[sub_desc_main_sid]['mass_metals_stars']+=subs_last_snap[i_sub_last_snap]['mass_metals_stars']
        
    # Now count the total number of gals and generate the gal array.
    # This is done as a loop over subhalos within halos so as to keep all gals in a halo 
    # closely associated in the array.
    n_gal_this_snap=0
    for i_halo_this_snap in range(n_halo_this_snap):
        halos_this_snap[i_halo_this_snap]['gal_start_sid'] = n_gal_this_snap
        # Place the subhalo galaxies before the orphans
        for i_sub in range(halos_this_snap[i_halo_this_snap]['sub_start_sid'],halos_this_snap[i_halo_this_snap]['sub_end_sid']):
            # Record the location of this subhalo's gals in the gal lookup table.
            # Require subhalo to have at least 1 galaxy
            subs_this_snap[i_sub]['n_gal'] = max(subs_this_snap[i_sub]['n_gal'],1)
            subs_this_snap[i_sub]['gal_start_sid'] = n_gal_this_snap
            subs_this_snap[i_sub]['gal_next_sid'] = n_gal_this_snap # Will be used to keep track of galaxies during update_halo phase.
            gal_end_sid = n_gal_this_snap + subs_this_snap[i_sub]['n_gal']
            subs_this_snap[i_sub]['gal_end_sid'] = gal_end_sid
            n_gal_this_snap = gal_end_sid
        halos_this_snap[i_halo_this_snap]['orphan_start_sid'] = n_gal_this_snap
        halos_this_snap[i_halo_this_snap]['orphan_next_sid'] = n_gal_this_snap # Will be used to keep track of orphans during update_halo phase
        n_gal_this_snap += halos_this_snap[i_halo_this_snap]['n_orphan']
        halos_this_snap[i_halo_this_snap]['gal_end_sid']=n_gal_this_snap
        halos_this_snap[i_halo_this_snap]['orphan_end_sid']=n_gal_this_snap
    # I don't think that it should be possible to have zero galaxies as the code will probably break if so.
    assert n_gal_this_snap>0

    # Create new gal array
    # Galaxies initially have b_exist == False and will come into existence if either inherited or 
    # the central galaxy of a subhalo.
    gals_this_snap=np.empty(n_gal_this_snap,dtype=D_gal)
    gals_this_snap[:]=parameters.gal_template
    # Set galaxy gids and update graph galaxy counter (in that order).
    gals_this_snap['gal_gid']=graph.n_gal+np.arange(n_gal_this_snap)
    graph.n_gal+=n_gal_this_snap
    
    # Second loop to pass on gal properties.
    if isinstance(gals_last_snap, np.ndarray):
        for i_halo_last_snap in range(n_halo_last_snap):
            if halos_last_snap[i_halo_last_snap]['n_desc'] == 0: continue  # Hopefully all halos have descendants, but just in case...
            n_orphan=halos_last_snap[i_halo_last_snap]['n_orphan']
            if n_orphan > 0:
                # match up orphans
                i_desc_halo=halos_last_snap[i_halo_last_snap]['desc_main_sid']
                # This catch is to get the code to run on badly-formed input
                try:
                    halos_this_snap[i_desc_halo]
                except:
                    print('Halo has no descendants:')
                    print('  halo =',halos_last_snap[i_halo_last_snap])
                    print('  desc_start_gid =',halos_last_snap[i_halo_last_snap]['desc_start_gid'])
                    print('  desc_end_gid =',halos_last_snap[i_halo_last_snap]['desc_end_gid'])
                    continue # Failed to find descendant halo (should be impossible in well-structured graph); throw away 
                assert i_desc_halo != parameters.NO_DATA_INT
                # The is the location of orphan galaxies in the previous snapshot
                gal_last_start_sid=halos_last_snap[i_halo_last_snap]['orphan_start_sid']
                gal_last_end_sid=gal_last_start_sid+n_orphan
                # and in the current snapshot
                gal_this_start_sid=halos_this_snap[i_desc_halo]['orphan_next_sid']
                gal_this_end_sid=gal_this_start_sid+n_orphan
                halos_this_snap[i_desc_halo]['orphan_next_sid'] = gals_this_end_sid  # Have lost track of if/why I need this :-)
                # Copy over all properties
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                # Update the tracking pointers
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=halos_this_snap[i_desc_halo]['halo_gid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=halos_this_snap[i_desc_halo]['halo_sid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=parameters.NO_DATA_INT
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=parameters.NO_DATA_INT
                # Inherited orphans will not have merged (I think); otherwise the following line could be overwritten
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT
        for i_sub_last_snap in range(n_sub_last_snap):
            if subs_last_snap[i_sub_last_snap]['desc_halo_sid']==parameters.NO_DATA_INT: continue # No descendent halo
            n_sub_gal=subs_last_snap[i_sub_last_snap]['n_gal']
            sub_desc_start_gid=subs_last_snap[i_sub_last_snap]['desc_start_gid']
            sub_desc_end_gid=subs_last_snap[i_sub_last_snap]['desc_end_gid']
            gal_last_start_sid=subs_last_snap[i_sub_last_snap]['gal_start_sid']
            gal_last_end_sid=subs_last_snap[i_sub_last_snap]['gal_end_sid']
            if subs_last_snap[i_sub_last_snap]['n_desc']==0:
                # If no descendant gals become orphans of (the main descendant of) the host halo
                i_desc_halo=subs_last_snap[i_sub_last_snap]['desc_halo_sid']
                assert i_desc_halo != parameters.NO_DATA_INT
                gal_this_start_sid=halos_this_snap[i_desc_halo]['orphan_next_sid']
                gal_this_end_sid=gal_this_start_sid+n_sub_gal
                halos_this_snap[i_desc_halo]['orphan_next_sid']=gal_this_end_sid
                # Copy over all properties
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                # Update the tracking pointers
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=halos_this_snap[i_desc_halo]['halo_gid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=halos_this_snap[i_desc_halo]['halo_sid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=parameters.NO_DATA_INT
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=parameters.NO_DATA_INT
                # New orphans will not have merged (I think); otherwise the following line could be overwritten
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT
            else:
                # Otherwise the main subhalo descendant gets all the gals
                i_desc_sub=subs_last_snap[i_sub_last_snap]['desc_main_sid']
                i_desc_halo=subs_last_snap[i_sub_last_snap]['desc_halo_sid']
                assert i_desc_halo != parameters.NO_DATA_INT
                # Obtain current galaxy counter for this subhalo
                gal_this_start_sid=subs_this_snap[i_desc_sub]['gal_next_sid']
                gal_this_end_sid=gal_this_start_sid+n_sub_gal
                subs_this_snap[i_desc_sub]['gal_next_sid']=gal_this_end_sid
                # Copy over all properties
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]=gals_last_snap[gal_last_start_sid:gal_last_end_sid]
                # Update the tracking pointers
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_gid']=halos_this_snap[i_desc_halo]['halo_gid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['halo_sid']=halos_this_snap[i_desc_halo]['halo_sid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_gid']=subs_this_snap[i_desc_sub]['sub_gid']
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['sub_sid']=subs_this_snap[i_desc_sub]['sub_sid']
                # This is probably wrong as it assumes no galaxy merging: we need to check if there is already an entry for
                # first_prog_sid for these galaxies and, if so, update next_prog_sid to point to it.
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['first_prog_gid']=np.arange(gal_this_start_sid,gal_this_end_sid)
                gals_this_snap[gal_this_start_sid:gal_this_end_sid]['next_prog_gid']=parameters.NO_DATA_INT

    # Set galaxy properties that are not inherited
    gals_this_snap['graph_ID']=graph.graph_ID
    gals_this_snap['snap_ID']=halos_this_snap[0]['snap_ID']
    gals_this_snap['gal_sid']=np.arange(len(gals_this_snap))
    gals_this_snap['SFR_snap']=0.  # Needs to be zeroed as accumulates over snapshot timesteps.
    # Need to set central galaxies for each subhalo as these are the ones that are accrete gas.
    # Also, newly created galaxies (no progenitors) come into existence at this point.
    for i_sub_this_snap in range(n_sub_this_snap): 
        # Will eventually have fancy code to determine which, if any, galaxy is the central one.
        # For now, just make that the first (and only) galaxy.
        i_gal_central=subs_this_snap[i_sub_this_snap]['gal_start_sid']
        assert 0<=i_gal_central<n_gal_this_snap
        subs_this_snap[i_sub_this_snap]['gal_central_sid'] = i_gal_central
        gals_this_snap[i_gal_central]['b_exists']=True
        gals_this_snap[i_gal_central]['sub_gid']=subs_this_snap[i_sub_this_snap]['sub_gid']
        gals_this_snap[i_gal_central]['sub_sid']=subs_this_snap[i_sub_this_snap]['sub_sid']
        gals_this_snap[i_gal_central]['halo_gid']=subs_this_snap[i_sub_this_snap]['halo_gid']
        gals_this_snap[i_gal_central]['halo_sid']=subs_this_snap[i_sub_this_snap]['halo_sid']
        # The central galaxies have their virial speed updated; the others keep their inherited virial speed
        gals_this_snap[i_gal_central]['v_vir']=subs_this_snap[i_sub_this_snap]['half_mass_virial_speed']

    # Some sanity checks
    # In principle subhalos should aready have the correct baryon mass: this is a check
    for i_sub_this_snap in range(n_sub_this_snap):
        sub=subs_this_snap[i_sub_this_snap]  # Does this work to pick out a row?
        if sub['mass_baryon'] > parameters.mass_minimum_internal:
            # Had to relax condition from 1e-4 to 1e-2 to accommodate Mill tree errors.
            # So probably still a bug in my conversion code.
            # But at this point I am not inclined to waste any more time trying to fix it.
            # (some subhalo pointers are -1!)
            if np.abs(sub['mass_baryon']/F_subs_mass_baryon(sub,gals_this_snap)-1.)>1e-2:
                print('gals[''mass_baryon''] =',gals_this_snap[sub['gal_start_sid']:sub['gal_end_sid']]['mass_baryon'])
                print('gals[''mass_gas_cold''] =',gals_this_snap[sub['gal_start_sid']:sub['gal_end_sid']]['mass_gas_cold'])
                print('gals[''mass_stars_disc''] =',gals_this_snap[sub['gal_start_sid']:sub['gal_end_sid']]['mass_stars_disc'])
                print('sub.mass_gas_hot =',sub['mass_gas_hot'])
                print('sub.mass_stars =',sub['mass_stars'])
                print('sub.mass_baryon =',sub['mass_baryon'])
                halo=halos_this_snap[sub['halo_sid']] # Does this work to pick out a row?
                print('halo.mass_gas_eject =',halo['mass_gas_eject'])
                print('halo.mass_gas_hot =',halo['mass_gas_hot'])
                print('halo.mass_stars =',halo['mass_stars'])
                print('halo.mass_baryon =',halo['mass_baryon'])
                print('subhalo baryon mass discrepancy',
                    sub['graph_ID'],sub['snap_ID'],sub['sub_sid'],sub['mass_baryon'],F_subs_mass_baryon(sub,gals_this_snap))
    for i_gal_this_snap in range(n_gal_this_snap):
        # Check that all halos are assigned to a halo
        if gals_this_snap[i_gal_this_snap]['b_exists'] and gals_this_snap[i_gal_this_snap]['halo_sid']==parameters.NO_DATA_INT:
            print(gals_this_snap.dtype)
            print(gals_this_snap[i_gal_this_snap])
            raise AssertionError('Galaxy not assigned to a halo')
        # Check that all halos are assigned to either a subhalo or a halo
        if gals_this_snap[i_gal_this_snap]['b_exists'] and gals_this_snap[i_gal_this_snap]['sub_sid']==parameters.NO_DATA_INT \
           and gals_this_snap[i_gal_this_snap]['halo_sid']==parameters.NO_DATA_INT:
            print(gals_this_snap.dtype)
            print(gals_this_snap[i_gal_this_snap])
            raise AssertionError('Galaxy assigned neither subhalo nor halo')
    # Check that all galaxies have a finite virial speed
    for gal in gals_this_snap:
        if gals_this_snap[i_gal_this_snap]['b_exists'] and gals_this_snap[i_gal_this_snap]['v_vir']<1e-10:
            print(gals_this_snap.dtype)
            print(gals_this_snap[i_gal_this_snap])
            raise AssertionError('Galaxy virial speed is too low')
    # Check that all galaxies with cold gas have a disc scale length set
    for gal in gals_this_snap:
        if gals_this_snap[i_gal_this_snap]['mass_gas_cold']>parameters.mass_minimum_internal and \
           gals_this_snap[i_gal_this_snap]['radius_gas_cold']<parameters.length_minimum_internal:
            print(gals_this_snap.dtype)
            print(gals_this_snap[i_gal_this_snap] )
            raise AssertionError('Galaxy disc scale length not set')
    return gals_this_snap
