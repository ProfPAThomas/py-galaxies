"""**High-level driver routines.**
"""

import numpy as np
from codetiming import Timer
from profiling import conditional_decorator

import commons
b_profile_cpu=commons.load('b_profile_cpu')
b_SFH=commons.load('b_SFH')

# Import astrophysics modules (could be run-time parameter dependent)
from cooling import F_halo as F_halo_cooling
from cooling import F_sub as F_sub_cooling
from gals import D_gal, F_gal_template
from star_formation_and_feedback import F_gal_form_stars, F_gal_SNR_feedback
from mergers import F_merge_gals as F_merge_gals
if b_SFH: from sfh import F_sfh_update_bins

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_halo_set_baryon_fraction',logger=None),b_profile_cpu)
def F_halo_set_baryon_fraction(halo,parameters):
    """
    Updates the baryon content to be the universal mean, or the sum of the baryon content from
    the progenitors, whichever is larger (so that baryons are not lost).

    Any excess baryons arrive in the form of base_metallicity hot gas.

    Arguments
    ---------
    halo : obj : C_halo
       The halo currently being processed.
    parameters : obj : C_parameters
       Instance of class containing global parameters

    Returns
    -------
    None
    """  
    delta_baryon=max(0.,parameters.baryon_fraction*max(halo.mass,halo.mass_from_progenitors)-halo.mass_baryon)
    halo.mass_baryon+=delta_baryon
    halo.mass_gas_hot+=delta_baryon
    halo.mass_metals_gas_hot+=delta_baryon*parameters.base_metallicity
    return None

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_halo_reincorporation',logger=None),b_profile_cpu)
def F_halo_reincorporation(halo,parameters):
    """
    Reincorporation of ejected gas.

    Currently just assumes Hen15 model.
    Might be better to pass dt_halo and c_reinc explicitly

    Arguments
    ---------
    halo : obj : C_halo
       The halo currently being processed.
    parameters : obj : C_parameters
       Instance of class containing global parameters

    Returns
    -------
    None
    """
    dt_halo=commons.load('dt_halo')
    
    t_reinc = parameters.c_Hen15_reinc/halo.mass
    mass_reinc = halo.mass_gas_eject * (1.-np.exp(-dt_halo/t_reinc))
    mass_metals_reinc = mass_reinc * (halo.mass_metals_gas_eject/halo.mass_gas_eject)
    halo.mass_gas_eject -= mass_reinc
    halo.mass_metals_gas_eject -= mass_metals_reinc
    halo.mass_gas_hot += mass_reinc
    halo.mass_metals_gas_hot += mass_metals_reinc
    return None

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_process_halos',logger=None),b_profile_cpu)
def F_process_halos(halos,subs,gals,graph,parameters):
    """
    This is the controlling routine for halo processing.

    Arguments
    ---------
    halos : obj : C_halo[n_halo]
       The halos currently being processed in this graph/snapshot
    halos : obj : C_sub[n_sub]
       The subhalos currently being processed in this graph/snapshot
    gals : obj : np.darray[n_gal]
       The galaxies currently being processed in this graph/snapshot.
    parameters : obj : C_parameters
       Instance of class containing global parameters

    Returns
    -------
    None
    """

    # Set flag for existencs (or otherwise) of galaxies.
    if isinstance(gals, np.ndarray):
        b_gals_exist = True
    else:
        b_gals_exist = False

    # Load the number of timesteps
    n_dt_halo=commons.load('n_dt_halo')
    n_dt_gal=commons.load('n_dt_gal')
    if b_SFH:
        sfh=parameters.sfh
        i_dt=commons.load('i_dt')
    
    for halo in halos:
        if parameters.verbosity>=4: print('Processing halo ',halo.halo_gid)
        if halo.b_done==True:
            raise RuntimeError('halo '+str(halo.halo_gid)+' in graph '+str(halo.graph_ID)+' already processed.')
        # Accretion onto halos.
        # First determine current baryon mass, then call function to set to universal value.
        # The rate could be calculated once at the beginning of a step, but probably very quick
        halo.set_mass_baryon(subs,gals)
        F_halo_set_baryon_fraction(halo,parameters)
        # Reincorporation of ejected gas
        if halo.mass_gas_eject > parameters.mass_minimum_internal: F_halo_reincorporation(halo,parameters)
        # Cooling of gas from halo onto central subhalo (or, in L-Galaxies mode, the most massive subhalo)
        # Cooling occurs only if a central subhalo exists.
        if halo.sub_central_sid != parameters.NO_DATA_INT: 
            sub_central=subs[halo.sub_central_sid]
            F_halo_cooling(halo,sub_central,parameters)
            # In l-galaxies mode the virial velocity of the subhalo may have changed, so need to reset that of the central galaxy also
            if parameters.b_lgalaxies:
                if sub_central.gal_central_sid != parameters.NO_DATA_INT:
                    gals[sub_central.gal_central_sid]['v_vir']=sub_central.half_mass_virial_speed
        halo.n_dt+=1
        if halo.n_dt==n_dt_halo: halo.b_done=True
    if subs != None:
        for sub in subs:
            if sub.b_done==True:
                raise RuntimeError('subhalo '+str(sub.sub_gid)+' in graph '+str(sub.graph_ID)+' already processed.')
            if sub.n_gal>1:
                # Initially assume instantaneous merging of galaxies in subhalos
                F_merge_gals(halos[sub.halo_sid],sub,gals[sub.gal_start_sid:sub.gal_end_sid],parameters)
            # Not all subhalos may have hot gas
            if sub.mass_gas_hot > parameters.mass_minimum_internal:
                # debug code
                gal=gals[sub.gal_central_sid]
                if gal['b_exists']==False:
                    print('gal.dtype =',gal.dtype)
                    print('gal =',gal)
                    print('sub.sub_gid, sub.sub_sid =',sub.sub_gid,sub.sub_sid)
                    print('sub.halo_gid, sub.halo_sid =',sub.halo_gid,sub.halo_sid)
                    print('sub.gal_start_sid, sub.gal_end_sid =',sub.gal_start_sid, sub.gal_end_sid)
                    raise AssertionError('Subhalo central galaxy does not exist')
                # Cooling of hot gas in subhalo onto galaxy
                # This also includes radio mode BH growth and feedback
                F_sub_cooling(sub,gal,parameters)
                pass
            sub.n_dt+=1
            if sub.n_dt==n_dt_halo: sub.b_done=True
    if b_gals_exist:
        for i_dt_gal in range(n_dt_gal):
            for gal in gals:
                if not gal['b_exists']: continue  #  Galaxies may have merged
                gal['SFR_dt'] = 0. # This will fail to capture mergers (done above in subs loop) until we have a proper merger time for them.
                if gal['mass_gas_cold'] > parameters.mass_minimum_internal: 
                    mass_stars = F_gal_form_stars(gal,parameters)
                    # If subhalo does not exist, use halo as proxy.  This will work here as only need access to hot gas phase.
                    sub_sid=gal['sub_sid']
                    halo_sid=gal['halo_sid']
                    if sub_sid==parameters.NO_DATA_INT:
                        F_gal_SNR_feedback(mass_stars,gal,halos[halo_sid],halos[halo_sid],parameters)
                    else:
                        F_gal_SNR_feedback(mass_stars,gal,subs[sub_sid],halos[halo_sid],parameters)
            if b_SFH:
                F_sfh_update_bins(gals,sfh,parameters)
                i_dt+=1
                commons.save('i_dt',i_dt)
    return None

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_set_central_galaxy',logger=None),b_profile_cpu)
def F_set_central_galaxy(sub,parameters):
    """
    Determine the central galaxy in a subhalo.

    Will eventually have fancy code to determine which, if any, galaxy is the central one.
    For now, just make that the first (and only) galaxy.

    Arguments
    ---------
    sub : obj : C_sub
       The subhalo currently being processed.
    parameters : obj : C_parameters
       Instance of class containing global parameters
    """
    
    sub.gal_central_sid = sub.gal_start_sid
    return None

#------------------------------------------------------------------------------------------------------

@conditional_decorator(Timer(name='F_update_halos',logger=None),b_profile_cpu)
def F_update_halos(halos_last_snap,halos_this_snap,subs_last_snap,subs_this_snap,
                   gals_last_snap,graph,parameters):
    """
    Propagate properties from progenitor halos to descendants.

    Done as a push rather than a pull because sharing determined by progenitor.
    First loop to push halo / subhalo properties; 
    then structured gal array needs to be generated;
    then second push of gal properties.

    Arguments
    ---------
    halos : obj : C_halo[n_halo]
       The halos currently being processed in this graph/snapshot
    halos : obj : C_sub[n_sub]
       The subhalos currently being processed in this graph/snapshot
    gals : obj : np.darray[n_gal]
       The galaxies currently being processed in this graph/snapshot.
    parameters : obj : C_parameters
       Instance of class containing global parameters

    Returns
    -------
    None
    """
    
    # These offsets give the first (sub)halo in this snapshot
    halo_offset=halos_this_snap[0].halo_gid
    if subs_this_snap != None: sub_offset=subs_this_snap[0].sub_gid
    if halos_last_snap != None: halo_offset_last=halos_last_snap[0].halo_gid
    if subs_last_snap != None: sub_offset_last=subs_last_snap[0].sub_gid
    if halos_last_snap != None:
        for halo in halos_last_snap:
            # First determine what fraction to give to each descendant
            desc_start_gid=halo.desc_start_gid
            desc_end_gid=halo.desc_end_gid
            if (halo.n_desc==0): 
                if parameters.verbosity >= 4: print('No descendants for halo:',halo,flush=True)
                halo.b_desc_exists = False
                # For now just skip this halo; might want in future to log these occurrences
                # Note that any orphan galaxies will cease to exist, gal['b_exists'] == False
                continue
            fractions=graph.halo_desc_contribution[desc_start_gid:desc_end_gid]/ \
                np.sum(graph.halo_desc_contribution[desc_start_gid:desc_end_gid])
            # The main descendant is the one that inherits the greatest contribution
            desc_main_sid=graph.halo_desc_IDs_gid[desc_start_gid+np.argmax(fractions)]-halo_offset
            halo.desc_main_sid=desc_main_sid
            # assert desc_main_sid < parameters.n_graph
            # All orphans gals go to main descendant so increase relevant orphan count
            halos_this_snap[desc_main_sid].n_orphan+=halo.n_orphan
            # Now loop over descendants transferring properties to them:
            for i_desc in range(desc_start_gid,desc_end_gid):
                desc_halo=halos_this_snap[graph.halo_desc_IDs_gid[i_desc]-halo_offset]
                # assert desc_halo_gid == desc_halo.halo_gid
                if parameters.verbosity>=5: print('Processing descendant',desc_halo.halo_gid)
                # Distribute mass to descendants in proportion to fractional contributions
                i_frac=i_desc-desc_start_gid # fraction index corresponding to descendent index i_desc
                desc_halo.mass_from_progenitors+=fractions[i_frac]*halo.mass
                desc_halo.mass_gas_hot+=fractions[i_frac]*halo.mass_gas_hot
                desc_halo.mass_metals_gas_hot+=fractions[i_frac]*halo.mass_metals_gas_hot
                desc_halo.mass_stars+=fractions[i_frac]*halo.mass_stars
                desc_halo.mass_metals_stars+=fractions[i_frac]*halo.mass_metals_stars
            
    # Now loop over the subhalos
    if subs_last_snap != None:
        for sub in subs_last_snap:
            sub_desc_start_gid=sub.desc_start_gid
            sub_desc_end_gid=sub.desc_end_gid
            halo_sid=sub.halo_sid
            desc_main_sid=halos_last_snap[halo_sid].desc_main_sid  # This possibly does not exist
            if desc_main_sid==parameters.NO_DATA_INT:
                #warn('No descendant for halo '+str(halo_sid)+': galaxies will be lost')
                print('No descendant for halo '+str(halo_sid)+': galaxies will be lost')
                continue
            sub.desc_halo_sid=desc_main_sid
            if sub.n_desc==0:
                # If no descendant, subhalo components get given to the (main descendant of) the host halo
                # and gals become orphans of that halo.   So add to relevant orphan count.
                # Note: it seems that this can result in a huge excess of baryons in the descendant halo.
                halos_this_snap[desc_main_sid].n_orphan+=sub.n_gal
                halos_this_snap[desc_main_sid].mass_gas_hot+=sub.mass_gas_hot
                halos_this_snap[desc_main_sid].mass_metals_gas_hot+=sub.mass_metals_gas_hot
                halos_this_snap[desc_main_sid].mass_stars+=sub.mass_stars
                halos_this_snap[desc_main_sid].mass_metals_stars+=sub.mass_metals_stars
            else:
                # Otherwise the main subhalo descendant gets all the gals and hot gas - 
                # i.e. assume that subhalos cannot split. 
                fractions=graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid]/ \
                    np.sum(graph.sub_desc_contribution[sub_desc_start_gid:sub_desc_end_gid])
                sub_desc_main_sid=graph.sub_desc_IDs_gid[sub_desc_start_gid+np.argmax(fractions)]-sub_offset
                sub.desc_main_sid=sub_desc_main_sid
                subs_this_snap[sub_desc_main_sid].n_gal+=sub.n_gal
                # So long as subhalos don't split, we can simply inherit the baryon mass
                subs_this_snap[sub_desc_main_sid].mass_baryon+=sub.mass_baryon
                subs_this_snap[sub_desc_main_sid].mass_gas_hot+=sub.mass_gas_hot
                subs_this_snap[sub_desc_main_sid].mass_metals_gas_hot+=sub.mass_metals_gas_hot
                subs_this_snap[sub_desc_main_sid].mass_stars+=sub.mass_stars
                subs_this_snap[sub_desc_main_sid].mass_metals_stars+=sub.mass_metals_stars
                # Now loop over descendants transferring properties to them
                # Only required if we decide that subhalos can split
                # for i_desc in range(sub_desc_start_gid,sub_desc_end_gid):
                #     desc_sub_gid=graph.sub_desc_IDs_gid[i_desc]
                #     desc_sub=subs_this_snap[desc_sub_gid-sub_gid_offset]
                #     assert desc_sub_gid == desc_sub.sub_gid
                #     sub_desc_halo.<quantity>+=fractions[i_desc-desc_start]*sub.<quantity>
        
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
    # Create new gal array
    # Galaxies initially have b_exist == False and will come into existence if either inherited or 
    # the central galaxy of a subhalo.
    gals_this_snap=np.empty(n_gal,dtype=D_gal)
    gals_this_snap[:]=parameters.gal_template
    # Set galaxy gids and update graph galaxy counter (in that order).
    gals_this_snap['gal_gid']=graph.n_gal+np.arange(n_gal)
    graph.n_gal+=n_gal
    
    # Second loop to pass on gal properties.
    if isinstance(gals_last_snap, np.ndarray):
        for halo in halos_last_snap:
            if halo.n_desc == 0: continue
            n_orphan=halo.n_orphan
            if n_orphan > 0:
                # match up orphans
                desc_halo=halos_this_snap[halo.desc_main_sid]
                assert desc_halo != parameters.NO_DATA_INT
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
                #gals_this_snap[gal_this_start_sid:gal_this_end_sid]['v_vir']=desc_halo.half_mass_virial_speed
        if subs_last_snap != None:
            for sub in subs_last_snap:
                if sub.desc_halo_sid==parameters.NO_DATA_INT: continue # No descendent halo
                n_sub_gal=sub.n_gal
                sub_desc_start_gid=sub.desc_start_gid
                sub_desc_end_gid=sub_desc_start_gid+sub.n_desc
                gal_last_start_sid=sub.gal_start_sid
                gal_last_end_sid=gal_last_start_sid+n_sub_gal
                if sub.n_desc==0:
                    # If no descendant gals become orphans of (the main descendant of) the host halo
                    desc_halo=halos_this_snap[sub.desc_halo_sid]
                    assert desc_halo != parameters.NO_DATA_INT
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
                    #gals_this_snap[gal_this_start_sid:gal_this_end_sid]['v_vir']=desc_halo.half_mass_virial_speed
                else:
                    # Otherwise the main subhalo descendant gets all the gals
                    desc_sub=subs_this_snap[sub.desc_main_sid]
                    desc_halo=halos_this_snap[sub.desc_halo_sid]
                    assert desc_halo != parameters.NO_DATA_INT
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
                    #gals_this_snap[gal_this_start_sid:gal_this_end_sid]['v_vir']=desc_sub.half_mass_virial_speed

    # Set galaxy properties that are not inherited
    gals_this_snap['graph_ID']=graph.graph_ID
    gals_this_snap['snap_ID']=halos_this_snap[0].snap_ID
    gals_this_snap['gal_sid']=np.arange(len(gals_this_snap))
    # Need to set central galaxies for each subhalo as these are the ones that are accrete gas.
    # Also, newly created galaxies (no progenitors) come into existance at this point.
    for sub in subs_this_snap: 
        F_set_central_galaxy(sub,parameters)    
        gal_central=gals_this_snap[sub.gal_central_sid]
        gal_central['b_exists']=True
        gal_central['sub_gid']=sub.sub_gid
        gal_central['sub_sid']=sub.sub_sid
        gal_central['halo_gid']=sub.halo_gid
        gal_central['halo_sid']=sub.halo_sid
        # The central galaxies have their virial speed updated; the others keep their inherited virial speed
        gals_this_snap[sub.gal_central_sid]['v_vir']=sub.half_mass_virial_speed
     
    # Some sanity checks
    # In principle subhalos should aready have the correct baryon mass: this is a check
    for sub in subs_this_snap:
        if sub.mass_baryon > parameters.mass_minimum_internal:
            if np.abs(sub.mass_baryon/sub.sum_mass_baryon(gals_this_snap)-1.)>1e-4:
                print('gals[''mass_baryon''] =',gals_this_snap[sub.gal_start_sid:sub.gal_end_sid]['mass_baryon'])
                print('sub.mass_gas_hot =',sub.mass_gas_hot)
                print('sub.mass_stars =',sub.mass_stars)
                halo=halos_this_snap[sub.halo_sid]
                print('halo.mass_gas_eject =',halo.mass_gas_eject)
                print('halo.mass_gas_hot =',halo.mass_gas_hot)
                print('halo.mass_stars =',halo.mass_stars)
                raise AssertionError('subhalo baryon mass discrepancy',
                    sub.graph_ID,sub.snap_ID,sub.sub_sid,sub.mass_baryon,sub.sum_mass_baryon(gals_this_snap))
    # Check that all halos are assigned to a halo
    for gal in gals_this_snap:
        if gal['b_exists'] and gal['halo_sid']==parameters.NO_DATA_INT:
            print(gal.dtype)
            print(gal)
            raise AssertionError('Galaxy not assigned to a halo')
    # Check that all halos are assigned to either a subhalo or a halo
    for gal in gals_this_snap:
        if gal['b_exists'] and gal['sub_sid']==parameters.NO_DATA_INT and gal['halo_sid']==parameters.NO_DATA_INT:
            print(gal.dtype)
            print(gal)
            raise AssertionError('Galaxy assigned neither subhalo nor halo')
    # Check that all galaxies have a finite virial speed
    for gal in gals_this_snap:
        if gal['b_exists'] and gal['v_vir']<1e-10:
            print(gal.dtype)
            print(gal)
            raise AssertionError('Galaxy virial speed is too low')
    # Check that all galaxies with cold gas have a disc scale length set
    for gal in gals_this_snap:
        if gal['mass_gas_cold']>parameters.mass_minimum_internal and \
           gal['radius_gas_cold']<parameters.length_minimum_internal:
            print(gal.dtype)
            print(gal)
            raise AssertionError('Galaxy disc scale length not set')
    return gals_this_snap
