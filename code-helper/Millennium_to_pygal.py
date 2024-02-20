#!/usr/bin/env python
# coding: utf-8

# # Script to convert Millenium tree files to py-gal input files
# 
# This is the guts of the io_tree.c read routine.
# Firstly `load_tree_table` reads the header info:
# ```
#   //read header on trees_** file
#   myfread(&Ntrees, 1, sizeof(int), tree_file);
#   myfread(&totNHalos, 1, sizeof(int), tree_file);
# 
#   TreeNHalos = mymalloc("TreeNHalos", sizeof(int) * Ntrees);
#   TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * Ntrees);
#   TreeNgals[0] = mymalloc("TreeNgals[n]", NOUT * sizeof(int) * Ntrees);
#   for(n = 1; n < NOUT; n++)
#     TreeNgals[n] = TreeNgals[n - 1] + Ntrees;
# 
#   myfread(TreeNHalos, Ntrees, sizeof(int), tree_file);
# 
#   if(Ntrees)
#     TreeFirstHalo[0] = 0;
#   /*Define a variable containing the number you have to jump to
#    * get from one firshalo to the next. */
#   for(i = 1; i < Ntrees; i++)
#     TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];
# ```
# In the above, `NOUT` is the number of snapshots
# 
# Then `load_tree` reads in the halo data:
# ```
#   Halo = mymalloc("Halo", sizeof(struct halo_data) * TreeNHalos[nr]);
#   myfseek(tree_file, sizeof(int) * (2 + Ntrees) + sizeof(struct halo_data) * TreeFirstHalo[nr], SEEK_SET);
#   myfread(Halo, TreeNHalos[nr], sizeof(struct halo_data), tree_file);
# ```
# where
# ```
# struct halo_data
# {
# 	/* merger tree pointers */
# 	int Descendant;
# 	int FirstProgenitor;
# 	int NextProgenitor;
# 	int FirstHaloInFOFgroup;
# 	int NextHaloInFOFgroup;
# 
#   /* properties of halo */
# 	int Len;
# 	float M_Mean200, M_Crit200, M_TopHat;
# 	float Pos[3];
# 	float Vel[3];
# 	float VelDisp;
# 	float Vmax;
# 	float Spin[3];
# 	long long MostBoundID;
# 
#   /* original position in subfind output */
# 	int SnapNum;
# 	int FileNr;
# 	int SubhaloIndex;
# 	float SubHalfMass;
# }
# ```

import h5py 
import numpy as np
from numpy.lib import recfunctions as rfn

infile='/Users/petert/data/MR/treedata/trees_063.5'
snapfile='/Users/petert/lgalaxies/Development_Branch/input/MRPlancksnaplist.txt'
outfile='/Users/petert/data/MR/treedata/pygal_063_5.hdf5'

# May want to restrict the number of trees to create small files for testing/development
n_tree_max=-1   # Use -1 to signify all trees
#n_tree_max=10  # or limit to a specific number

# These are parameters/attributes for the Hen15 version of the MR (W1_Planck)
Hubble_h=0.673
Omega_m=0.315
Omega_L=1-Omega_m
baryon_fraction=0.155
unit_length='Mpc/h'
unit_mass='1e10 Msun/h'
unit_speed='km/s'
unit_time='yr'   # For snapshot table; not otherwise needed

# We need to convert from mass to radius at some stage to estimate the sizes of halos,
# which is not given in the data.  We can use the formula R ~ (GM_crit200/100H^2)^(1/3)
# Then if R is in Mpc/h, M is in 1e10 Msun/h and G is in SI units, we have
# R = ((G kg s^2/m^3) (Msun/kg) / (100 (Mpc/m)) )^(1/3) (M/(Omega_L+Omega_m*(1+z)^3)^(1/3) 
# R = C_R (M/(Omega_L+Omega_m*(1+z)^3)^(1/3) where
C_R = (6.67e-11 * 2e30 / (100 * 3.086e22) )**(1./3.)

# Note that inputs get rescaled according to cosmology.
PartMass=0.0961104   # After scaling
ScalePos=0.960558
ScaleMass=1.11671
# Then presumably velocity dispersion should scale as sqrt(M/R)
ScaleVelDisp=np.sqrt(ScaleMass/ScalePos)

halo_data_dtype=([
    ('Descendant',np.int32),
    ('FirstProgenitor',np.int32),
    ('NextProgenitor',np.int32),
    ('FirstHaloInFOFgroup',np.int32),
    ('NextHaloInFOFgroup',np.int32),
    ('Len',np.int32),
    ('M_Mean200',np.float32),
    ('M_Crit200',np.float32),
    ('M_TopHat',np.float32),
    ('Pos',(np.float32,3)),
    ('Vel',(np.float32,3)),
    ('VelDisp',np.float32),
    ('Vmax',np.float32),
    ('Spin',(np.float32,3)),    
    ('MostBoundID',np.int64),
    ('SnapNum',np.int32),
    ('FileNr',np.int32),
    ('SubhaloIndex',np.int32),
    ('SubHalfMass',np.float32)
])

# Read in tree data
f=open(infile,'rb')
n_tree=np.fromfile(f,dtype=np.int32,count=1)[0]
n_halo=np.fromfile(f,dtype=np.int32,count=1)[0]
n_halo_in_tree=np.fromfile(f,dtype=np.int32,count=n_tree)
# first_halo_in_tree=np.cumsum(N_halos_in_tree) # Actually last halo
# first_halo_in_tree[1:]=first_halo_in_tree[:-1]
# first_halo_in_tree[0]=0
halos=np.fromfile(f,dtype=halo_data_dtype,count=n_halo)
# Check that at end of file: should get empty entry
dummy=np.fromfile(f,dtype=np.int32,count=1)
assert len(dummy)==0
f.close()

# Determine the location of the halos in each tree
print(n_tree)
print(n_halo)
print(n_halo_in_tree)
i_first_halo_in_tree=np.zeros(n_tree,int)
i_first_halo_in_tree[1:]=np.cumsum(n_halo_in_tree)[:-1]
n_graph=n_tree # Synonyms; use either depending upon context.

# For testing, restrict number of trees
if n_tree_max>-1: n_tree=min(n_tree_max,n_tree)
n_graph=n_tree
n_halo_in_tree=n_halo_in_tree[:n_tree]
print(n_halo_in_tree[58])
n_halo=np.sum(n_halo_in_tree)
halos=halos[:n_halo]

# We will need to know the location of each halo within its tree.
# This may seem trivial, but we are going to loop in a different order
halos=rfn.append_fields(halos,'loc',np.full(len(halos),-1),dtypes=np.int32,usemask=False)

# Add column to store location of main halo (not always the same as FirstHaloInFOFgroup)
halos=rfn.append_fields(halos,'main_halo',np.full(len(halos),-1),dtypes=np.int32,usemask=False)

# We will also need information from the snap table
snap_table=np.loadtxt(snapfile,usecols=[0,2,4],skiprows=1,
    dtype=[('snap_ID',np.int32),('redshift',np.float32),('time_in_years',np.float32)])
n_snap=len(snap_table)

# Because the FirstHaloInFOFgroup is not always the main_halo (defined as the one with the greatest Len)
# then we need to preprocess to store the latter.  We don't need to do that here.
# However, this is a test here to see if trees are well-formed before we hack them.

# Loop over creating groups for each tree/graph
for i_graph in range(n_graph):
    i_tree=i_graph # graph and tree are synonymous here; I will use both to emhasise context
    halos_in_tree=halos[i_first_halo_in_tree[i_tree]:i_first_halo_in_tree[i_tree]+n_halo_in_tree[i_tree]]
    halos_in_tree['loc']=np.arange(len(halos_in_tree),dtype=int)
    for i_snap in range(n_snap):
        halos_in_snap=halos_in_tree[halos_in_tree['SnapNum']==i_snap]
        if len(halos_in_snap)>0: 
            # Each FOF group has a distinct first subfind halo:
            i_first_halos=set(halos_in_snap['FirstHaloInFOFgroup'])
            assert len(i_first_halos)>0
            for i_first_halo in i_first_halos:
                # Find all the members of the FOF group
                i_halo=i_first_halo
                i_halos=[i_halo]
                i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                while i_halo != -1:
                    assert halos_in_tree['FirstHaloInFOFgroup'][i_halo]==i_first_halo
                    i_halos.append(i_halo)   # Hopefully not too slow here (no massive FOF groups)
                    i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                halos_in_FOF_group=halos_in_tree[i_halos]   # These are the Subfind (sub)halos in this Mega halo
                # It seems that the most massive halo is NOT always the first one.
                # So we need here to identify that
                main_halo=halos_in_FOF_group[np.argmax(halos_in_FOF_group['Len'])]
                halos_in_tree['main_halo'][halos_in_FOF_group['loc']]=main_halo['loc']
assert (halos['FirstHaloInFOFgroup']!=-1).all()

# Test that halos skip at most 1 generation in their progenitors
for i_tree in range(n_tree):
    i_offset=i_first_halo_in_tree[i_tree]
    for i_halo in i_offset+range(n_halo_in_tree[i_tree]):
        i_prog=halos['FirstProgenitor'][i_halo]  # This is relative to the graph
        if i_prog!=-1:
            if halos['SnapNum'][i_prog+i_offset]==halos['SnapNum'][i_halo]-1 or \
               halos['SnapNum'][i_prog+i_offset]==halos['SnapNum'][i_halo]-2: # Skipped a snapshot
                continue
            else:
                print(i_tree,halos[i_halo][['loc','FirstProgenitor','SnapNum']])

# Test that Descendant of Progenitor is the same halo, and fix if necessary
for i_tree in range(n_tree):
    i_offset=i_first_halo_in_tree[i_tree]
    for i_halo in i_offset+range(n_halo_in_tree[i_tree]):
        i_prog=halos['FirstProgenitor'][i_halo]  # This is relative to tree
        while i_prog!=-1:
            i_prog+=i_offset # This is now relative to entire halos array
            if halos['Descendant'][i_prog]!=halos['loc'][i_halo]:
                print(i_tree,halos[i_proj][['loc','Descendant','SnapNum']],halos[i_halo][['loc','SnapNum']])
                halos[i_prog]['Descendant']=halos[i_halo]['loc']
            i_prog=halos[i_prog]['NextProgenitor'] # Again relative to tree

# Test that halos skip at most 1 generation in their descendants
for i_tree in range(n_tree):
    i_offset=i_first_halo_in_tree[i_tree]
    for i_halo in i_offset+range(n_halo_in_tree[i_tree]):
        i_desc=halos['Descendant'][i_halo]  # This is relative to the graph
        if i_desc!=-1:
            if halos['SnapNum'][i_desc+i_offset]==halos['SnapNum'][i_halo]+1 or \
               halos['SnapNum'][i_desc+i_offset]==halos['SnapNum'][i_halo]+2: # Skipped a snapshot
                continue
            else:
                print(i_tree,halos[i_halo][['loc','Descendant','SnapNum']])

# Next we have to find and correct nonsense NextHaloInFOFgroup entries
for i_tree in range(n_tree):
    count1=0
    count2=0
    i_offset=i_first_halo_in_tree[i_tree]
    for i_halo in i_offset+range(n_halo_in_tree[i_tree]):
        i_next=halos['NextHaloInFOFgroup'][i_halo]  # This is relative to the graph
        while i_next!=-1:
            i_next+=i_offset # This is now relative to entire halos array
            if i_next>=n_halo:
                count1+=1
                halos[i_halo]['NextHaloInFOFgroup']=-1
                break
            if halos['SnapNum'][i_next]!=halos['SnapNum'][i_halo]:
                count2+=1
                halos[i_halo]['NextHaloInFOFgroup']=-1
                break
            i_next=halos[i_next]['NextHaloInFOFgroup'] # Again relative to tree
    if count1+count2>0: print(i_tree,count1,count2)

# # Descendants can skip a generation
# # This cell is just a test to show that happens.
# count=0
# count_100=0
# count_first=0
# len_first_max=0
# for i_halo in range(n_halo):
#     i_desc=halos['Descendant'][i_halo]  # This is relative to the graph
#     i_offset=i_first_halo_in_tree[np.where(i_first_halo_in_tree<=i_halo)[0][-1]]
#     if i_desc!=-1:
#         if halos['SnapNum'][i_desc+i_offset]!=halos['SnapNum'][i_halo]+1:
#             count+=1
#             #print('i_halo_snap, i_desc_snap, Len_halo, Len_desc =',halos['SnapNum'][i_halo],halos['SnapNum'][i_desc+i_offset],halos['Len'][i_halo],halos['Len'][i_desc+i_offset])
#             #if count ==1000: break
#             i_parent=halos['FirstHaloInFOFgroup'][i_halo]  # This is relative to the graph
#             i_parent_desc=halos['Descendant'][i_parent+i_offset]  # This is relative to the graph
#             while halos['SnapNum'][i_parent_desc+i_offset]!=halos['SnapNum'][i_parent+i_offset]+1:
#                 i_parent_next=halos['NextHaloInFOFgroup'][i_parent+i_offset]
#                 if i_parent_next==-1:
#                     count_first+=1
#                     if halos['Len'][i_halo]>=100: count_100+=1
#                     len_first_max=max(len_first_max,halos['Len'][i_parent+i_offset])
#                     break
#                 else:
#                     i_parent=i_parent_next
#                     i_parent_desc=halos['Descendant'][i_parent+i_offset]
# print('n_halo, count, fraction =',n_halo,count,count/n_halo)
# print('n_halo, count_first, fraction =',n_halo,count_first,count_first/n_halo)
# print('n_halo, count_100, fraction =',n_halo,count_100,count_100/n_halo)
# print('len_first_max =',len_first_max)

# So we need a way to fix this.
# 
# The only way that I can think of that has any merit is to copy the L-Galaxies data, graph by graph, adding in intermediate halos along the way (the (sub)halos in each tree need to be consecutive in the data set).
# 
# We will need to place the new halos in FOF groups: either the descendant of that of their progenitor, if it exists, or in one of their own.  We will take the other properties of the intermediate haloes to be the mean of the progenitor and descendant.
# 
# The addition (and elimination) of halos will require a relabelling of pointers.  This is achieved by creating a mapping from old pointers to new ones.  These pointers are relative to the tree.

# Allocate more than enough space for new halos (will truncate at end)
halos_new=np.empty(int(1.5*n_halo),dtype=halos.dtype)
n_halo_in_tree_new=np.zeros(n_tree,dtype=np.int32)
pointer_offset=2*n_halo # Used to distinguish pointers relative to new tree from those relative to old tree

i_halo_new=0
for i_tree in range(n_tree):
    # Mappings from location in old tree to new one
    # Set to be -1 just in case reference is made to a halo that is eliminated.  I don't think it should be!
    pointer_new=np.full(n_halo_in_tree[i_tree]+1,-2,dtype=np.int32)
    pointer_new[-1]=-1 # Because an index of -1 is the last entry in the array!
    tree_offset=i_halo_new
    i_halo_in_tree=0 # Will increment as add new halos
    # The following line appears to give a view even though ipython returns `False` on an `is` assertion.
    halos_in_tree=halos[i_first_halo_in_tree[i_tree]:i_first_halo_in_tree[i_tree]+n_halo_in_tree[i_tree]]
    halos_in_tree['loc']=np.arange(n_halo_in_tree[i_tree],dtype=np.int32)
    # Run through all the (sub)halos, checking whether they skip a generation.
    # Also, eliminate all halos that have no descendants (except those in the final snapshot);
    # We can only do that if we process in **reverse** snapshot order
    for i_snap in range(n_snap-1,-1,-1):
        # This time the following statement appears to give a copy rather than a view, which is not
        # surprising given the complex nature of the indexing.
        halos_in_snap=halos_in_tree[halos_in_tree['SnapNum']==i_snap]
        n_halo_in_snap=len(halos_in_snap)
        if n_halo_in_snap>0:
            if i_snap==n_snap-1: # last snapshot; accept all halos
                halos_new[i_halo_new:i_halo_new+n_halo_in_snap]=halos_in_snap
                pointer_new[halos_in_snap['loc']]=np.arange(i_halo_in_tree,i_halo_in_tree+n_halo_in_snap,dtype=int)
                i_halo_new+=n_halo_in_snap
                i_halo_in_tree+=n_halo_in_snap
            else:
                for halo in halos_in_snap:
                    if halo['Descendant']==-1: # Halo has no descendant; eliminate it from consideration
                        if halo['Len']>=100:
                            print('Warning: halo of particle number {:d} has no descendant'.format(halo['Len']))
                    else:
                        halo_desc=halos_in_tree[halo['Descendant']] # Because pointers are relative to the tree
                        if pointer_new[halo_desc['loc']]==-2: # Descendant has been eliminated from new tree
                            if halo['Len']>=100:
                                print('Warning: halo of particle number {:d} has eliminated descendant'.format(halo['Len']))
                        elif halo_desc['SnapNum']==halo['SnapNum']+1: # All is as it should be
                            halos_new[i_halo_new]=halo
                            pointer_new[halo['loc']]=i_halo_in_tree
                            i_halo_new+=1
                            i_halo_in_tree+=1
                        else: # We have skipped a snapshot; introduce an intermediate halo
                            # First add in existing halo
                            halos_new[i_halo_new]=halo
                            pointer_new[halo['loc']]=i_halo_in_tree 
                            descendant=halo['Descendant'] # Will need to remember old descendant location in tree
                            # Recast descendant link.
                            # This is relative to the new tree, so give offset to distinguish that: correct later.
                            halos_new[i_halo_new]['Descendant']=pointer_offset+i_halo_in_tree+1
                            i_halo_new+=1
                            i_halo_in_tree+=1
                            # We need now to determine the location in the halo list of the descendant;
                            # Remember that 'descendant' is location in the original tree
                            descendant_halo=halos_in_tree[descendant] # This is the existing halo, not the new copy.
                            # Now reset the progenitor of the original descendant
                            # *** Not yet done as not needed for py-gal ***
                            # The complication is to know where the descendant is in the new listing
                            # Next create the new one.  Note that it does not matter that this may break the snapshot ordering: that will be corrected when we generate the pygal dataset
                            halos_new[i_halo_new]['Descendant']=descendant
                            halos_new[i_halo_new]['FirstProgenitor']=-1 # Not following
                            halos_new[i_halo_new]['NextProgenitor']=-1 # Not following
                            halos_new[i_halo_new]['FirstHaloInFOFgroup']=pointer_offset+i_halo_in_tree # Itself
                            halos_new[i_halo_new]['NextHaloInFOFgroup']=-1 # Place in its own FOFgroup
                            halos_new[i_halo_new]['Len']=(halo['Len']+descendant_halo['Len'])//2
                            halos_new[i_halo_new]['M_Mean200']=(halo['M_Mean200']+descendant_halo['M_Mean200'])/2.
                            halos_new[i_halo_new]['M_Crit200']=(halo['M_Crit200']+descendant_halo['M_Crit200'])/2.
                            halos_new[i_halo_new]['M_TopHat']=(halo['M_TopHat']+descendant_halo['M_TopHat'])/2.
                            halos_new[i_halo_new]['Pos'][:]=(halo['Pos'][:]+descendant_halo['Pos'][:])/2.
                            halos_new[i_halo_new]['Vel'][:]=(halo['Vel'][:]+descendant_halo['Vel'][:])/2.
                            # Velocity dispersion and velocity should go as sqrt(M/r) or M^{1/3}
                            halos_new[i_halo_new]['VelDisp']=((halo['VelDisp']**3+descendant_halo['VelDisp']**3)/2.)**(1./3.)
                            halos_new[i_halo_new]['Vmax']=((halo['Vmax']**3+descendant_halo['Vmax']**3)/2.)**(1./3.)
                            # Take a mass-weighted average spin
                            halos_new[i_halo_new]['Spin'][:]=(halo['Spin'][:]*halo['Len']+descendant_halo['Spin'][:]*descendant_halo['Len'])/ \
                                                             (halo['Len']+descendant_halo['Len'])
                            halos_new[i_halo_new]['MostBoundID']=-1
                            halos_new[i_halo_new]['SnapNum']=i_snap+1
                            halos_new[i_halo_new]['FileNr']=halo['FileNr']
                            halos_new[i_halo_new]['SubhaloIndex']=-1 # Pointer to location in subfind output
                            halos_new[i_halo_new]['SubHalfMass']=(halo['SubHalfMass']+descendant_halo['SubHalfMass'])/2.
                            i_halo_new+=1
                            i_halo_in_tree+=1
            # At this point in the new tree, some halos will have been omitted.
            # That means that some of the FOF pointers may be pointing to non-existant halos.
            # So we need to look through all the (old) FOF groups and remake them eliminating the missing halos.
            # Each FOF group has a distinct first subfind halo:
            i_first_halos=set(halos_in_snap['FirstHaloInFOFgroup'])
            assert len(i_first_halos)>0
            for i_first_halo in i_first_halos:
                # Find the first true halo in the FOF group
                i_first_halo_corrected=i_first_halo
                while pointer_new[i_first_halo_corrected]==-2:
                    i_first_halo_corrected=halos_in_tree[i_first_halo_corrected]['NextHaloInFOFgroup']
                pointer_new[i_first_halo]=pointer_new[i_first_halo_corrected] # Have corrected first halo
                # Now run through rest of tree updating NextHaloInFOFgroup pointers to exclude eliminated halos
                i_next_halo=halos_in_tree[i_first_halo]['NextHaloInFOFgroup']
                while i_next_halo != -1:
                    i_next_halo_corrected=i_next_halo
                    while pointer_new[i_next_halo_corrected]==-2:
                        i_next_halo_corrected=halos_in_tree[i_next_halo_corrected]['NextHaloInFOFgroup']
                    pointer_new[i_next_halo]=pointer_new[i_next_halo_corrected]
                    i_next_halo=halos_in_tree[i_next_halo_corrected]['NextHaloInFOFgroup']
    n_halo_in_tree_new[i_tree]=i_halo_in_tree
    # Now need to run through and correct all the pointers
    for halo in halos_new[tree_offset:i_halo_new]:
        pointer_old=halo['Descendant']
        if pointer_old==-1:
            halo['Descendant']=-1
        elif pointer_old<pointer_offset:
            assert pointer_new[pointer_old]!=-2
            halo['Descendant']=pointer_new[pointer_old]
        else:
            halo['Descendant']=pointer_old-pointer_offset
        halo['FirstProgenitor']=-1
        halo['NextProgenitor']=-1
        pointer_old=halo['FirstHaloInFOFgroup']
        if pointer_old==-1:
            raise valueError("halo['FirstHaloInFOFgroup']=-1")
        elif pointer_old<pointer_offset:
            assert pointer_new[pointer_old]!=-2
            halo['FirstHaloInFOFgroup']=pointer_new[pointer_old]
        else:
            halo['FirstHaloInFOFgroup']=pointer_old-pointer_offset
        pointer_old=halo['NextHaloInFOFgroup']
        if pointer_old==-1:
            halo['NextHaloInFOFgroup']=-1
        elif pointer_old<pointer_offset:
            assert pointer_new[pointer_old]!=-2
            halo['NextHaloInFOFgroup']=pointer_new[pointer_old]
        else:
            halo['NextHaloInFOFgroup']=pointer_old-pointer_offset
n_halo_new=i_halo_new
print('n_halo, n_halo_new =',n_halo,n_halo_new)
halos_new=halos_new[:n_halo_new]
assert np.sum(n_halo_in_tree_new)==n_halo_new
assert (halos_new['FirstHaloInFOFgroup']!=-1).all()

# Done in separate cell so don't always have to execute
# halos_old=halos.copy()
halos=halos_new
n_halo=len(halos)
n_halo_in_tree=n_halo_in_tree_new
i_first_halo_in_tree=np.zeros(n_tree,int)
i_first_halo_in_tree[1:]=np.cumsum(n_halo_in_tree)[:-1]

# Now we have some work to do to turn this into halos and subhalos.
# We will take the halo to be (a copy of) the main (most massive) halo in each FOF group; an alternative would be to use the sum of the halos withing the FOF group, but that leads to major issues with determining properties.  This also matches the current use within L-Galaxies.
# All L-Galaxies halos are subhalos within py-gal.
# 
# The `FirstHaloInFOFgroup` array labels halos from 0 within each tree.
# We will loop over snapshots (from early times to the present).
# Within each snapshot, we will identify halos (as opposed to subhalos) by the number of `FirstHaloInFOFgroup` entries.
# Note that pointers in the L-Galaxies data set seem to be relative to halos in the same tree, which is good as that is what we want for graphs also.

# # This cell is just to make a few plots to check that numbers look sensible

# plt.figure(figsize=[12,8])
# plt.loglog(halos['M_Mean200']**(1./3.),halos['VelDisp'],'.')
# plt.plot(halos_new['M_Mean200']**(1./3.),halos_new['VelDisp'],'+')
# plt.xlabel(r'$M_\mathrm{Mean,200}/10^{10}$M$_\odot$')
# plt.ylabel(r'$\sigma/$km$\,$s$^{-1}$')

# Now need to convert to py-gal format.
# Have yet to decide what that should be, so let's define it here!
# It is an HDF5 file so I will list Groups, Attributes and Datasets
# 
# Note: I will assume that the format requires every halo to have a subhalo;
# if that does not exist then simply copy the halo entry.  This is much easier
# to do here than in the guts of py-galaxies because of the need to renumber if
# subhalos are created on the fly.
# 
# Descendant pointers are tricky because we reorder the (sub)halos.  To keep track, we have to introduce an extra property to keep track of the halo's position in the L-Galaxies listing.  We first point to these old positions, then relabel the pointers at the end, once the new locations in the pygal listing are known.
# 
# Group: '/'
# * Attributes:
#   - n_graph : number of graphs in file
#   - n_snap  : number of shapshots in file
#   - Omega_m : density parameter
#   - Hubble_h : Hubble parameter at z=0 in units of 100 km/s/Mpc
#   - baryon_fraction : Omega_b/Omega_m
#   - unit_length : description of length unit (to be used to set value in input.yaml)
#   - unit_mass : description of mass unit (ditto)
#   - unit_speed : description of speed unit (ditto)
#   - unit_time : description of time unit (ditto)
#   - other_attributes : anything that may prove of use [optional]
# * Datasets:
#   - snap_table[n_snap,3] : IDs, redshift and times of available snapshots
# * Groups:
#   - graph_#, # in range(n_graph) : containing the graph data
#     
# Group: 'graph_#'
# * Attributes:
#   - first_snap : ID of first snap to be populated with halo(s)
#   - n_halo : number of halos in graph
#   - n_halo_desc : total number of descendant halos
#   - n_sub : number of subhalos in graph
#   - n_sub_desc : total number of descendant subhalos
# * Datasests:
#   - snap_n_halo[n_snap] : number of halos in each snapshot
#   - snap_first_halo[n_snap] : first halo in each snapshot (redundant)
#   - halo_n_sub[n_halo] : Number of subhalos of this halo
#   - halo_first_sub[n_halo] : location in subhalo arrays of first subhalo
#   - halo_n_desc[n_halo] : Number of descendants
#   - halo_first_desc[n_halo] : location in desc arrays of first descendant
#   - halo_desc_contribution[n_halo_desc] : contribution of halo to descendant
#   - halo_desc_halo[n_halo_desc] : halo corresponding to each descendant
#   - halo_n_prog[n_halo] : Number of progenitors [optional]
#   - halo_first_prog[n_halo] : location in graph arrays of first progenitor [optional]
#   - halo_mass[n_halo] : only one measure allowed!
#   - halo_pos[n_halo] : comoving mean position (weighted by mass)
#   - halo_vel[n_halo] : mean velocity (weighted by mass)
#   - halo_rms_radius[n_halo] : halo rms radius [optional]
#   - halo_rms_speed[n_halo] : halo rms speed
#   - halo_half_mass_radius[n_halo] : radius enclosing half the mass
#   - halo_v_max[n_halo] : maximum circular velocity
#   - halo_spin[n_halo,3] : specific angular momentum
#   - snap_n_sub[n_snap] : number of subhalos in each snapshot
#   - snap_first_sub[n_snap] : first subhalo in each snapshot (redundant)
#   - sub_n_desc[n_sub] : Number of descendant subhalos
#   - sub_first_desc[n_sub] : location in desc arrays of first descendant
#   - sub_desc_contribution[n_sub_desc] : contribution of sub to descendant
#   - sub_desc_sub[n_sub_desc] : sub corresponding to each descendant
#   - sub_mass[n_sub] : only one measure allowed!
#   - sub_pos[n_sub] : comoving mean position (weighted by mass)
#   - sub_vel[n_sub] : mean velocity (weighted by mass)
#   - sub_rms_radius[n_sub] : sub rms radius [optional]
#   - sub_rms_speed[n_sub] : sub rms speed
#   - sub_half_mass_radius[n_sub] : radius enclosing half the mass
#   - sub_v_max[n_sub] : maximum circular velocity
#   - sub_spin[n_sub,3] : specific angular momentum

# Because the FirstHaloInFOFgroup is not always the main_halo (defined as the one with the greatest Len)
# then we need to preprocess to store the latter
# Loop over creating groups for each tree/graph
for i_graph in range(n_graph):
    i_tree=i_graph # graph and tree are synonymous here; I will use both to emhasise context
    halos_in_tree=halos[i_first_halo_in_tree[i_tree]:i_first_halo_in_tree[i_tree]+n_halo_in_tree[i_tree]]
    halos_in_tree['loc']=np.arange(len(halos_in_tree),dtype=int)
    for i_snap in range(n_snap):
        halos_in_snap=halos_in_tree[halos_in_tree['SnapNum']==i_snap]
        if len(halos_in_snap)>0: 
            # Each FOF group has a distinct first subfind halo:
            i_first_halos=set(halos_in_snap['FirstHaloInFOFgroup'])
            assert len(i_first_halos)>0
            for i_first_halo in i_first_halos:
                # Find all the members of the FOF group
                i_halo=i_first_halo
                i_halos=[i_halo]
                i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                while i_halo != -1:
                    assert halos_in_tree['FirstHaloInFOFgroup'][i_halo]==i_first_halo
                    i_halos.append(i_halo)   # Hopefully not too slow here (no massive FOF groups)
                    i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                halos_in_FOF_group=halos_in_tree[i_halos]   # These are the Subfind (sub)halos in this Mega halo
                # It seems that the most massive halo is NOT always the first one.
                # So we need here to identify that
                main_halo=halos_in_FOF_group[np.argmax(halos_in_FOF_group['Len'])]
                halos_in_tree['main_halo'][halos_in_FOF_group['loc']]=main_halo['loc']
assert (halos['FirstHaloInFOFgroup']!=-1).all()

# Open HDF5 file for output
try:
    f.close()
except:
    pass
f=h5py.File(outfile,'w')

# Write header attributes: these are specific to the Millennium simulation
f.attrs['n_graph']=n_tree
f.attrs['n_snap']=n_snap
f.attrs['Hubble_h']=Hubble_h
f.attrs['Omega_m']=Omega_m
f.attrs['baryon_fraction']=baryon_fraction
f.attrs['unit_length']=unit_length
f.attrs['unit_mass']=unit_mass
f.attrs['unit_speed']=unit_speed
f.attrs['unit_time']=unit_time

# Create snap_table dataset
f.create_dataset('snap_table',data=snap_table)

# Loop over creating groups for each tree/graph
for i_graph in range(n_graph):
    i_tree=i_graph # graph and tree are synonymous here; I will use both to emhasise context
    g=f.create_group('graph_'+str(i_graph))
    i_halo_in_graph=0
    i_sub_in_graph=0
    i_halo_desc=0
    i_sub_desc=0
    halos_in_tree=halos[i_first_halo_in_tree[i_tree]:i_first_halo_in_tree[i_tree]+n_halo_in_tree[i_tree]]
    n_halo_in_graph=len(set(halos_in_tree['FirstHaloInFOFgroup']))
    n_sub_in_graph=len(halos_in_tree)
    g.attrs['n_halo']=n_halo_in_graph
    g.attrs['n_sub']=n_halo_in_tree[i_graph]
    snap_n_halo=np.full(n_snap,0,int)
    snap_n_sub=np.full(n_snap,0,int)
    snap_first_halo=np.full(n_snap,-1,int)
    snap_first_sub=np.full(n_snap,-1,int)
    halo_n_sub=np.full(n_halo_in_graph,0,int)
    halo_first_sub=np.full(n_halo_in_graph,-1,int)
    # In Subfind, there is a single descendant for each (sub)halo
    n_halo_desc=n_halo_in_graph  # Will be too long, truncate later
    halo_first_desc=np.full(n_halo_in_graph,-1,int)
    halo_n_desc=np.full(n_halo_in_graph,1,int)
    halo_desc_contribution=np.full(n_halo_desc,1.) 
    halo_desc_halo=np.full(n_halo_desc,-1,int)
    halo_loc=np.full(n_halo_in_graph,-1,int)  # to enable relabelling of pointers; location in new list
    halo_mass=np.full(n_halo_in_graph,0.)
    halo_pos=np.full([n_halo_in_graph,3],0.)
    halo_vel=np.full([n_halo_in_graph,3],0.)
    halo_rms_speed=np.full(n_halo_in_graph,0.)
    halo_rms_radius=np.full(n_halo_in_graph,0.)
    halo_half_mass_radius=np.full(n_halo_in_graph,0.)
    halo_v_max=np.full(n_halo_in_graph,0.)
    halo_spin=np.full([n_halo_in_graph,3],0.)
    # Now subhalos
    # In Subfind, there is a single descendant for each (sub)halo
    n_sub_desc=n_sub_in_graph  # Will be too long, truncate later
    sub_first_desc=np.full(n_sub_in_graph,-1,int)
    sub_n_desc=np.full(n_sub_in_graph,1,int)
    sub_desc_contribution=np.full(n_sub_desc,1.) 
    sub_desc_sub=np.full(n_sub_desc,-1)
    sub_loc=np.full(n_sub_in_graph,-1,int)  # To enable relabelling of pointers
    sub_mass=np.full(n_sub_in_graph,0.)
    sub_pos=np.full([n_sub_in_graph,3],0.)
    sub_vel=np.full([n_sub_in_graph,3],0.)
    sub_rms_speed=np.full(n_sub_in_graph,0.)
    sub_rms_radius=np.full(n_sub_in_graph,0.)
    sub_half_mass_radius=np.full(n_sub_in_graph,0.)
    sub_v_max=np.full(n_sub_in_graph,0.)
    sub_spin=np.full([n_sub_in_graph,3],0.)
    for i_snap in range(n_snap):
        C_R_snap=C_R*(Omega_L+Omega_m*(1+snap_table[i_snap]['redshift'])**3)**(-1./3.)
        halos_in_snap=halos_in_tree[halos_in_tree['SnapNum']==i_snap]
        n_sub_in_snap=len(halos_in_snap)
        if n_sub_in_snap==0: 
            first_snap=i_snap+1  # Do not accumulate, Just in case there is a gap in the snaplist
        else:
            snap_n_sub[i_snap]=n_sub_in_snap
            snap_first_halo[i_snap]=i_halo_in_graph
            snap_first_sub[i_snap]=i_sub_in_graph
            # Each FOF group has a distinct first subfind halo:
            i_first_halos=set(halos_in_snap['FirstHaloInFOFgroup'])
            n_halo_in_snap=len(i_first_halos)
            assert n_halo_in_snap>0
            snap_n_halo[i_snap]=n_halo_in_snap
            for i_first_halo in i_first_halos:
                # Find all the members of the FOF group
                i_halo=i_first_halo
                i_halos=[i_halo]
                i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                while i_halo != -1:
                    assert halos_in_tree['FirstHaloInFOFgroup'][i_halo]==i_first_halo
                    i_halos.append(i_halo)   # Hopefully not too slow here (no massive FOF groups)
                    i_halo=halos_in_tree['NextHaloInFOFgroup'][i_halo]
                n_sub=len(i_halos)
                assert n_sub>0
                halo_n_sub[i_halo_in_graph]=n_sub
                halo_first_sub[i_halo_in_graph]=i_sub_in_graph
                halos_in_FOF_group=halos_in_tree[i_halos]   # These are the Subfind (sub)halos in this Mega halo
                # Pointers to descendants: this is tricky as we have not generated them yet!
                # Instead we remember the location in the original halo list, then relabel later
                main_halo=halos_in_tree[halos_in_FOF_group['main_halo'][0]]
                i_desc=main_halo['Descendant']
                if i_desc == -1:
                    pass
#                     if i_snap != n_snap-1:
#                         print('***Warning: main halo has no descendant***')
#                         print('i_graph, i_snap, i_halo_in_graph, i_desc =',i_graph, i_snap, i_halo_in_graph, i_desc)
                else:
                    halo_first_desc[i_halo_in_graph]=i_halo_desc
                    i_desc_main=halos_in_tree[i_desc]['main_halo']
#                     if i_desc_main!=i_desc:
#                         print('***Warning: descendant of main halo may not be a main halo***')
#                         print('i_graph, i_halo_in_graph, i_halo_desc, i_desc, i_desc_first] =',\
#                                i_graph, i_halo_in_graph, i_halo_desc, i_desc, i_desc_first)
                    # As this is a descendant halo, not a subhalo, we are free to redefine
                    halo_desc_halo[i_halo_desc]=halos_in_tree[i_desc_main]['loc']  # To be updated later
                    i_halo_desc+=1
                # Properties
                # In L-Galaxies the main subhalo is declared to be the enclosing FOF halo
                halo_loc[i_halo_in_graph]=main_halo['loc']
                halo_mass[i_halo_in_graph]=main_halo['Len']*PartMass
                halo_pos[i_halo_in_graph]=main_halo['Pos']              # Should this be scaled?
                halo_vel[i_halo_in_graph]=main_halo['Vel']              # Should this be scaled?
                # We don't have a radius for halos in the L-Galaxies data, so determine it from the mass
                # and mean density (assumed to be 200 times the critical density) - see start of code
                radius=C_R_snap*halo_mass[i_halo_in_graph]**(1./3.)
                halo_rms_radius[i_halo_in_graph]=radius/np.sqrt(3.)
                # Assuming an isothermal sphere, then half mass radius is half the total radius
                halo_half_mass_radius[i_halo_in_graph]=0.5*radius
                halo_rms_speed[i_halo_in_graph]=np.sqrt(3.)*main_halo['VelDisp']
                halo_v_max[i_halo_in_graph]=main_halo['Vmax']
                halo_spin[i_halo_in_graph]=main_halo['Spin']
                i_halo_in_graph+=1
                # Loop over (sub) halos
                assert len(halos_in_FOF_group)==n_sub
                for sub in halos_in_FOF_group:
                    # Pointers
                    sub_loc[i_sub_in_graph]=sub['loc']
                    i_desc=sub['Descendant']
                    if i_desc != -1:
                        sub_first_desc[i_sub_in_graph]=i_sub_desc
                        sub_desc_sub[i_sub_desc]=halos_in_tree[i_desc]['loc']  # To be updated later
                        i_sub_desc+=1
                    # Properties
                    sub_mass[i_sub_in_graph]=sub['Len']*PartMass
                    sub_pos[i_sub_in_graph]=sub['Pos']
                    sub_vel[i_sub_in_graph]=sub['Vel']
                    # We don't have a radius for halos in the L-Galaxies data, so determine it from the mass
                    # and mean density (assumed to be 200 times the critical density) - see start of code
                    radius=C_R_snap*sub_mass[i_sub_in_graph]**(1./3.)
                    sub_rms_radius[i_sub_in_graph]=radius/np.sqrt(3.)
                    # Assuming an isothermal sphere, then half mass radius is half the total radius
                    sub_half_mass_radius[i_sub_in_graph]=0.5*radius
                    sub_rms_speed[i_sub_in_graph]=np.sqrt(3.)*sub['VelDisp']
                    sub_v_max[i_sub_in_graph]=sub['Vmax']
                    sub_spin[i_sub_in_graph]=sub['Spin']
                    i_sub_in_graph+=1
    # Sanity checks
    assert i_halo_in_graph==n_halo_in_graph
    assert i_sub_in_graph==n_sub_in_graph
    # Now have to relabel the descendant pointers
    assert i_halo_desc <= n_halo_desc
    n_halo_desc=i_halo_desc
    g.attrs['n_halo_desc']=n_halo_desc
    # We have a list of subfind locations for descendant halos; now we need to relabel to new locations
    for i_halo_desc in range(n_halo_desc):
        locs=np.where(halo_loc==halo_desc_halo[i_halo_desc])[0]
        assert len(locs)==1
        halo_desc_halo[i_halo_desc]=locs[0]
    assert i_sub_desc <= n_sub_desc
    n_sub_desc=i_sub_desc
    g.attrs['n_sub_desc']=n_sub_desc
    for i_sub_desc in range(n_sub_desc):
        locs=np.where(sub_loc==sub_desc_sub[i_sub_desc])[0]
        assert len(locs)==1
        sub_desc_sub[i_sub_desc]=locs[0]
    # Save results to HDF5 file
    g.attrs['first_snap']=first_snap
    g.create_dataset('snap_n_halo',data=snap_n_halo)
    g.create_dataset('snap_first_halo',data=snap_first_halo)
    g.create_dataset('snap_n_sub',data=snap_n_sub)
    g.create_dataset('snap_first_sub',data=snap_first_sub)
    g.create_dataset('halo_n_desc',data=halo_n_desc)
    g.create_dataset('halo_first_desc',data=halo_first_desc)
    g.create_dataset('halo_desc_contribution',data=halo_desc_contribution[:n_halo_desc])
    g.create_dataset('halo_desc_halo',data=halo_desc_halo[:n_halo_desc])
    g.create_dataset('halo_mass',data=halo_mass)
    g.create_dataset('halo_pos',data=halo_pos)
    g.create_dataset('halo_vel',data=halo_vel)
    g.create_dataset('halo_rms_radius',data=halo_rms_radius)
    g.create_dataset('halo_half_mass_radius',data=halo_half_mass_radius)
    g.create_dataset('halo_rms_speed',data=halo_rms_speed)
    g.create_dataset('halo_v_max',data=halo_v_max)
    g.create_dataset('halo_spin',data=halo_spin)
    g.create_dataset('halo_n_sub',data=halo_n_sub)
    g.create_dataset('halo_first_sub',data=halo_first_sub)
    g.create_dataset('sub_n_desc',data=sub_n_desc)
    g.create_dataset('sub_first_desc',data=sub_first_desc)
    g.create_dataset('sub_desc_contribution',data=sub_desc_contribution[:n_sub_desc])
    g.create_dataset('sub_desc_sub',data=sub_desc_sub[:n_sub_desc])
    g.create_dataset('sub_mass',data=sub_mass)
    g.create_dataset('sub_rms_radius',data=sub_rms_radius)
    g.create_dataset('sub_half_mass_radius',data=sub_half_mass_radius)
    g.create_dataset('sub_pos',data=sub_pos)
    g.create_dataset('sub_vel',data=sub_vel)
    g.create_dataset('sub_rms_speed',data=sub_rms_speed)
    g.create_dataset('sub_v_max',data=sub_v_max)
    g.create_dataset('sub_spin',data=sub_spin)
f.close()

assert (halos['FirstHaloInFOFgroup']!=-1).any()

