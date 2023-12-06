#!/usr/bin/env python
# coding: utf-8

# # Python implementation of L-Galaxies
# 
# This is a playground to test out the possibility of using `python` as an interface into L-Galaxies.
# 
# Abbreviations used:
# * desc   – descendent
# * gal(s) – galax(y|ies)
# * prog   – progenitor
# * sub(s) – subhalo(s)
# 
# List index indentifiers:
# * _gid – relative to the graph
# * _sid - relative to the snap
# 
# – galaxies/orphans do not need an index identifier as they are numpy arrays defined per snap
# 
# Ptyhon type identifiers:
# * C_ - class
# * D_ - numpy dtype
# * F_ - function
# * b_ - boolean variable
# * c_ - constant (value may be set during parameter initialisation)
# * [ijk]_ - variable counter (integer)
# * n_ - total counter (integer)

# In[1]:


# Imports of generic python routines

#get_ipython().run_line_magic('load_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')

import astropy.constants as c
import astropy.units as u
#import gc
import h5py
#h5py.enable_ipython_completer()
import numpy as np
np.seterr(all='raise')
import pickle
import sys
# Switch on traceback for warnings
import traceback
import warnings

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


# In[2]:


# Location of code
C_DIR='code-C'
PYTHON_DIR='code-python'

# Development limiter
n_GRAPH=np.inf
#n_GRAPH=1000       # Change output files to 'test' to avoid over-writing!
n_GRAPH_START=0
#n_GRAPH_START=58

# Verbosity
VERBOSITY=0 # 1 - Major program steps only; 1/2 - Major/minor Counters; 3/4/5 - Debugging diags.

# List of runtime parameters
FILE_PARAMETERS='input/input.yml'


# ### Imports of parameter and data classes

# In[3]:


# Imports of py-galaxies python routines
sys.path.insert(1,PYTHON_DIR)

# Place to store common variables
# Use sparingly as hides dependencies in code, but sometimes necessary
import commons

# The parameter class, used to store run-time parameters
from parameters import C_parameters

# Read in parameters from yaml input files
parameters=C_parameters(FILE_PARAMETERS)
parameters.verbosity=VERBOSITY
commons.save('verbosity',parameters.verbosity)

# Open graph input file: needs to come before F_set_dt and F_update_parameters
graph_file=h5py.File(parameters.graph_file,'r')
# Update parameters with attributes from graph_file
parameters.F_update_parameters(graph_file)
n_graph=parameters.n_graph
n_snap=parameters.n_snap

# Need to set the timesteps now, because that information is needed to determine the structure of
# the halo & subhalo classes and especially the galaxy arrays.
# This loads in the snapshot times and determines the number of timesteps required.
from misc import F_set_dt
F_set_dt(parameters)
if VERBOSITY >=2:
    for i_snap in range(n_snap):
        print('i_snap, n_dt_halo, dt_halo, n_dt_gal, dt_gal')
        print(i_snap, parameters.n_dt_halo[i_snap], parameters.dt_halo[i_snap], parameters.n_dt_gal[i_snap], parameters.dt_gal[i_snap])
# This generates a class instance that holds the structure of the SFH arrays at all different timesteps
b_SFH=parameters.b_SFH
if b_SFH:
    from sfh import C_sfh
    sfh=C_sfh(parameters)
    parameters.sfh=sfh     # Useful to avoid having to optionally pass sfh as a function argument

# These parameters are needed at import, so save to commons here
commons.save('b_SFH',parameters.b_SFH)
if b_SFH: commons.save('sfh.n_bin',sfh.n_bin)

# The conditional decorator and profilers
from profiling import conditional_decorator

# The graph class, used to store graphs for processing
from graphs import C_graph

# The halo class, used to store halo properties
from halos import C_halo

# The halo_output class and methods used to output halos
from halos import C_halo_output

# The subhalo class, used to store subhalo properties
from subs import C_sub

# The subhalo_output class, used to output subhalos
from subs import C_sub_output

# The galaxy_output class, used to output galaxies
from gals import C_gal_output



# ### Initialisation of SAM

# In[4]:


# Start CPU profiling, if required.
b_profile_cpu=parameters.b_profile_cpu
commons.save('b_profile_cpu',b_profile_cpu)
if b_profile_cpu:
    from codetiming import Timer
    from profiling import C_timer
    timer = C_timer()
    timer.start('Initialisation')

# Create counter to locate graphs within the galaxy output file
n_graph_start=n_GRAPH_START
n_graph=min(n_graph-n_GRAPH_START,n_GRAPH)
n_gal_graph_start=np.full(n_graph_start+n_graph,parameters.NO_DATA_INT,dtype=np.int32)
n_gal=0

# Create cooling table
import cooling
cooling_table = cooling.C_cooling(parameters)
# Not sure if this is the best way to do it, but for now store all globals in parameters
parameters.cooling_table=cooling_table
    
# Create galaxy template
from gals import F_gal_template
gal_template=F_gal_template(parameters)
# Not sure if this is the best way to do it, but for now store all globals in parameters
parameters.gal_template=gal_template

# Create output buffers
halo_output=C_halo_output(parameters)
sub_output=C_sub_output(parameters)
gal_output=C_gal_output(parameters)

# Import driver routines
from driver import F_update_halos     # Propagates info from last snapshot to current one
from driver import F_process_halos    # Does all the astrophysics

if b_profile_cpu: timer.stop('Initialisation')
    
# Start Memory profiling, if required
b_profile_mem=parameters.b_profile_mem
if b_profile_mem: 
    import tracemalloc as tm
    tm.start()
    # These lines attempt to filter out profiling memory itself; not clear how well it does that.
    tm.Filter(False, tm.__file__)
    if b_profile_cpu: tm.Filter(False, codetiming.__file__)
    from profiling import C_mem
    mem = C_mem()


# ### Loop over graphs, snapshots, halos, implementing the SAM
# 
# Note that loops over graphs can be done in parallel.
# 
# Also, F_update_halos needs to be serial, but all halos can be processed in parallel in F_process_halos.

# In[5]:


# Loop over graphs
for i_graph in range(n_graph_start,n_graph_start+n_graph):
    if VERBOSITY >= 1: print('Processing graph',i_graph,flush=True)
    graph_str = 'Graph{:03d}'.format(i_graph)
    if b_profile_mem: mem.start(graph_str)
    if b_profile_cpu: timer.start(graph_str)
    graph = C_graph(i_graph,graph_file,parameters)
    
    # Keep track of location in galaxy output file
    n_gal_graph_start[i_graph]=n_gal
    
    # Loop over snapshots
    halos_last_snap = None
    subs_last_snap = None
    gals_last_snap = None
    for i_snap in range(parameters.n_snap):
        n_halo_in_snap=graph.snap_n_halo[i_snap]
        if n_halo_in_snap == 0:
            assert halos_last_snap == None
            continue
        if VERBOSITY >= 2: print('Processing snapshot',i_snap,flush=True)
            
        # These timestepping parameters will be needed in the processing, so save them to commons
        commons.save('dt_snap',parameters.dt_snap[i_snap])
        commons.save('dt_halo',parameters.dt_halo[i_snap])
        commons.save('dt_gal',parameters.dt_gal[i_snap])
        commons.save('n_dt_halo',parameters.n_dt_halo[i_snap])
        commons.save('n_dt_gal',parameters.n_dt_gal[i_snap])
            
        # This is the ministep, needed to track star formation histories
        if b_SFH:
            i_dt=sfh.i_dt_snap[i_snap]-1   # This gives the ministep, BEFORE updating
            commons.save('i_dt',i_dt)
            commons.save('i_bin_sfh',sfh.i_bin[i_dt])
        
        # Initialise halo and subhalo properties.
        # This returns a list of halo and subhalo instances
        # This may be slow: an alternative would be to use np arrays.
        halos_this_snap = [C_halo(i_graph,i_snap,i_halo,graph,parameters) for i_halo in 
                          graph.snap_first_halo_gid[i_snap]+range(graph.snap_n_halo[i_snap])]
        subs_this_snap = [C_sub(i_graph,i_snap,i_sub,graph,parameters) for i_sub in 
                         graph.snap_first_sub_gid[i_snap]+range(graph.snap_n_sub[i_snap])]
        
        # Propagate information from progenitors to this generation
        # Done as a push rather than a pull because sharing determined by progenitor
        # Have to do this even if no progenitors in order to initialise galaxy array
        gals_this_snap=F_update_halos(halos_last_snap,halos_this_snap,subs_last_snap,
                                          subs_this_snap,gals_last_snap,graph,parameters)
        del halos_last_snap
        del subs_last_snap
        del gals_last_snap
        #gc.collect() # garbage collection -- safe but very slow.
        
        # Process the halos
        for i_dt_halo in range(parameters.n_dt_halo[i_snap]):
            commons.save('i_dt_halo',i_dt_halo)
            F_process_halos(halos_this_snap,subs_this_snap,gals_this_snap,graph,parameters)
            
        # Once all halos have been done, output results
        # This could instead be done on a halo-by-halo basis in F_process_halos
        halo_output.append(halos_this_snap,parameters)
        if subs_this_snap != None: sub_output.append(subs_this_snap,parameters)
        if isinstance(gals_this_snap, np.ndarray):
            gal_output.append(gals_this_snap,parameters)
            n_gal+=len(gals_this_snap)
            
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
    if b_profile_cpu: 
        timer.stop(graph_str)
        if VERBOSITY>1: print('graph processing time / s =',timer.timers[graph_str])
    if b_profile_mem:
        mem.stop(graph_str)
        if VERBOSITY>1: print('graph memory usage / bytes =',mem.mem[graph_str])


# ###  Tidy up and exit

# In[6]:


if b_profile_mem: tm_snap = tm.take_snapshot()
if b_profile_cpu: timer.start('Finalise')
        
# Flush buffers, close files and exit
graph_file.close()
halo_output.close()
sub_output.close()
gal_output.close()

# Reopen galaxy output file to add graph start locations
# Don't simply do this before originally closing the output file, as the close flushes output buffers
gal_output=h5py.File(parameters.galaxy_file,'r+')
# Sanity check that have correct number of galaxies
assert n_gal==len(gal_output['Galaxies'])
# Add graph start locations as new dataset
dset=gal_output.create_dataset('Graph_start_locations',data=n_gal_graph_start,compression='gzip')
# And close
gal_output.close()

if b_profile_cpu: 
    timer.stop('Finalise')
    print(timer.timers['Finalise'])
    print(timer)
    timer.dump(parameters.profile_cpu_file)
    for key in Timer.timers:
        print('{:s}: {:g}'.format(key,Timer.timers[key]))
    with open(parameters.profile_cpu_file, 'ab') as f:
        pickle.dump(Timer.timers, f)
if b_profile_mem: 
    mem.dump(parameters.profile_mem_file)
    with open(parameters.profile_mem_file, 'ab') as f:
        pickle.dump(tm_snap, f)

