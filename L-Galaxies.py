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

# Imports of generic python routines

import astropy.constants as c
import astropy.units as u
from codetiming import Timer
import gc
import h5py
import numpy as np
import pickle
import sys
import tracemalloc as tm

# Switch on traceback for warnings
import traceback
import warnings

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback

# Location of code
C_DIR='code-C'
PYTHON_DIR='code-python'

# Development limiter
n_GRAPH=np.inf
#n_GRAPH=10       # Change output files to 'test' to avoid over-writing!

# Verbosity
VERBOSITY=1 # 1 - Major program steps only; 1/2 - Major/minor Counters; 3/4/5 - Debugging diags.

# List of runtime parameters
FILE_PARAMETERS='input/input.yml'


# ### Imports of parameter and data classes

# Imports of py-galaxies python routines
sys.path.insert(1,PYTHON_DIR)

# Place to store common variables
# Use sparingly as hides dependencies in code, but sometimes necessary
import commons

# The parameter class, used to store run-time parameters
from parameters import C_parameters

# The conditional decorator and profilers
# These may not be needed but we don't know until we read in the parameters
from profiling import conditional_decorator
from profiling import C_mem, C_timer

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

# Read in parameters from yaml input files
parameters=C_parameters(FILE_PARAMETERS,)
parameters.verbosity=VERBOSITY
commons.save('verbosity',parameters.verbosity)

b_profile_cpu=parameters.b_profile_cpu
commons.save('b_profile_cpu',b_profile_cpu)
if b_profile_cpu: 
    timer = C_timer()
    timer.start('Initialisation')

# Open graph input file: needs to come before F_update_parameters
graph_file=h5py.File(parameters.graph_file,'r')
# Update parameters with attributes from graph_file
parameters.F_update_parameters(graph_file)

# Create counter to locate graphs within the galaxy output file
n_graph=min(parameters.n_graph,n_GRAPH)
n_gal_graph_start=np.full(n_graph,parameters.NO_DATA_INT,dtype=np.int32)
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
    tm.start()
    # These lines attempt to filter out profiling memory itself; not clear how well it does that.
    tm.Filter(False, tm.__file__)
    if b_profile_cpu: tm.Filter(False, codetiming.__file__)
    mem = C_mem()

# ### Loop over graphs, snapshots, halos, implementing the SAM
# 
# Note that loops over graphs can be done in parallel.
# 
# Also, F_update_halos needs to be serial, but all halos can be processed in parallel in F_process_halos.

# Loop over graphs
for i_graph in range(n_graph):
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
    for i_snap in graph.snap_ID:
        if i_snap == parameters.NO_DATA_INT: 
            assert halos_last_snap == None
            continue
        if VERBOSITY >= 2: print('Processing snapshot',i_snap,flush=True)
            
        # Initialise halo and subhalo properties.
        # This returns a list of halo and subhalo instances
        # This may be slow: an alternative would be to use np arrays.
        halos_this_snap = [C_halo(i_graph,i_snap,i_halo,graph,parameters) for i_halo in 
                         graph.halo_start_gid[i_snap]+range(graph.n_halo_snap[i_snap])]
        subs_this_snap = None
        if graph.n_sub > 0:
            if graph.n_sub_snap[i_snap] > 0:
                subs_this_snap = [C_sub(i_graph,i_snap,i_sub,graph,parameters) 
                                     for i_sub in graph.sub_start_gid[i_snap]+range(graph.n_sub_snap[i_snap])]
        
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
        # The determination of timesteps could be done at initialisation (in update_parameters)
        dt_snap=((parameters.snap_table['time_in_years'][i_snap]- \
            parameters.snap_table['time_in_years'][i_snap-1]) * u.yr / parameters.units_time_internal).value
        parameters.dt_snap = dt_snap
        n_dt_halo=int(dt_snap*1.000001/parameters.timestep_halo_internal)+1
        parameters.n_dt_halo=n_dt_halo
        parameters.dt_halo=dt_snap/n_dt_halo
        if VERBOSITY >= 2: print('t_snap, n_dt_halo, dt_halo =',t_snap, n_dt_halo,parameters.dt_halo)
        for i_dt in range(n_dt_halo):
            F_process_halos(halos_this_snap,subs_this_snap,gals_this_snap,graph,parameters)
            
        # Once all halos have been done, output results
        # This could instead be done on a halo-by-halo basis in F_process_halos
        halo_output.append(halos_this_snap,parameters)
        if subs_this_snap != None: sub_output.append(subs_this_snap,parameters)
        if isinstance(gals_this_snap, np.ndarray):
            gals_this_snap 
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
    gc.collect()
    if b_profile_cpu: 
        timer.stop(graph_str)
        if VERBOSITY>1: print(timer.timers[graph_str])
    if b_profile_mem:
        mem.stop(graph_str)
        if VERBOSITY>1: print(mem.mem[graph_str])

# ###  Tidy up and exit

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
