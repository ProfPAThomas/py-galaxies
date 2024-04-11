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

import astropy.constants as c
import astropy.units as u
#import gc
import h5py
import numpy as np
np.seterr(all='raise')
import pickle
import subprocess
import sys

# Switch on traceback for warnings
import traceback
import warnings
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
warnings.showwarning = warn_with_traceback

# Add some functionality for ipython
def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')
if is_interactive():
    h5py.enable_ipython_completer()
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')


# In[2]:


# Ideally the following parameters should not be changed, but retained here as parameters for flexibility

# Location of code
C_DIR='code-C'
PYTHON_DIR='code-python'

# Default location of runtime parameters - can be over-ridden by command line parameter
FILE_PARAMETERS='input/input.yml'


# ### Imports of parameter and data classes, and of functions

# In[3]:


# Imports of py-galaxies python routines
sys.path.insert(1,PYTHON_DIR)

# Commons is a module used to store run-time parameters that are used during other module imports.
# It is a bit OTT as it contains very few entries and it is deleted after initialisation.
import commons

# The parameter class, used to store run-time parameters.
from parameters import C_parameters

# Read in parameters from yaml input files
# If running interactively this is input/input.par; otherwise add option to read from command line.
if is_interactive() or len(sys.argv)==1:
    parameters=C_parameters(FILE_PARAMETERS)
elif len(sys.argv)==2:
    parameters=C_parameters(sys.argv[1])
else:
    raise ValueError('Usage: python3 L-Galaxies.py <input.yml>')
# Define how verbose we want the terminal output to be
verbosity=parameters.verbosity
#commons.add('verbosity',verbosity)
# Define which graphs do we want to process
n_graph_start=parameters.n_graph_start
n_graph_end=parameters.n_graph_end
if n_graph_end==-1: n_graph_end=np.inf

# Open graph input file: needs to come before F_set_dt and F_update_parameters
graph_file=h5py.File(parameters.graph_file,'r')
# Update parameters with attributes from graph_file
parameters.F_update_parameters(graph_file)
if n_graph_end == np.inf:
    n_graph_end=parameters.n_graph
elif n_graph_end>parameters.n_graph: 
    raise ValueError('Not enough graphs in input file')
n_snap=parameters.n_snap

# Dictionary of run-time variables that it would be too messy to pass individually.
# It is permitted to add items to this dictionary in response to run-time flags.
# It will be converted to a ctypes structure at the end of the initialisation phase.
variables_dict=dict({
    'dt_snap':0.,      # Time between snaps
    'dt_halo':0.,      # Halo timestep
    'dt_gal':0.,       # Galaxy timestep
    'n_dt_halo':-1,    # Number of halo timesteps per snapshot
    'n_dt_gal':-1      # Number of galaxy timesteps per halo timestep
})

# Need to set the timesteps now, because that information is needed to determine the structure of
# the halo & subhalo classes and especially the galaxy arrays.
# This loads in the snapshot times and determines the number of timesteps required.
from misc import F_misc_set_dt
F_misc_set_dt(parameters)
if verbosity >=2:
    for i_snap in range(n_snap):
        print('i_snap, n_dt_halo, dt_halo, n_dt_gal, dt_gal')
        print(i_snap, parameters.n_dt_halo[i_snap], parameters.dt_halo[i_snap], parameters.n_dt_gal[i_snap], parameters.dt_gal[i_snap])
# This generates a class instance that holds the structure of the SFH arrays at all different timesteps
b_SFH=parameters.b_SFH
if b_SFH:
    from sfh import C_sfh
    sfh=C_sfh(parameters)
    sfh.F_create_header_file(parameters)
    
# These parameters are needed at import, so save to commons here
commons.update('b_SFH',parameters.b_SFH)
if b_SFH: 
    commons.update('sfh_n_bin',sfh.n_bin)
    # These counters need to be added to variables before saving as a C struct
    variables_dict['i_dt_sfh']=-1
    variables_dict['i_bin_sfh']=-1
    
# The conditional decorator and profilers
from profiling import conditional_decorator

# The graph class, used to store graphs for processing
from graphs import C_graph

# Halos
from halos import F_halos_create_header_file
F_halos_create_header_file()
from halos import F_halos_initialise
from halos import F_halos_template
halo_template=F_halos_template(parameters)
from halos import C_halo_output
halo_output=C_halo_output(parameters)

# Subhalos
from subs import F_subs_create_header_file
F_subs_create_header_file()
from subs import F_subs_initialise
from subs import F_subs_template
sub_template=F_subs_template(parameters)
from subs import C_sub_output
sub_output=C_sub_output(parameters)

# Galaxies
from gals import F_gals_template
gal_template=F_gals_template(parameters)
parameters.gal_template=gal_template
from gals import F_gals_create_header_file
F_gals_create_header_file()
from gals import C_gal_output
gal_output=C_gal_output(parameters)

# Now create the C parameters.h file.
# *** DO NOT MODIFY PARAMETERS AFTER THIS POINT ***
parameters.F_create_header_file()

# Create the variables structure to be used for passing variables to C
# *** DO NOT ADD VARIABLES AFTER THIS POINT ***
from misc import F_misc_create_variables_structure_and_header_file
C_variables=F_misc_create_variables_structure_and_header_file(variables_dict)
variables=C_variables(variables_dict)
del variables_dict

# Create a header file containing the list of all other header files
from misc import F_misc_create_all_headers_header_file
F_misc_create_all_headers_header_file()

# Import driver routine
from push_snap import F_push_snap     # Propagates info from last snapshot to current one


# ### Initialisation of SAM

# In[4]:


# Start CPU profiling, if required.
b_profile_cpu=parameters.b_profile_cpu
commons.update('b_profile_cpu',b_profile_cpu)
if b_profile_cpu:
    from codetiming import Timer
    from profiling import C_timer
    timer = C_timer()
    timer.start('Initialisation')

# Create counter to locate graphs within the galaxy output file
n_gal_graph_start=np.full(n_graph_end,parameters.NO_DATA_INT,dtype=np.int32)
n_gal=0  # This is the total number of galaxies over all graphs

# Create cooling table
import cooling
cooling_table = cooling.C_cooling(parameters)
cooling_table.F_create_header_file()

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


# ### C integration
# 
# Needs to be done after creation of all the header files, including for example that which stores the cooling tables.

# In[5]:


# C integration
import ctypes

# Make and load C-library
from misc import F_misc_create_Makefile
F_misc_create_Makefile(parameters)
subprocess.run(['make'],cwd=C_DIR)
L_C=ctypes.CDLL(C_DIR+'/lib_C.so')

# Define interfaces to C functions: currently only one of these
L_C.F_process_snap.argtypes=(np.ctypeslib.ndpointer(halo_template.dtype),ctypes.c_int,
                             np.ctypeslib.ndpointer(sub_template.dtype),ctypes.c_int,
                             np.ctypeslib.ndpointer(gal_template.dtype),ctypes.c_int,
                             C_variables)
L_C.F_process_snap.restype=None


# ### Loop over graphs, snapshots, halos, implementing the SAM
# 
# Note that loops over graphs can be done in parallel.
# 
# Also, F_push_snap needs to be serial, but all halos can be processed in parallel in F_process_snap.

# In[ ]:


# Loop over graphs
for i_graph in range(n_graph_start,n_graph_end):
    if verbosity >= 1: print('Processing graph',i_graph,flush=True)
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
        n_halo_this_snap=graph.snap_n_halo[i_snap]
        if n_halo_this_snap == 0:
            assert halos_last_snap == None
            continue
        if verbosity >= 2: print('Processing snapshot',i_snap,flush=True)
            
        # These timestepping parameters will be needed in the processing
        variables.dt_snap=parameters.dt_snap[i_snap]
        variables.dt_halo=parameters.dt_halo[i_snap]
        variables.dt_gal=parameters.dt_gal[i_snap]
        variables.n_dt_halo=parameters.n_dt_halo[i_snap]
        variables.n_dt_gal=parameters.n_dt_gal[i_snap]
            
        # This is the ministep, needed to track star formation histories
        if b_SFH:
            i_dt_sfh=sfh.i_dt_snap[i_snap]-1   # This gives the ministep, BEFORE updating
            variables.i_dt_sfh=i_dt_sfh
            variables.i_bin_sfh=sfh.i_bin[i_dt_sfh]
        
        # Initialise halo and subhalo properties.
        halos_this_snap = np.full(n_halo_this_snap,halo_template)
        F_halos_initialise(halos_this_snap,i_graph,i_snap,graph,parameters)
        n_sub_this_snap=graph.snap_n_sub[i_snap]
        if n_sub_this_snap>0:
            subs_this_snap = np.full(n_sub_this_snap,sub_template)
            F_subs_initialise(subs_this_snap,i_graph,i_snap,graph,parameters)        
        
        # Propagate information from progenitors to this generation.
        # Done as a push rather than a pull because sharing determined by progenitor.
        # Have to do this even if no progenitors in order to initialise galaxy array.
        # This is in python rather than C because of lots of graph lookups and variable-length arrays.
        gals_this_snap=F_push_snap(halos_last_snap,halos_this_snap,subs_last_snap,subs_this_snap,
                                   gals_last_snap,graph,parameters,variables)
        del halos_last_snap
        del subs_last_snap
        del gals_last_snap
        
        # Perform the astrophysics
        L_C.F_process_snap(halos_this_snap,n_halo_this_snap,subs_this_snap,n_sub_this_snap,
                           gals_this_snap,len(gals_this_snap),variables)

        # Once all halos have been done, output results
        # This could instead be done on a halo-by-halo basis in F_process_halos
        halo_output.append(halos_this_snap,parameters)
        if n_sub_this_snap >0: sub_output.append(subs_this_snap,parameters)
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
        if verbosity>1: print('graph processing time / s =',timer.timers[graph_str])
    if b_profile_mem:
        mem.stop(graph_str)
        if verbosity>1: print('graph memory usage / bytes =',mem.mem[graph_str])


# ###  Tidy up and exit

# In[ ]:


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

