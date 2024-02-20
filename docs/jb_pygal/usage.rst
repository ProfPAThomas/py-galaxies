Usage
=====

Installation
------------

The best way to obtain a working version of py-galaxies is to download from the `github repository <https://github.com/ProfPAThomas/py-galaxies>`_.

You will need an up-to-date version of python3 with the following modules installed:

* astropy
* codetiming (if you wish to switch on the built-in CPU profiling)
* cProfile (if you wish to do CPU profiling on the python part of the code)
* ctypes
* h5py
* jupyter
* matplotlib (if you wish to use the supplied plotting routines)
* numpy
* pickle
* profile (if you wish to use the built-in memory profiling)
* pyYAML
* subprocess
* seaborn (if you wish to use the supplied plotting routines)
* sys
* tracemalloc (if you wish to use the built-in memory profiling)
* warnings

The C code has been compiled and tested only using gcc (Apple clang version 14.0.3) but will hopefully work with any up to date C compiler.

Directory structure
-------------------

Not every directory/file is listed below; these are the important ones.

* L-Galaxies.ipynb - jupyter notebook version of L-Galaxies.py, very useful for code development.
* L-Galaxies.py - the main calling routine.
* MergerGraphs/ - contains an example merger graph; place here links to other input merger graph files.
* README.md - github front page.
* anal/ - contains routines to analyse and plot the results of runs.
* aux_files/ - contains an example cooling function; place here links to alternative cooling functions.  Will also contain yield tables once that is implemented.
* code-C/ - contains the C functions that perform the astrophysics.
* code-helper/ - python routines to generate files of use when running py-galaxies.
* code-python - contains python modules that are used during initialisation, and routines to set up the data structures for halos, subhalos and galaxies.  Up to v0.2 also contains the main python bookkeeping routine driver.py.
* docs/ - contains the documentaion, mostly written in rst notation for use with ReadTheDocs.
* input - contains input YAML files that specify the parameters of the run.
* output - contains some (small) example output files; place here links to destination directories/folders for the output.
* scratch - test code used during development; can be ignored.

MergerGraphs
^^^^^^^^^^^^

* Mill5_100.hdf5 - an example merger graph input file containing the first 100 trees from Millennium file 5.

anal
^^^^

* Comparison.ipynb/comparison.yml - compares the output of two runs; useful only for code development.
* GalaxyPlots.ipynb/galaxy.yml - makes plots of galaxy-related quantities.
* HaloPlots.ipynb/halo.yml - makes plots of halo-related quantities.
* SFHPlots.ipynb/sfh.yml - plots of star formation histories.
* SubhaloPlots.ipynb/subhalo.yml - makes plots of subhalo-related quantities.

aux-files
^^^^^^^^^

* Lambda_table.npz - the default cooling tables.

code-C
^^^^^^

* bh_agn.c - Functions related to black holes and agn activity.
* cooling.c - Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy.
* halo_infall.c - Functions for accretion of primordial gas, and of reincorporation of ejected gas onto halos.
* mergers.c - Functions related to galaxy mergers.
* sfh.c - Function to merge SFH bins if required.
* star_formation_and_feedback.c - Functions to make stars from galaxies and to provide feedback of cold gas from SN.

In addition to the above, there may be a Makefile and various header (.h) files that are automatically generated when py-galaxies is run.

code-helper
^^^^^^^^^^^

* generate_cooling_table.py/Lambda_table.hdf5 - generates a version of the cooling table for use with py-galaxies from that used in the standard L-Galaxies distribution.
* Millennium_to_pygal.py - converts Millennium tree files to py-galaxies input format (generates 1 imperfect tree on file 5, so may contain a minor bug).

code-python
^^^^^^^^^^^

* commons.py - used to hold common (ie global) variables during initialisation phase of the run.
* cooling.py - during initialisation, converts the cooling table to a dimensionless format and writes out to code-c/cooling.h
* driver.py - contains the two main workhorse routines:
  - F_update_halos - propagates halos, subhalos and galaxies from one snapshot to the next.
  - F_process_halos - perform all the astrophysics on halos.
* gals.py - set up the galaxy data structures.
* graphs.py - set up the graph data structures, plus associated methods.
* halo.py - set up the halo data structures, plus associated methods.
* misc.py - miscellaneous helper routines, mostly associated with writing out the C Makefile and header files.
* parameters.py - set up the parameters class, plus methods to read in the parameters from input.yml.
* sfh.py - set up the star formation history binning.
* subs.py - set up the subhalo data structures, plus associated methods.

docs
^^^^

* Makefile - :code:`make html` will generate the documentation in ReadTheDocs format.
* _build/ - contains the documentation after building.
* docs/ - miscellaneous documents, mostly historical.
* index.rst - main rst page listing all the documentation to be generated.
* jb_pygal/ - contains all the documentation in rst format.
* requirements.txt - a list of python modules that are needed for documentation generation.

input
^^^^^

* input.yml - a pointer to the YAML file that contains all the runtime parameters.  This is the file that is read in py-galaxies.
* input_test.yml - an example input file that processes just 10 graphs.
* input_timing.yml - an example input file that processes 1000 graphs.  Unfortunately the input data file is too large for github and so will need to be regenerated using code-helper/Millennium_to_pygal.py and the original Millennium File 5 tree file.

output
^^^^^^

* v0.1_10_*.hdf5 - output files from a test run on version 0.1 of the code (python only).
* v0.2_10_*.hdf5 - output files from a test run on version 0.2 of the code (initial C integration).
  

Installation testing
--------------------

By default, the installation should be set up to run a small test using just 10 merger trees from Millennium file 5.  To check the installation, perform the following steps:

* :code:`cd <Installation directory>` - change directory to the top level py-galaxies directory.
* :code:`python3 L-Galaxies.py` - by default this should run py-galaxies using the input file input/input_test.yml.  It should create the following output files in output/:

  - test_10_SFH.hdf5 - containing the star formation history bin structure;
  - test_10_GalaxyOutput.hdf5 - galaxy properties;
  - test_10_SubhaloOutput.hdf5 - subhalo properties;
  - test_10_HaloOutput.hdf5 - halo properties.

* Start a jupyter notebook and navigate to anal/Comparison.ipynb.  Then run all cells.  That will compare the output in test_10_GalaxyOutput.hdf5 with that in v_0.1_10_GalaxyOutput.hdf5 and the two should be identical to within close to machine precision (different machine architecture and compilers may cause very minor differences).  You should see a set of 4 plots which essentially show a 1-1 relationship between galaxy properties in the two runs.

Analysis of results
-------------------

There are various jupyter notebooks in anal/ to create analysis plots from the data.  Each has its own YAML file to control which data you want to analyse and what plots to make:

* GalaxyPlots.ipynb / galaxy.yml - this is the most useful one;
* SubhaloPlots.ipynb / subhalo.yml;
* HaloPlots.ipynb / halo.yml.

As yet, these plots are pretty rudimentary and do not contain any observational data - it would be a very useful exercise to add those in.

Code development
----------------

I find it very useful to develop code using a jupyter notebook as that allows for easy inspection of variables should any errors be encountered.  As shipped, the code in L-Galaxies.ipynb and that in L-Galaxies.py should be identical: the latter is a "File: Download as Python" version of the former.

Using 10 graphs is sufficient in most cases and allows for a quick comparison with previous results.  However, the first graph of significant size in Millenniumm file 5 is number 58, so that is a useful one to run also.

If I am making changes that should not affect the results (e.g. when converting to C) then I run Comparison.ipynb frequently.  Note that Comparison.ipynb simply compares galaxies 1-by-1 in each of 2 files: it does not attempt to match galaxies up via halo or subhalo (although it could be modified to do so).  Therefore if the number of galaxies generated by the code changes, then the comparison will break down.

Timing tests
------------

To run timing tests one needs to do a reasonable number of graphs and so I recommend using the full dataset supplied with the code, ie all 1000 graphs in MergerGraphs/Mill5_1000.hdf5.  This will also maintain consistency over the whole of the code development.

Timing should be done with optimisation (python3 -O) switched on.

Production runs
---------------

The code requires input graphs/trees to be in a particular format, as described in the "Input graph file structure" section in the "Developer notes".  Currently, there is one helper application to convert Millennium tree files into this format (code-helper/Millennium_to_pygal.py), so this will need to be run once for each tree file to convert it into the appropriate format.  During that conversion, trees that end prematurely (i.e. have no descendant prior to the final snap) are removed, and, where a descendant has skipped a snapshot, an intermediate halo is created to fill in the gap.  [Note that File 5 contained one, and only one, broken tree link after conversion, so I am not absolutely certain that the conversion code is bug free.]

It is best to keep large data sets distinct from the main py-galaxies distribution.  I find it useful to organise things as follows, but you are of course free to do as you wish.  Create a main simulation directory location <Simulation dir> with subdirectories output and anal.  Then create the following symbolic links:

* :code:`MergerGraphs: ln -s <Merger graph location> <Simulation name>`
* :code:`output: ln -s <Simulation dir>/output <Simulation name>`
* :code:`anal/data: ln -s <Simulation dir>/output <Simulation name>`
* :code:`anal/figs: ln -s <Simulation_dir>/anal <Simulation name>`

One the run is complete, if they are not already in the correct location, then copy over the input YAML file and any log files to the <Simulation dir> directory.

Production runs should be undertaken with the python (.py) script and not with the jupyter notebook, which is very resource hungry:

:code:`python3 -O L-Galaxies.py <input YAML file> > <log file> &`.

If the input YAML file is omitted from the above line, then it will default to input/input.yml.
If the verbosity is greater than 0 then it is not a good idea to copy stdout to the console as simply printing out a list of graphs that have been processed can eat up a lot of time.

