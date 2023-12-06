Profiling
=========

This page outlines the profiling of the code and documents the attempts to improve both CPU and memory usage.

I am really not sure what tools to use for the best here.  One thing is certain, this needs to be done at the command line to avoid bloating the code: in particular `jupyter notebook` is horrendously memory guzzling.  The inbuilt profiling tools are controlled by runtime boolean parameters `b_profile_cpu` and `b_profile_mem`.


CPU time
--------

See `https://realpython.com/python-timer/#the-python-timer-code <https://realpython.com/python-timer/#the-python-timer-code>`_ for hints.

Of the three methods described below, I have found `cProfile` to be the most useful: it does not require editing of the code and it produces a modest-sized profile file.

cProfile
^^^^^^^^

This is a command line tool that can tell you in which functions the code spends most of its time.  It is easy to use as it requires no modification to the code itself: :code:`python -m cProfile -o output_dir/L-Galaxies.prof L-Galaxies.py`.

Then analyse using :code:`python -m pstats output_dir/L-Galaxies.prof` which pops you into an interactive interface.  The recommended basic commands are:

* :code:`strip` – cleans the output;
* :code:`sort (tot|cum)time` – sorts by time spent in functions excluding (tot) or including (cum) subfunctions
* :code:`stats #` – shows the top #  entries.
  
Once you know which function/method is taking most of the time, then you can use a line profiler to determine which individual lines of code are at fault.  This requires an :code:`@profile` directive to be placed in front of the function to be profiled.

* :code:`pip install line_profiler`
* :code:`kernprof -lv -o output_dir/L-Galaxies.py.lprof L-Galaxies.py`
* :code:`python -m line_profiler output_dir/L-Galaxies.py.lprof`

As of 6-Dec-23, prior to any attempt to move any code to C, these were the timings for Millennium File 5 on a MacBook Pro, 2.3 GHz 8-Core Intel Core i9, 16 GB

         594090760 function calls (593190020 primitive calls) in 2382.652 seconds

.. list-table::
   :widths: 10 10 10 10 10 50
   :header-rows: 1
		 
   * - ncalls
     - tottime
     - percall
     - cumtime
     - percall
     - filename:lineno(function)
   * - 19101255
     - 401.648
     - 0.000
     - 459.901
     - 0.000
     - star_formation_and_feedback.py:15(F_gal_form_stars)
   * - 10704577
     - 400.221
     - 0.000
     - 483.569
     - 0.000
     - sfh.py:203(F_sfh_update_bins)
   * - 3128618
     - 210.500
     - 0.000
     - 1704.717
     - 0.001
     - driver.py:55(F_process_halos)
   * - 19114008
     - 160.633
     - 0.000
     - 240.673
     - 0.000
     - star_formation_and_feedback.py:155(F_gal_SNR_feedback)
   * - 5244472
     - 128.201
     - 0.000
     - 246.350
     - 0.000
     - cooling.py:130(F_sub)
   * - 789720
     - 116.177
     - 0.000
     - 140.783
     - 0.000
     - gals.py:237(append)
   * - 789720
     - 106.502
     - 0.000
     - 192.182
     - 0.000
     - driver.py:170(F_update_halos)
   * - 782652
     - 69.800
     - 0.000
     - 113.590
     - 0.000
     - group.py:348(__getitem__)
   * - 5244472
     - 59.468
     - 0.000
     - 80.984
     - 0.000
     - cooling.py:270(F_get_metaldependent_cooling_rate)
   * - 19101255
     - 47.920
     - 0.000
     - 47.920
     - 0.000
     - star_formation_and_feedback.py:112(F_star_formation_unresolved)

The line profiler (not shown here) revealed the following:

* Accidental leaving of units within the main body of the code is very slow.
* Use of column headings in numpy structured arrays takes up most of the time -- could be effectively eliminated by interfacing with C.
* Manipulation of halo and subhalo class instances is slow but does not seem to dominate.

That suggests that significant speed-up could be achieved by:

* Defining a (static?) C struct to hold the internal parameters that are currently attributes of `parameters.py`.
* Defining a (static?) C struct to hold the internal variabless that are currently stored in `commons.py`.
* Associating the `gals` numpy structured array with an equivalent array of C structures (possibly only one at a time rather than an array).
* Changing all the routines below `driver.py` (ie most of them!) to C equivalents.

C_time class
^^^^^^^^^^^^

A simple class that holds a dictionary of numpy records:

* `key` – name of the record
* `value` – n_start, n_stop, cpu_time_start, cpu_time_total
  
with methods:

* `__repr__` – prints out the dictionary.
* `dump(filename)` – saves the dictionary as a pickle file `filename`.
* `start(name)` – adds an entry to the dictionary with key `name`, or reopens an existing one.
* `stop(name)` – accumulates the time spent in cpu_time_total.
  
This routine should be relatively lightweight.  It is used to track the time taken to process each graph.  The following plots show the 




codetiming.Timer
^^^^^^^^^^^^^^^^

`https://pypi.org/project/codetiming/ <https://pypi.org/project/codetiming/>`_

This can be used as a decorator to profile individual python functions.  It is useful but the output seems incredibly bloated.  For example, on processing just 1000 halos it produces an output file that is 300MB in size.  



