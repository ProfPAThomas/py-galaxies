"""
Functions to permit profiling of the code.
"""

import numpy as np
import pickle
from time import perf_counter
import tracemalloc as tm

def conditional_decorator(dec, condition):
   """
   Applies a decorator iff the condition is true.
   """
   def decorator(func):
      if not condition:
         # Return the function unchanged, not decorated.
         return func
      return dec(func)
   return decorator

class C_timer:
   """
   Explicit starting and stopping of timing of code segments.

   Creates a new dictionary entry each time it is called with a different keyword
   """

   def __init__(self):
      # The dictionary to store the timers
      self.timers={}

      # # Dictionaries are stupidly large, so use numpy array to store timing data
      # - this turned out to be exactly the same size.
      # May be possible to do something using np.savez after converting to a 2-D array/
      dtype=np.dtype([
         ('n_start',np.int32),
         ('n_stop',np.int32),
         ('cpu_time_start',np.float32),
         ('cpu_time_total',np.float32)
      ])
      self.template=np.empty(1,dtype)
      self.template['n_start']=0
      self.template['n_stop']=0
      self.template['cpu_time_start']=0.
      self.template['cpu_time_total']=0.
   
      return None

   def __repr__(self):
      for key in self.timers:
         print(key,': ',self.timers[key])
      return ''

   def dump(self,filename):
      with open(filename, 'wb') as f:
         pickle.dump(self.timers, f)
      return None
   
   def start(self,name):
      if name=='': raise AttributeError('No name supplied to C_timer.start()')
      try:
         entry=self.timers[name]
         # Note that the following line will not halt the code as contained in a try statement!
         if entry['n_start']>entry['n_stop']: raise RuntimeError('timer '+name+' already started')
         entry[n_start] += 1
      except:
         # Create new dictionary entry and initialise
         self.timers[name]=self.template.copy()
         #self.timers[name]={}
         entry=self.timers[name]
         entry['n_start']=1
         entry['n_stop']=0
         entry['cpu_time_total']=0.
      entry['cpu_time_start']=perf_counter()
      return None

   def stop(self,name):
      if name=='': raise AttributeError('No name supplied to C_timer.stop()')
      try:
         entry=self.timers[name]
         if entry['n_start']==entry['n_stop']: raise RuntimeError('timer '+name+' already stopped')
         entry['n_stop'] +=1
      except:
         raise ValueError('timer '+name+' does not exist')
      entry['cpu_time_total'] += perf_counter()-entry['cpu_time_start']
      return None

class C_mem:
   """
   Keep track of memory usage.
   All this does is track memory size and peak memory size.

   Creates a new dictionary entry each time it is called with a different keyword
   """

   def __init__(self):
      # The dictionary to store the timers
      self.mem={}

      dtype=np.dtype([
         ('mem_size',np.float32),
         ('mem_peak',np.float32)
      ])
      self.template=np.empty(1,dtype)
      self.template['mem_size']=0.
      self.template['mem_peak']=0.

      # Start profiler
      tm.start()
   
      return None

   def __repr__(self):
      for key in self.mem:
         print(key,': ',self.mem[key])
      return ''

   def dump(self,filename):
      with open(filename, 'wb') as f:
         pickle.dump(self.mem, f)
      return None
   
   def start(self,name):
      if name=='': raise AttributeError('No name supplied to C_mem.start()')
      try:
         entry=self.mem[name]
      except:
         # Create new dictionary entry and initialise
         self.mem[name]=self.template.copy()
         #self.mem[name]={}
      tm.reset_peak()

   def stop(self,name):
      if name=='': raise AttributeError('No name supplied to C_mem.stop()')
      try:
         entry=self.mem[name]
      except:
         raise ValueError('timer '+name+' does not exist')
      tm_size, tm_peak = tm.get_traced_memory()
      entry['mem_size'] = tm_size
      entry['mem_peak'] = tm_peak
      return None
  
