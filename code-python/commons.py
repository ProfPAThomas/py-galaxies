# Module to save common variables that it would be too messy/nit-picky to pass exlicitly.
# Also, some parameters that are needed on import (actually, only b_profile_cpu).
# Other parameters (ie things fixed during the run) are stored in C_parameters.

common={}

def save(key,value):
   common[key]=value
   return None

def load(key):
   try:
      return common[key]
   except:
      raise ValueError('name '+key+' does not exist in commons')

def list(key='_all'):
   if key == '_all':
      for key,value in common.items():
         print(key,': ',value)
   else:
      print(key,': ',common[key])
   return None

# Add some entries in an attempt to prevent documentation build failing

common['b_profile_cpu']=False
common['b_profile_mem']=False
common['b_SFH']=True
common['sfh.n_bin']=1000
