"""
Module to save common variables that it would be too messy/nit-picky to pass exlicitly.
Also, some parameters that are needed on import (actually, only b_profile_cpu).
Other parameters (ie things fixed during the run) are stored in C_parameters.
"""

common={}

def save(key,value):
   """
   Saves quantities in a dictionary.

   Attributes
   ----------
   key : str
      Name of the quantity to be saved.
   value : obj (various)
      Value of the quantity to be stored.
   """
   common[key]=value
   return None

def load(key):
   """
   Reads quantities from a dictionary.

   Attributes
   ----------
   key : str
      Name of the quantity to be retrieved.
   """
   try:
      return common[key]
   except:
      raise ValueError('name '+key+' does not exist in commons')

def list(key='_all'):
   """
   Lists one or all quantities stored in the common dictionary.

   Attributes
   ----------
   key : str
      Name of the quantity to be listed (or all if none specified).
   """
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
