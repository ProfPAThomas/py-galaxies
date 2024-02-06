"""
Module to save common parameters that are needed during module importing.

I have deliberately separated adding new entries from updating existing ones to prevent bugs.
The 'lock' key can be set to True to prevent the addition on new entries.

This module is a bit OTT but formerly had more uses in the code than it does at present.
"""

common={'lock': False}

# Add some entries in an attempt to prevent documentation build failing

common['b_profile_cpu']=False
common['b_profile_mem']=False
common['b_SFH']=True
common['sfh_n_bin']=1000

def add(key,value):
   """
   Saves new quantities in a dictionary.

   Attributes
   ----------
   key : str
      Name of the quantity to be saved.
   value : obj (various)
      Value of the quantity to be stored.
   """
   if common['lock']:
      raise ValueError('common directory is locked to new entries (because C struct has already been generated)')
   elif key in common:
      raise ValueError('key '+key+' already present in common dictionary: use update() for existing entries')
   else:
      common[key]=value
   return None

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

def lock(b_lock=True):
   """
   Locks/unlocks the common directory to the addition of new entries.
   Needed to prevent accidental addition of entries after the corresponding C struct is generated.

   Attributes
   ----------
   b_lock : bool
      Whether to lock or unlock the common dictionary to new entries
   """
   common['lock']=b_lock
   return None

def update(key,value):
   """
   Updates existing quantity in a dictionary.

   Attributes
   ----------
   key : str
      Name of the quantity to be saved.
   value : obj (various)
      Value of the quantity to be stored.
   """
   if key in common:
      common[key]=value
   else:
      raise ValueError('key '+key+' not present in common dictionary: use add() for new entries')
   return None
