# module to save common variables needed on import

common={}

def save(name,value):
   common[name]=value
   return None

def load(name):
   try:
      return common[name]
   except:
      raise ValueError('name '+name+' does not exist in commons')

