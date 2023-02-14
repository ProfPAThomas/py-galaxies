"""
Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy
Calling sequences:
F_halo(halo,subhalo,dt,parameters):
    Inputs:
       halo : a halo instance
       subhalo: a subhalo instance
          The central subhalo of this halo
       dt : float
          The length of time for which the cooling occurs
       paramters : an instance of the parameter class
          Contains all the necessary runtime contants and units conversions
    Outputs:
       There are no explicit outputs but the halo and subhalo instances will be updated appropriately 
       to reflect the amount of cooling from halo gas_hot to subhalo gas_hot.
F_sub(subhalo,galaxy,dt,parameters):
    Inputs:
       subhalo: a subhalo instance
          The central subhalo of this halo
       galaxy:  a row of the galaxy nparray
          The galaxy associated with the central subhalo
       dt : float
          The length of time for which the cooling occurs
       parameters : an instance of the parameter class
          Contains all the necessary runtime contants and units conversions
    Outputs:
       There are no explicit outputs but the subhalo instance and galaxy row will be updated appropriately
       to reflect the amount of cooling from subhalo gas_hot to galaxy cold_gas.
F_cooling_<cooling_model>(...):
    Separated off so that same model can be used for halos and subhalos.  Passes in values rather than (sub)halo instances.
    Would be easy to convert to a C-routine.
    Inputs:
    Outputs:
       mass_cooled: float
          The amount of gas cooled from the (sub)halo.
"""

def F_halo(halo,sub,dt,parameters):
    """
    Cooling of halo onto subhalo.
    """
    if parameters.b_lgalaxies:
        """
        Simple model to mimic that of L-Galaxies.  
        The hot gas in halos and subhalos is regarded as a single entity with temperature equal to that of the halo.
        So all we do here is give the hot gas to the subhalo, and
        set the virial temperature of the subhalo to equal that of the halo.
        """
        sub.mass_gas_hot += halo.mass_gas_hot
        sub.mass_metals_gas_hot += halo.mass_metals_gas_hot
        sub.mass_baryons += halo.mass_gas_hot
        halo.mass_gas_hot = 0.
        halo.mass_metals_gas_hot = 0.
        sub.temperature = halo.temperature
    elif parameters.cooling_model == 'SIS':
        mass_cooled=F_cooling_SIS(halo.mass_gas_hot,halo.mass_metals_gas_hot,halo.temperature,sub.temperature,dt,parameters.c_cooling)
        mass_metals_cooled  = (mass_cooled/halo.mass_gas_hot) * halo.mass_metals_gas_hot
        halo.mass_gas_hot -= mass_cooled
        halo.mass_metals_hot -= mass_metals_cooled
        sub.mass_gas_hot += mass_cooled
        sub.mass_metals_gas_hot += mass_metals_cooled
    else:
        raise valueError('cooling.F_halo: cooling model '+parameters.cooling_model+' not implemented.')
    return None

def F_sub(sub,gal,dt,parameters):
    """
    Cooling of subhalo onto galaxy
    """
    if sub.mass_gas_hot <= parameters.mass_minimum:
        return None
    if parameters.cooling_model == 'SIS':
        gal_temperature=1e4  # Cool down to 1e4 K
        mass_cooled=F_cooling_SIS(sub.mass_gas_hot,sub.mass_metals_gas_hot,sub.temperature,gal_temperature,dt,parameters.c_cooling)
        mass_metals_cooled  = (mass_cooled/sub.mass_gas_hot) * sub.mass_metals_gas_hot
        sub.mass_gas_hot -= mass_cooled
        sub.mass_metals_gas_hot -= mass_metals_cooled
        gal['mass_gas_cold'] += mass_cooled
        gal['mass_metals_gas_cold'] += mass_metals_cooled
    else:
        raise valueError('cooling.F_sub: cooling model '+parameters.cooling_model+' not implemented.')
    return None

def F_cooling_SIS(mass_halo,tau_dyn,mass_gas,mass_metals_gas,temp_start,temp_end,dt):
    """
    Implements the isothermal cooling model as used in L-Galaxies and many other SAMs.
    """
    # Not sure if subhalo virial temperature can ever exceed that of the halo that it is in.
    # If it can, trap out before call to this subroutine, so raise error here
    assert temp_end <= temp_start
    
    # Determine cooling function
    logZ = np.max(-10.,np.log10(mass_metals_gas/mass_gas))
    # Cooling rate per unit density of electrons & ions
    lambda = F_get_metaldependent_cooling_rate(np.log10(temp_start),logZ)
    tau_cool = (temp_start-temp_end)/lambda
    
    # Cooling at constant temperature, but allowing density to vary: see documentation.
    fg0 = mass_gas/mass_halo;
    dt_ratio=dt/tau_dyn;
    tau_ratio=tau_dyn*fg0/tau_cool;
    if tau_ratio <=1:
	  fg=fg0/(1+0.5*np.sqrt(tau_ratio)*dt_ratio)**2
      else:
	  teq_ratio=np.log(tau_ratio)
	  if dt_ratio <= teq_ratio:
	      fg=fg0*np.exp(-dt_ratio)
	  else:
	      fg=fg0/(tau_ratio*(1+0.5*(dt_ratio-teq_ratio))**2)
    
    return fg*mass_halo

def F_get_metaldependent_cooling_rate(logT,logZ):
    """
    Returns the cooling function, ie the cooling rate per unit density of electrons and ions.
    Assumes that the cooling function is tabulated in code units.
    """
    # Needs to read in/know cooling function
    # Needs to read in/know table limits
    # Fix tables so that cannot fall outside range to save checks here
    # The following indices are at the bottom end of the range
    i_T=np.argmax(np.where(logT>logT_table))
    fracT=(logT-logT_table[i_T])/(logT_table[i_T+1]-logT_table[i_T})
    assert: 0 <= frac <=1
    i_Z=np.argmax(np.where(logZ>logZ_table))
    fracZ=(logZ-logZ_table[i_Z])/(logZ_table[i_Z+1]-logZ_table[i_Z})
    # Interpolate in 2-d
    lambda0 = fractT*lambda_table[i_T+1,i_Z]+(1-fracT)*lambda_table[i_T,i_Z]
    lambda1 = fractT*lambda_table[i_T+1,i_Z+1]+(1-fracT)*lambda_table[i_T,i_Z+1]
    lambda = fracZ*lambda1+(1-fracZ)*lambda0
    return lambda

    
