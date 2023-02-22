import numpy as np
import astropy.units as u
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
   
class C_cooling:
    """
    Class to generate and store the cooling tables.
    """
    def __init__(self,parameters):
        """
        Reads in the cooling function and tabulates it.
        Normalises it such that tau_cool=T/Lambda in code units.
        Also, the temperature lookup is adjusted to code units.
        Metallicities are always assumed to be absolute (ie not relative to solar).
        """
        import astropy.constants as c
        import astropy.units as u
        data=np.load(parameters.cooling_function_file)
        self.log10_T_table=data['log10_T']+np.log10(u.K/parameters.units_temperature_internal)  # Converted to code units
        self.log10_Z_table=data['log10_Z']
        # The first correction term converts to code units, and the second provides a dimensionless factor to reduce the cooling time
        # formula to dt=T/Lambda
        self.log10_Lambda_table=data['log10_Lambda'] + np.log10(
            u.erg*u.cm**3*parameters.units_time_internal / \
            (u.s*parameters.units_energy_internal*parameters.units_length_internal**3) ) - \
            np.log10(parameters.c_cooling)
        print('C_cooling verbosity =',parameters.verbosity)
        if parameters.verbosity >=4: print('log10_Lambda_table =',self.log10_Lambda_table)

def F_halo(halo,sub,dt,cooling_table,parameters):
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
        sub.mass_baryon += halo.mass_gas_hot
        halo.mass_gas_hot = parameters.mass_minimum # Just so it appears in the plots as a non-zero value.
        halo.mass_metals_gas_hot = parameters.mass_minimum*parameters.base_metallicity
        sub.temperature = halo.temperature
    elif parameters.cooling_model == 'SIS':
        mass_cooled=F_cooling_SIS(halo.mass,halo.tau_dyn,halo.half_mass_radius,halo.mass_gas_hot,halo.mass_metals_gas_hot,
                                  halo.temperature,sub.temperature,dt,cooling_table)
        mass_metals_cooled  = (mass_cooled/halo.mass_gas_hot) * halo.mass_metals_gas_hot
        halo.mass_gas_hot -= mass_cooled
        halo.mass_metals_hot -= mass_metals_cooled
        sub.mass_gas_hot += mass_cooled
        sub.mass_metals_gas_hot += mass_metals_cooled
    else:
        raise valueError('cooling.F_halo: cooling model '+parameters.cooling_model+' not implemented.')
    return None

def F_sub(sub,gal,dt,cooling_table,parameters):
    """
    Cooling of subhalo onto galaxy
    """
    if sub.mass_gas_hot <= parameters.mass_minimum:
        return None
    if parameters.cooling_model == 'SIS':
        gal_temperature=1e4*u.K/parameters.units_temperature_internal  # Cool down to 1e4 K
        mass_cooled=F_cooling_SIS(sub.mass,sub.tau_dyn,sub.half_mass_radius,sub.mass_gas_hot,sub.mass_metals_gas_hot,
                                  sub.temperature,gal_temperature,dt,cooling_table)
        mass_metals_cooled  = (mass_cooled/sub.mass_gas_hot) * sub.mass_metals_gas_hot
        sub.mass_gas_hot -= mass_cooled
        sub.mass_metals_gas_hot -= mass_metals_cooled
        gal['mass_gas_cold'] += mass_cooled
        gal['mass_metals_gas_cold'] += mass_metals_cooled
    else:
        raise valueError('cooling.F_sub: cooling model '+parameters.cooling_model+' not implemented.')
    return None

def F_cooling_SIS(mass,tau_dyn,half_mass_radius,mass_gas,mass_metals_gas,temp_start,temp_end,dt,cooling_table):
    """
    Implements the isothermal cooling model as used in L-Galaxies and many other SAMs.
    """
    # Not sure if subhalo virial temperature can ever exceed that of the halo that it is in.
    # If it can, trap out before call to this subroutine, so raise error here
    if temp_end >= temp_start:
        print('cooling:F_cooling_SIS: mass, temp_start, temp_end =',mass,temp_start,temp_end)
        return mass_gas
    
    # Determine cooling function
    log10_Z = max(-10.,np.log10(mass_metals_gas/mass_gas))
    # Cooling rate per unit density of electrons & ions
    Lambda = F_get_metaldependent_cooling_rate(np.log10(temp_start),log10_Z,cooling_table)
    tau_cool = half_mass_radius**3*(temp_start-temp_end)/(mass*Lambda)
    
    # Cooling at constant temperature, but allowing density to vary: see documentation.
    # The gas fraction here is relative to the halo mass (= 200 times the critical density in Millennium/SIS)
    fg0 = mass_gas/mass;
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
    #if 100<mass<110.: print('cooling:F_cooling_SIS: temp_start, tau_dyn, tau_cool, dt_ratio, fg0, fg =',
    #      temp_start,tau_dyn,tau_cool,dt_ratio,fg0,fg)
    
    return fg*mass

def F_get_metaldependent_cooling_rate(log10_T,log10_Z,cooling_table):
    """
    Returns the cooling function, ie the cooling rate per unit density of electrons and ions.
    Assumes that the cooling function is tabulated in code units.
    """
    # Needs to read in/know cooling function
    # Needs to read in/know table limits
    # Fix tables so that cannot fall outside range to save checks here
    # The following indices are at the bottom end of the range
    log10_T_table = cooling_table.log10_T_table
    i_T=np.argmax(np.where(log10_T>log10_T_table))
    fracT=(log10_T-log10_T_table[i_T])/(log10_T_table[i_T+1]-log10_T_table[i_T])
    assert 0 <= fracT <=1
    log10_Z_table = cooling_table.log10_Z_table
    i_Z=np.argmax(np.where(log10_Z>log10_Z_table))
    fracZ=(log10_Z-log10_Z_table[i_Z])/(log10_Z_table[i_Z+1]-log10_Z_table[i_Z])
    if not 0 <= fracZ <=1:
        print('log10_Z_table =',log10_Z_table)
        print('log10_Z =',log10_Z)
        print('i_Z =',i_Z)
    assert 0 <= fracZ <=1, 'base_metallicity below minimum of cooling tables'
    # Interpolate in 2-d
    log10_Lambda_table = cooling_table.log10_Lambda_table
    log10_Lambda0 = fracT*log10_Lambda_table[i_Z,i_T+1]+(1-fracT)*log10_Lambda_table[i_Z,i_T]
    log10_Lambda1 = fracT*log10_Lambda_table[i_Z+1,i_T+1]+(1-fracT)*log10_Lambda_table[i_Z+1,i_T]
    log10_Lambda = fracZ*log10_Lambda1+(1-fracZ)*log10_Lambda0
    return 10.**log10_Lambda

    
