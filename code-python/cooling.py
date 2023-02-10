"""
Functions to cool gas from halos onto the central subhalo, and from the central subhalo onto the galaxy
Calling sequence:
halo_<cooling_model>(halo,subhalo,dt,parameters):
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
sub_<cooling_model>(subhalo,galaxy,dt,parameters):
    Inputs:
    subhalo: a subhalo instance
        The central subhalo of this halo
    galaxy:  a row of the galaxy nparray
        The galaxy associated with the central subhalo
    dt : float
        The length of time for which the cooling occurs
    paramters : an instance of the parameter class
        Contains all the necessary runtime contants and units conversions
    Outputs:
    There are no explicit outputs but the subhalo instance and galaxy row will be updated appropriately
    to reflect the amount of cooling from subhalo gas_hot to galaxy cold_gas.
"""

def halo_lgalaxies(halo,subhalo,dt,parameters):
    """
    Simple model to mimic that of L-Galaxies.  
    The hot gas in halos and subhalos is regarded as a single entity with temperature equal to that of the halo.
    So all we do here is give the hot gas to the subhalo, and
    set the virial temperature of the subhalo to equal that of the halo.
    """
    subhalo.mass_gas_hot += halo.mass_gas_hot
    subhalo.mass_baryons += halo.mass_gas_hot
    halo.mass_gas_hot = 0.
    subhalo.temperature = halo.temperature
    return None

def sub_isothermal(subhalo,galaxy,dt,paramters):
    """
    Isothermal cooling as in the orginal L-galaxies and many other SAMs.
    Currently just a dummy entry.
    """
    return None


