# Kristi, sent/received March 3rd, 2022

import numpy as np

__all__ = ['sfr_to_mwa', 'logsfr_ratios_to_masses', 'logsfr_ratios_to_sfrs',
           'build_agebins', 'zfrac_to_sfrac', 'zfrac_to_masses',
           'zfrac_to_sfr', 'sfrac_to_logsfr_ratios']

def logsfr_ratios_to_masses(logmass=0, logsfr_ratios=None, agebins=None, **extras) :
    '''This converts from an array of log_10(SFR_j / SFR_{j+1}) and a value of
    log10(\Sum_i M_i) to values of M_i.  j=0 is the most recent bin in lookback
    time.
    '''
    nbins = agebins.shape[0]
    sratios = 10**np.clip(logsfr_ratios, -100, 100) # numerical issues...
    dt = (10**agebins[:, 1] - 10**agebins[:, 0])

    if sratios.ndim>1:
        coeffs = np.array([ (1. / np.prod(sratios[:,:i], axis=1)) * (np.prod(dt[1: i+1]) / np.prod(dt[: i])) for i in range(nbins)])
        m1 = (10**logmass) / np.sum(coeffs, axis=0)
        try:
            m1 = (m1.T * coeffs.T).T
        except:
            m1 = m1.T * coeffs
    else:
        coeffs = np.array([ (1. / np.prod(sratios[:i]))           * (np.prod(dt[1: i+1]) / np.prod(dt[: i])) for i in range(nbins)])
        m1 = (10**logmass) / coeffs.sum()
        m1 = m1 * coeffs

    return m1

def logsfr_ratios_to_sfrs(logmass=0, logsfr_ratios=None, agebins=None, **extras) :
    '''Convenience function
    '''
    masses = logsfr_ratios_to_masses(logmass=logmass, logsfr_ratios=logsfr_ratios, agebins=agebins)
    dt = (10**agebins[:, 1] - 10**agebins[:, 0])

    if masses.ndim>1:
        nx,ny = np.shape(masses)[1], len(agebins)
        dt = np.tile(dt, (nx,1)).T
        sfrs = masses / dt
        if np.shape(logsfr_ratios)+(0,1) != np.shape(sfrs) : sfrs = sfrs.T
        return sfrs
    else:
        return masses / dt


# Define fixed age bins
def build_agebins(redshift=1.2, ncomp=10, tuniv=None, tlims_first=[0.03,0.1,0.5,1.], **extras) :
    '''
    Define fixed age bins
    redshift = 1.2 is the median redshift of the GOGREEN spectroscopic sample
    Age bins are spaced:
        0 < t < 30 Myr
        30 < t < 100 Myr
        100 < t < 500 Myr
        500 Myr < t < 1 Gyr
        ncomp-4 linearly spaced bins
        0.95*t_universe < t < t_universe
    '''
    from prospect.sources.constants import cosmo # In my case WMAP9
    if tuniv is None: tuniv = cosmo.age(redshift).value
    tbinmax = (tuniv * 0.95)
    agelims = [1e-9] + tlims_first[:-1] + np.linspace(tlims_first[-1], tbinmax, ncomp-len(tlims_first)).tolist() + [tuniv]

    agelims = np.log10(np.array(agelims) * 1e9)
    agebins = np.array([agelims[:-1], agelims[1:]])
    agebins = agebins.T

    #agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
    return agebins

# Functions for the SFH prior
def zfrac_to_sfrac(z_fraction=None, **extras) :
    z_fraction = np.atleast_2d(z_fraction)
    shape = list(z_fraction.shape)
    shape[-1] += 1
    sfr_fraction = np.zeros(shape)
    sfr_fraction[:, 0] = 1. - z_fraction[:, 0]
    for i in range(1, shape[-1]-1) :
        sfr_fraction[:,i] = np.prod(z_fraction[:,:i], axis=-1) * (1. - z_fraction[:,i])
    sfr_fraction[:,-1] = 1. - np.sum(sfr_fraction[:,:-1], axis=-1)
    sfr_fraction = np.squeeze(sfr_fraction)
    return sfr_fraction

# functions copied from prospect.models.transforms (and adapeted for 2D arrays)
def zfrac_to_masses(z_fraction=None, total_mass=1, agebins=None, sfr_fraction=None, **extras) :
    if sfr_fraction is None: sfrac = zfrac_to_sfrac(z_fraction)
#     sfr_fraction = np.atleast_2d(np.copy(sfr_fraction))
    else: sfrac = np.copy(sfr_fraction)

    # convert to mass fractions
    time_per_bin = np.diff(np.power(10, agebins), axis=-1)[:,0]
    sfrac *= time_per_bin
    mtot = np.atleast_1d(sfrac.sum(axis=-1))
    mass_fraction = sfrac / mtot[:,None]
    masses = np.atleast_2d(total_mass) * mass_fraction.T
    return masses.T

def zfrac_to_sfr(total_mass=1, z_fraction=None, agebins=None, sfr_fraction=None, **extras) :
    time_per_bin = np.diff(np.power(10, agebins), axis=-1)[:,0]
    masses = zfrac_to_masses(total_mass=total_mass, z_fraction=z_fraction, agebins=agebins, sfr_fraction=sfr_fraction)
    return masses / time_per_bin

def sfrac_to_logsfr_ratios(sfrac) :
    N = len(sfrac)
    logsfr_ratios = np.full(N-1, np.nan)

    for i in range(N-1) :
        logsfr_ratios[i] = np.log10(sfrac[i]/sfrac[i+1])
    return logsfr_ratios

def sfr_to_mwa(agebins, sfrs, **extras) :
    '''
    t_mw = (int_0^t_obs t SFR(t) dt) / (int_0^t_obs SFR(t) dt)
    '''
    agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
    mts = np.median(agebins_Gyr, axis=1)
    dts = np.diff(agebins_Gyr)[:,0]
    if sfrs.ndim>1:
        t_mws = np.sum(sfrs * mts * dts, axis=1) / np.sum(sfrs * dts, axis=1)
    else:
        t_mws = np.sum(sfrs * mts * dts) / np.sum(sfrs * dts)
    return t_mws

def dust2_to_Av(dust2) :
    return 2.5*np.log10(np.exp(1)) * dust2

def chain_to_sfr(chain, theta_index, agebins, norm_by_mass=True, mass=None, **extras) :
    if 'logsfr_ratios' in theta_index.keys() :
        logsfr_ratios = chain[:, theta_index['logsfr_ratios']]
        if norm_by_mass: logmass = 0
        elif mass is None: logmass = (chain[:,theta_index['logmass']]).T
        else: logmass = mass
        sfrs = logsfr_ratios_to_sfrs(logmass=logmass, logsfr_ratios=logsfr_ratios, agebins=agebins)

    elif 'z_fraction' in theta_index.keys() :
        z_fraction = chain[:, theta_index['z_fraction']]
        if norm_by_mass: total_mass = 1
        elif mass is None: total_mass = (chain[:,theta_index['total_mass']]).T
        else: total_mass = mass
        sfrs = zfrac_to_sfr(total_mass=total_mass, z_fraction=z_fraction, agebins=agebins)

    return sfrs

def chain_to_mwa(chain, theta_index, agebins, **extras) :
    sfrs = chain_to_sfr(chain, theta_index, agebins, **extras)
    return sfr_to_mwa(agebins, sfrs, **extras)

def chain_to_masses(chain, theta_index, agebins, **extras) :
    if 'logsfr_ratios' in theta_index.keys() :
        logsfr_ratios = chain[:,theta_index['logsfr_ratios']]
        masses = logsfr_ratios_to_masses(logmass=0, logsfr_ratios=logsfr_ratios, agebins=agebins).T

    elif 'z_fraction' in theta_index.keys() :
        z_fraction = chain[:,theta_index['z_fraction']]
        masses = zfrac_to_masses(total_mass=1, z_fraction=z_fraction, agebins=agebins)

    return masses

def sfr_dt(ssfrs, time_Gyr, agebins, **extras) :
    agebins_Gyr = np.power(10, agebins-9)
    dt = np.diff(agebins_Gyr)
    itime = np.argmin(np.abs(agebins_Gyr[:,0] - time_Gyr))
    x = np.average(ssfrs[:,:itime], weights=dt[:itime,0], axis=1)
    return x

def cumulate_masses(masses) :
    return np.cumsum(masses[:,::-1], axis=1)[:,::-1]
