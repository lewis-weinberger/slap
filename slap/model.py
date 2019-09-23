#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate an LAE population from a dark matter halo population, based on 
the methods outlined in Weinberger et al. (2019).

Examples:
--------
    Import the module:
        >>> import model

    Generate LAE population:
        >>> lae_pop = slap.LAEModel(z, delta_t, halo_mass, rew_min, lya_min)

    LAE properties:
        >>> halo_ids = lae_pop.halo_ids
        >>> halo_mass = lae_pop.halo_mass
        >>> uv_lum = lae_pop.uv_lum
        >>> lya_lum = lae_pop.lya_lum
        >>> lya_rew = lae_pop.lya_rew

"""

import warnings
import numpy as np
from numpy import trapz
from scipy.integrate import quad
from scipy.optimize import bisect
from scipy.interpolate import interp1d
from hmf import MassFunction

#------------------------------------------------------------------------------
# PARAMETERS

# CONSTANTS AND CONVERSIONS
YR = 365.25*24*60*60      # seconds in a year
MYR = 1e6*YR              # seconds in a Myr
PC = 3.086e16             # metres in a parsec
MPC = 1e6*PC              # metres in a Mpc
tconv = (MPC/1e3)/MYR     # convert inverse Hubble to Myr

# COSMOLOGY
omega_L = 0.692           # Cosmological constant
omega_m = 0.308           # Matter
smallh = 0.678            # Hubble constant h=H0/100

# z=6 LBG LUMINOSITY FUNCTION
# See Bouwens et al. (2015) "all fields"
phi0_z6 = 0.50e-3         # Mpc^-3
mag0_z6 = -20.94          # AB magnitude
alpha_z6 = -1.87          # dimensionless

# LAE EW DISTRIBUTION
l_alpha = 1216            # Lya wavelength, Angstroms
l_UV = 1600               # UV wavelength, Angstroms
nu_alpha = 2.47e15        # Lya frequency, Hz
betaUV = -1.7             # UV slope
dijC = (nu_alpha/l_alpha)*(l_UV/l_alpha)**(-betaUV-2)

# z=6 HALO MASS FUNCTION
Mmin = np.log10(5e8*smallh)
Mmax = np.log10(1e15*smallh)
SMT_HMF = MassFunction(Mmin=Mmin, Mmax=Mmax, z=6, hmf_model="SMT",
    cosmo_model="Planck13")

#------------------------------------------------------------------------------

def vectorize_if_needed(func, x):
    """Helper function to vectorise functions on array arguments.
    
    Args:
    -----
        func (function): input function.
        x (scalar or array_like): function argument(s).

    Returns:
    --------
        out (scalar or array_like): function applied to argument(s).

    """

    if np.isscalar(x):
        return func(x)

    vfunc = np.vectorize(func)
    return vfunc(x)

    
def Hub(z):
    """Hubble parameter.

    Args:
    -----
        z (float or array_like): redshift.

    Returns:
    --------
        Hub (float or array_like): Hubble parameter in km/s/Mpc.

    """

    omega_k = 1 - omega_m - omega_L
    H0 = 100*smallh
    H_z = H0 * np.sqrt( (omega_m)*(1 + z)**3 + omega_L + omega_k*(1 + z)**2 )

    return H_z 


def tHub(z):
    """Cosmological age.

    Args:
    -----
        z (float or array_like): redshift.

    Returns:
    --------
        tHub (float or array_like): age in Myr.
    
    """
    
    func = lambda x: 1/((1+x)*Hub(x))
    t_z_func = lambda zp: quad(func, zp, np.inf)[0]
    t_z = vectorize_if_needed(t_z_func, z)

    return t_z*tconv


def dtHub(z_end, z_start):
    """Cosmological time interval.

    Args:
    -----
        z_end (float or array_like): end redshift.
        z_start (float or array_like): start redshift.

    Returns:
    --------
        dtHub (float): time interval in Myr
    
    """
    
    func = lambda x: 1/((1+x)*Hub(x))
    dt_func = lambda ze: quad(func, z_start, ze)[0]
    dt = vectorize_if_needed(dt_func, z_end)
    
    return dt*tconv


def smt_hmf(znew, log=False):
    """Sheth-Mo-Tormen halo mass function.

    Args:
    -----
        znew (float): redshift for halo mass function.
        log (bool): if True use logarithmic mass bins.

    Returns:
    --------
        mass (array_like): mass bins.
        hmf (array_like): halo mass function for mass bins.

    """
    
    SMT_HMF.update(z=znew)
    m = SMT_HMF.m                # Msun/h
    if log:
        hmf = SMT_HMF.dndlog10m  # h^3 Mpc^-3
    else:
        hmf = SMT_HMF.dndm       # h^4 Msun^-1 Mpc^-3

    return m, hmf


def schecter_M(mag, phi0, mag0, alpha):
    """Schecter luminosity function.
    Note: function of magnitudes not luminosities.
    
    Args:
    -----
        mag (float or array_like): magnitude bin(s).
        phi0 (float): magnitude reference value.
        mag0 (float): luminosity function reference value.
        alpha (float): slope.

    Returns:
    --------
        n_M (float or array_like): number density in Mpc^-3.

    """
    magpow = 10**(0.4*(mag0 - mag))
    n_M = 0.4*np.log(10)*phi0 * magpow**(alpha+1) * np.exp(-magpow)

    return n_M


def edc(Mh, z, delta_t=50):
    """Effective duty cycle, see Trenti et al. 2010.
    
    Args:
    -----
        Mh (float or array_like): halo mass(es) in Msun/h.
        z (float): redshift.
        delta_t (float): duty cycle parameter.

    Returns:
    --------
        eps (float or array_like): duty cycle fraction.

    """
    
    guess = [z-1,z+5]
    test = lambda zdt: dtHub(zdt,z) - delta_t
    z_dt = bisect(test, guess[0], guess[1])
   
    m, hmf = smt_hmf(z) # Msun/h and  (Msun/h)^-1 (Mpc/h)^-3
    m_dt, hmf_dt = smt_hmf(z_dt)
    numer_func = lambda M: (trapz(y=hmf[m>=M],x=m[m>=M])
                          - trapz(y=hmf_dt[m_dt>=M],x=m_dt[m_dt>=M]))
    denom_func = lambda M: trapz(y=hmf[m>=M],x=m[m>=M])

    numer = vectorize_if_needed(numer_func, Mh)
    denom = vectorize_if_needed(denom_func, Mh)
    denom = np.where(denom <= 0, 1e-20, denom)
    eps = np.divide(numer, denom)

    return eps


def edc_quick(Mh, z, delta_t):
    """Interpolate effective duty cycle from sparser sampling.
    
    Args:
    -----
        Mh (float or array_like): halo mass(es) in Msun/h
        z (float): redshift
        delta_t (float): duty cycle parameter

    Returns:
    --------
        eps (float or array_like): duty cycle fraction

    """
    
    mass, _ = smt_hmf(z) # Msun/h and  (Msun/h)^-1 (Mpc/h)^-3
    edc_func = interp1d(mass, edc(mass, z, delta_t), 
        fill_value='extrapolate')

    return edc_func(Mh)


def REWc(Muv,z):
    """Characteristic REW.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        Muv (array_like): UV magnitude.
        z (float): redshift.

    Returns:
    --------
        REWc (array_like): characteristic REW.

    """
    
    REWc0 = 23
    deltaM = Muv + 21.9
    deltaz = z - 4
    dREW_dM = 7
    dREW_dz = 6
    REWc = REWc0 + dREW_dM*deltaM + dREW_dz*deltaz
    REWc_19 = REWc0 + dREW_dM*(-19+21.9) + dREW_dz*deltaz
    REWc = np.where(Muv > -19, REWc_19, REWc)

    return REWc


def a1(Muv):
    """Minimum REW.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        Muv (array_like): UV magnitudes.

    Returns:
    --------
        a1 (array_like): minimum REWs.

    """
    
    a1 = 20*np.ones(shape=Muv.shape)
    a1 = np.where(Muv >= -21.5, 20 - 6*(Muv + 21.5)**2, a1)
    a1 = np.where(Muv > -19, -17.5, a1) 

    return a1


def rewnorm(Muv, z, rewmax=400):
    """REW distribution normalisation.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        Muv (array_like): UV magnitudes.
        z (float): redshift.
        rewmax (float): maximum REW.

    Returns:
    --------
        N (array_like): normalisations.

    """
    
    rewc = REWc(Muv, z)
    N = (np.exp(a1(Muv)/rewc) - np.exp(-rewmax/rewc))**(-1)
    N /= rewc
    
    return N


def invCDF_ew(P, Muv, z):
    """Inverse cumulative probability density function.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        P (array_like): probability density.
        Muv (array_like): UV magnitudes.
        z (float): redshift.

    Returns:
    --------
        REW (array_like).

    """
    
    rewc = REWc(Muv, z)
    N = rewnorm(Muv, z)
    C = N*rewc
    C1 = np.exp(a1(Muv)/rewc)

    return -rewc*np.log(C1-P/C)


def PDF_ew(W, Muv, z):
    """REW probability density function.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        W (array_like): REWs.
        Muv (array_like): UV magnitudes.
        z (float): redshift.

    Returns:
    --------
        P (array_like): PDF(W|Muv).

    """
    
    assert(W.shape==Muv.shape)
    rewc = REWc(Muv, z)
    N = rewnorm(Muv, z)
    
    return N*np.exp(-W/rewc)


def CDF_ew(W, Muv, z):
    """Cumulative probability density function.
    See Dijkstra & Wyithe (2010).

    Args:
    -----
        W (array_like): REWs.
        Muv (array_like): UV magnitudes.
        z (float): redshift.

    Returns:
    --------
        P (array_like): CDF(W|Muv).

    """
    
    assert(W.shape==Muv.shape)
    rewc = REWc(Muv, z)
    N = rewnorm(Muv, z)
    C = np.exp(a1(Muv)/rewc)
    
    return rewc*N*(C-np.exp(-W/rewc))


def sample_dist(Muv, z):
    """Sampling of probability distribution.

    Args:
    -----
        Muv (float or array_like): UV magnitude(s).
        z (float): redshift

    Returns:
    --------
        REW sample (float or array_like): sample of REW(s).
    
    """
    
    if np.isscalar(Muv):
        urand = np.random.uniform(low=0, high=1, size=1)
    else:
        npoints = Muv.size
        urand = np.random.uniform(low=0, high=1, size=npoints)
    sample = invCDF_ew(urand, Muv, z)

    return sample


def abundance_match(delta_t):
    """Abundance match halo mass function to UV luminosity function at z=6.

    Args:
    -----
        delta_t (float): duty cycle parameter in Myr.

    Returns:
    --------
        lum_func (func): mapping from halo mass to LBG UV luminosity.

    """

    mass, hmf = smt_hmf(5.9) # Msun/h and (Msun/h)^-1 (Mpc/h)^-3
    edc_mass = edc(mass, 5.9, delta_t)
    lum = np.zeros(shape=mass.shape)
    for m in range(mass.size):
        M = mass[m]
        edc_M = edc_mass[m]
        mass_integral = edc_M*np.trapz(y=hmf[mass>=M]*smallh**3, x=mass[mass>=M]) 
        # units = Mpc^-3

        def find_map(desired_mag):
            mag_space = np.linspace(-27, desired_mag, num=100)
            phi_space = schecter_M(mag_space, phi0_z6, mag0_z6, alpha_z6)
            lum_integral = np.trapz(y=phi_space, x=mag_space) # units Mpc^-3

            return lum_integral - mass_integral

        guess = [-10, -27]
        lum_mass = bisect(find_map, guess[0], guess[1], disp=True)
        lum[m] = lum_mass

    mass /= smallh # units of Msun
    lum_func = interp1d(mass, lum, fill_value='extrapolate')
    # maps from [Msun] to [AB magnitude]

    return lum_func


def generate_LBGs(lum_func, halo_mass, redshift, delta_t):
    """Use UV-Mh mapping to generate LBG population.

    Args:
    -----
        lum_func (func): UV luminosity mapping.
        halo_mass (array_like): halo masses in Msun.
        redshift (float): redshift.
        delta_t (float): duty cycle parameter.

    Returns:
    --------
        UV_lum (array_like): UV luminosities of LBG sample.
        halo_mass (array_like): halo masses of LBG sample.
        halo_ids (array_like): halo ids of LBG sample.

    """

    halo_ids = np.arange(halo_mass.size)
    UV_lum = lum_func(halo_mass)

    edcs = edc_quick(halo_mass*smallh, redshift, delta_t)

    duty = np.ones_like(halo_mass,dtype=np.int32)
    for n in range(duty.size):
        duty[n] = np.random.choice(2, p=[1-edcs[n],edcs[n]])
    duty_mask = np.where(duty==1, True, False)

    halo_mass = halo_mass[duty_mask]
    UV_lum = UV_lum[duty_mask]
    halo_ids = halo_ids[duty_mask]

    return UV_lum, halo_mass, halo_ids


def generate_LAEs(UV_lum, halo_mass, halo_ids, redshift, rew_min, lya_min):
    """Use REW distribution to generate LAE population.

    Args:
    -----
        lum_func (func): UV luminosity mapping.
        halo_mass (array_like): halo masses in Msun.
        redshift (float): redshift.
        delta_t (float): duty cycle parameter.

    Returns:
    --------
        UV_lum (array_like): UV luminosities of LBG sample.
        halo_mass (array_like): halo masses of LBG sample.
        halo_ids (array_like): halo ids of LBG sample.

    """

    # Generate EW distribution
    ew_dist = sample_dist(UV_lum, redshift)

    # Calculate Lya luminosity
    flya = dijC * ew_dist
    lya_lum = flya * 10**(-0.4*(UV_lum-51.6))     # erg/s

    # Observational thresholds
    mask = np.logical_and(ew_dist > rew_min, lya_lum > lya_min)
    lya_lum = lya_lum[mask]
    halo_mass = halo_mass[mask]
    halo_ids = halo_ids[mask] 
    UV_lum = UV_lum[mask]
    ew_dist = ew_dist[mask]

    return halo_ids, halo_mass, UV_lum, ew_dist, lya_lum


class LAEModel:
    """The object that generates an LAE population.
    
    Attributes:
    -----------
        halo_ids (array_like): halo IDs of the LAE population.
        halo_mass (array_like): host halo masses in Msun.
        uv_lum (array_like): UV luminosities in erg/s.
        lya_rew (array_like): Lya REWs in Angstroms.
        lya_lum (array_like): Lya luminosities in erg/s.
    
    """

    def __init__(self, redshift, delta_t, mass_h, rew_min=0.0, lya_min=1.0e42):
        """Initialise the LAE model.
        
        Args:
        -----
            redshift (float): mean redshift of the LAE population.
            delta_t (float): duty cycle parameter in Myr.
            mass_h (array_like): array of halo masses in Msun.
            rew_min (float): minimum observationally detectable REW in 
                Angstroms.
            lya_min (float): minimum observationally detectable luminosity
                in erg/s.

        """
        if mass_h.min() < (10**Mmin)/smallh:
            print("Note: UV abundance matching for haloes",
                  "above M = {:.3e} Msun/h".format(10**Mmin))
            warnings.warn("Input halo mass below abundance matching regime!")
        
        uv_func = abundance_match(delta_t)
        print("Abundance matching at z=6 completed!")

        lbg_lum, lbg_mass, lbg_ids = generate_LBGs(uv_func, mass_h,
                                                   redshift, delta_t)
        print("LBG population generated!")

        lae_data = generate_LAEs(lbg_lum, lbg_mass, lbg_ids, redshift,
                                 rew_min, lya_min)
        print("LAE population generated!")

        self.halo_ids = lae_data[0]
        self.halo_mass = lae_data[1]
        self.uv_lum = lae_data[2]
        self.lya_rew = lae_data[3]
        self.lya_lum = lae_data[4]
