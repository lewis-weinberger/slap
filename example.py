#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Example use of slap.py to generate a mock LAE 
population and find the luminosity function.

    Example:
    --------
        $ python example.py

"""

import numpy as np
from slap import LAEModel
from hmf import MassFunction
from hmf.sample import sample_mf
import matplotlib.pyplot as plt

# Parameters
z = 7.3           # redshift
t_igm = 0.15      # toy IGM transmission fraction
delta_t = 50      # duty cycle parameter
rew_min = 0       # observational minimum Lya REW
lya_min = 2.4e42  # observational minimum Lya luminosity
vol = 1.0e7       # comoving mock survey volume Mpc^3
smallh = 0.678    # dimensionless Hubble parameter

# Create halo population in vol, using SMT halo function
hmf = MassFunction(hmf_model="SMT")
n_halo = int(smallh**3 * vol * np.trapz(x=hmf.m, y=hmf.dndm))
print("\nGenerating", n_halo, "halos in parent population")
halo_mass, _ = sample_mf(n_halo, 10, z=z, hmf_model="SMT")
print("Halo mass: min = {:.1e}, max = {:.1e} Msun/h".format(
    halo_mass.min(), halo_mass.max()))
halo_mass /= smallh # Msun

# Generate LAE population from halo population
lae_pop = LAEModel(z, delta_t, halo_mass, rew_min, lya_min)
n_LAE = lae_pop.lya_lum.size
print(n_LAE, "LAEs generated.")

# Simulated luminosity function (intrinsic)
lf, edges = np.histogram(np.log10(lae_pop.lya_lum), bins=10)
dl = np.diff(edges)
lf = lf / dl / vol
lbins = 0.5 * (edges[1:] + edges[:-1])
lbins = 10**lbins

# Apply observational selection
lya_lum_obs = t_igm * lae_pop.lya_lum
lya_rew_obs = t_igm * lae_pop.lya_rew

mask = np.logical_and(lya_lum_obs > lya_min, lya_rew_obs > rew_min)
lya_lum_obs = lya_lum_obs[mask]
print("LAE luminosity: min = {:.1e}, max = {:.1e} erg/s".format(
    lya_lum_obs.min(), lya_lum_obs.max()))
n_LAE_obs = lya_lum_obs.size
print(n_LAE_obs, "LAEs observed.\n")

# Simulated luminosity function (observed)
lf_obs, edges = np.histogram(np.log10(lya_lum_obs), bins=edges)
lf_obs = lf_obs / dl / vol

# Observed luminosity function from Konno et al. (2014)
k14_x = np.array([42.5, 42.7, 42.9])
k14_y = np.array([-3.6411, -3.9905, -4.6702])
k14_yerr_up = np.array([0.9742, 0.9785, 2.3074])
k14_yerr_down = np.array([0.5422, 0.5434, 0.8442])

k14_x = 10**k14_x
k14_y = (10**k14_y)
k14_yerr_up = k14_yerr_up*k14_y
k14_yerr_down = k14_yerr_down*k14_y

# Plot comparison
fig, ax = plt.subplots(1, 1, figsize=(5,5))

ax.plot(lbins, lf, ls="--", c="grey", label="Model, intrinsic")
ax.plot(lbins, lf_obs, ls="-", c="r", label="Model, T={:.2f}".format(t_igm))
ax.errorbar(k14_x, k14_y, yerr=[k14_yerr_down, k14_yerr_up], ls="",
            marker="o", c="k", label="z=7.3 (Konno et al. 2014)")

ax.legend(loc="upper right")

ax.set_xlabel("Lya luminosity [erg/s]")
ax.set_ylabel("Luminosity function")

ax.set_xscale("log")
ax.set_yscale("log")

fig.tight_layout()
plt.savefig("example.png", bbox_inches="tight")
