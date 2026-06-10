#!/usr/bin/env python3
"""v67 bath-balance order-of-magnitude estimate [estimate].

Measured (v66/FINDINGS.md section 4, eta=0.5, omega=1.39 seed, Q=482, E=692,
r_Q=3.71, r_half=4.41): peak dQ/dt = -1.72, dE/dt ~ -2.4.
Radiated quanta consistency: (dE/dt)/(dQ/dt) = 1.41 ~ omega_core = 1.39
=> theta radiation carries E/Q = omega (massless quanta at k=omega).

Absorption from an isotropic massless theta bath (c=1): P_abs = sigma_abs * e_bath.
Balance: e_bath* = |dE/dt| / sigma_abs.
"""
import numpy as np

dEdt, dQdt = 2.4, 1.72
omega, E, Q = 1.39, 691.9, 482.2
r_half, r_Q = 4.41, 3.71

print(f"radiation quanta check: (dE/dt)/(dQ/dt) = {dEdt/dQdt:.3f} vs omega = {omega}")
print(f"quality factor Q_fac = omega*E/P_em = {omega*E/dEdt:.0f}  (weak coupling)")
print(f"k_theta = omega = {omega}, lambda_theta = {2*np.pi/omega:.3f}, "
      f"omega*r_Q = {omega*r_Q:.2f} (>1: wavelength-scale ball)")
print()
# (i) geometric upper bound (black sphere)
sig_geo = np.pi*r_half**2
# (ii) low partial-wave unitarity: emission source is l=1-type (THEORY v66 sec.5);
#      sigma_max = pi(2l+1)/k^2 summed l<=1, ~2 effective transverse polarizations
sig_uni = 2*np.pi*(1+3)/omega**2
# (iii) occupation criterion: spontaneous vs stimulated balance when bath mode
#       occupation n~1 at the emission frequencies:
#       e* ~ g * omega^3 * dOmega / (2 pi^2), g = effective real polarizations (2-6),
#       dOmega = emission bandwidth (0.1-1 omega)
e_occ_lo = 2*omega**3*(0.1*omega)/(2*np.pi**2)
e_occ_hi = 6*omega**3*(1.0*omega)/(2*np.pi**2)

for name, sig in (("geometric  (upper bound on sigma)", sig_geo),
                  ("unitarity l<=1, 2 pol           ", sig_uni)):
    print(f"sigma_abs {name} = {sig:6.2f}  -> e_bath* = {dEdt/sig:.4f}")
print(f"occupation criterion (n~1)        -> e_bath* in [{e_occ_lo:.3f}, {e_occ_hi:.3f}]")
lo = dEdt/sig_geo; hi = e_occ_hi
cen = np.sqrt(lo*hi)
print(f"\ncombined bracket: e_bath* in [{lo:.3f}, {hi:.2f}], geometric-mean center ~ {cen:.2f}")
print(f"context: ball core e(0) = 2.25, condensate e_0(A=0.40) = 0.83")
print(f"\nrecommended scan ladder (log-spaced, brackets center):")
print(f"  e_bath in {{0.05, 0.2, 0.8}}  (+ e_bath=0 control)")
print("interpolate the zero crossing of dE/dt(e_bath); slope near 0 gives sigma_abs directly")
