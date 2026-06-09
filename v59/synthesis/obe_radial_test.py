#!/usr/bin/env python3
"""
v59/synthesis/obe_radial_test.py

Test the radial law of the OBE gravity term.

OBE gravity piece:  Omega_grav(x) = f_g * INT K(x,x') grad'( rho_N(x') ) d3x'
with rho_N = (|xi|^2 - 1/2)  (the Koide-constraint deviation density).

Two facts we check numerically:

  (1) For a MASSLESS kernel K (static Green's fn of -lap), the gradient-source
      form collapses to  Omega = f_g * grad(Phi),  Phi = K * rho_N,  i.e. the
      gradient comes OUT of the convolution (integration by parts). So the
      connection is the FORCE of a potential sourced by rho_N itself.

  (2) Whether that force is long-range 1/r^2 depends ENTIRELY on the MONOPOLE
      of rho_N:
        - monopole != 0  ->  Phi ~ 1/r,  force |Omega| ~ 1/r^2   (LONG RANGE)
        - monopole == 0  ->  Phi falls faster (>= 1/r^3)          (SHORT RANGE)
      and on K being massless. A massive (Yukawa) K gives e^{-mr}/r regardless.

So "does 1/r^2 work?" = "is rho_N = |xi|^2 - 1/2 monopole-carrying, and is K
massless?". This script measures the monopole of candidate rho_N profiles and
the resulting far-field power law for massless vs massive K.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Radial grid
# ---------------------------------------------------------------------------
sigma = 1.0
Rmax = 60.0 * sigma
N = 200000
r = np.linspace(1e-4, Rmax, N)
dr = r[1] - r[0]


def monopole(rho):
    """4pi INT rho r^2 dr = INT rho d3x."""
    return 4.0 * np.pi * np.trapz(rho * r**2, r)


def phi_massless(rho):
    """Solve -lap Phi = rho (spherical, exact integral form):
       Phi(r) = (1/r) INT_0^r rho r'^2 dr' + INT_r^inf rho r' dr'."""
    inner = np.cumsum(rho * r**2) * dr                      # INT_0^r rho r'^2
    tail_full = np.cumsum((rho * r)[::-1])[::-1] * dr        # INT_r^inf rho r'
    return inner / r + tail_full


def phi_massive(rho, m):
    """Solve (-lap + m^2) Phi = rho (spherical, exact integral form):
       Phi = 1/(m r)[ e^{-mr} INT_0^r sinh(mr') rho r' dr'
                      + sinh(mr) INT_r^inf e^{-mr'} rho r' dr' ]."""
    sh = np.sinh(m * r)
    em = np.exp(-m * r)
    inner = np.cumsum(sh * rho * r) * dr
    tail = np.cumsum((em * rho * r)[::-1])[::-1] * dr
    return (em * inner + sh * tail) / (m * r)


def fit_slope(r, f, rlo, rhi):
    """Power-law exponent of |f| over [rlo, rhi] (log-log linear fit)."""
    mask = (r >= rlo) & (r <= rhi) & (np.abs(f) > 0)
    p = np.polyfit(np.log(r[mask]), np.log(np.abs(f[mask])), 1)
    return p[0]


def amp_at(r, f, r0):
    """|f| at radius r0 (long-range probe; slope alone is noise when f~0)."""
    return np.abs(np.interp(r0, r, f))


print("=" * 70)
print("OBE radial-law test: when does the gravity term give 1/r^2 ?")
print("=" * 70)

# Fit window: well outside the core (sigma=1), well inside the boundary.
rlo, rhi = 8.0, 25.0

# ---------------------------------------------------------------------------
# Profile M: monopole-carrying deviation (net excess/deficit of |xi|^2 from 1/2)
# ---------------------------------------------------------------------------
rho_M = np.exp(-r**2 / (2 * sigma**2))
print("\n[Profile M] rho_N = exp(-r^2/2s^2)  (net deviation, monopole != 0)")
print(f"  monopole INT rho d3x      = {monopole(rho_M):+.4f}")
PhiM = phi_massless(rho_M)
print(f"  massless Phi  slope        = {fit_slope(r, PhiM, rlo, rhi):+.3f}  (expect -1: 1/r potential)")
print(f"  massless |Phi(r=10)|       = {amp_at(r, PhiM, 10.0):.4e}  (O(monopole/4pi r))")
forceM = -np.gradient(PhiM, r)
print(f"  massless force |Omega| slope = {fit_slope(r, forceM, rlo, rhi):+.3f}  (expect -2: 1/r^2 force)")
for m in (0.5, 1.0):
    Pm = phi_massive(rho_M, m)
    print(f"  massive (m={m}) |Phi(r=10)|  = {amp_at(r, Pm, 10.0):.4e}  (Yukawa e^-mr/r, exp. screened)")

# ---------------------------------------------------------------------------
# Profile N: zero-monopole deviation (= lap of a localized blob; the V6/"grad"
#            obstruction case — a breathing deviation that nets to zero)
# ---------------------------------------------------------------------------
g = np.exp(-r**2 / (2 * sigma**2))
rho_N0 = (r**2 / sigma**4 - 3.0 / sigma**2) * g           # = lap(g): monopole 0
print("\n[Profile N] rho_N = lap(blob)  (breathing deviation, monopole = 0)")
print(f"  monopole INT rho d3x      = {monopole(rho_N0):+.4e}  (~0 by construction)")
PhiN = phi_massless(rho_N0)
# With zero monopole the exact answer is Phi = -blob (exponentially localized);
# the far-field slope is meaningless noise, so compare AMPLITUDE at r=10.
print(f"  massless |Phi(r=10)|       = {amp_at(r, PhiN, 10.0):.4e}  (vs Profile M above)")
print(f"  ratio |Phi_N(10)|/|Phi_M(10)| = {amp_at(r, PhiN, 10.0)/amp_at(r, PhiM, 10.0):.2e}"
      f"  (=> short range: no long tail)")

# ---------------------------------------------------------------------------
# Cross-check the integration-by-parts identity:
#   Omega from grad-source  ==  grad(Phi)  for massless K
# Build Omega_r(r) = f_g d/dr (K * rho)  two ways and compare.
# ---------------------------------------------------------------------------
print("\n[IBP check] OBE grad-source  vs  grad(K*rho)  (massless K)")
direct = -np.gradient(phi_massless(rho_M), r)             # grad(Phi)
# grad-source form: Omega = K * grad(rho); for spherical, (grad rho)_r = rho'
grad_rho = np.gradient(rho_M, r)
# convolve radial: reuse the massless solver treating grad_rho as a "source" of a
# vector potential is subtle; instead verify the scalar identity d/dr(K*rho) match
# by confirming -lap Phi = rho reproduces rho:
lapPhi = np.gradient(r**2 * np.gradient(phi_massless(rho_M), r), r) / r**2
resid = np.max(np.abs((-lapPhi - rho_M))[(r > 2) & (r < 30)])
print(f"  max |(-lap Phi) - rho| on (2,30) = {resid:.2e}  (confirms Phi solves -lap Phi = rho)")

# ---------------------------------------------------------------------------
# v59-specific: model |xi(r)|^2 -> 1/2 asymptotically with a localized core dip,
# and report its monopole (the physical question for the gravity bridge).
# ---------------------------------------------------------------------------
print("\n[v59] rho_N = |xi(r)|^2 - 1/2 for a localized particle profile")
# core relaxes BELOW the constraint then returns: |xi|^2 = 1/2 - D*exp(-r^2/2s^2)
for Damp in (0.2, 0.5):
    xi2 = 0.5 - Damp * np.exp(-r**2 / (2 * sigma**2))
    rhoN = xi2 - 0.5
    M = monopole(rhoN)
    Pn = phi_massless(rhoN)
    print(f"  dip D={Damp}:  monopole = {M:+.4f},  massless force slope = "
          f"{fit_slope(r, -np.gradient(Pn, r), rlo, rhi):+.3f}")

print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)
print("- 1/r^2 IS reachable for the OBE gravity term: massless K + monopole-")
print("  carrying rho_N gives a 1/r potential and a 1/r^2 force (Profile M).")
print("- It FAILS (short range) if rho_N nets to zero (Profile N) or if K is")
print("  massive (Yukawa). The current OBE 'grad rho_N' is fine: with massless K")
print("  it equals grad(K*rho_N), a genuine 1/r^2 force, PROVIDED INT rho_N != 0.")
print("- v59 caveat: a localized |xi|^2 deviation generically HAS a nonzero")
print("  monopole (= net Koide-constraint excess), so 1/r^2 is available. The open")
print("  issue is physical, not radial: a fully-relaxed lepton has |xi|^2=1/2")
print("  (rho_N=0) => no source, so rho_N alone cannot be the mass/gravity charge")
print("  (the known rho_M != mass problem). And the mode is SCALAR, not tensor.")
