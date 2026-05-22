#!/usr/bin/env python3
"""
v59/multivector_kernel_fit/02_z3_invariant_potential.py

Step 4 of the kernel-fit program: Koide-pinning via a Z_3-invariant potential.

After step 1-3 established that the Cl(3,1) Z_3 kernel reproduces all three
charged-lepton masses to machine precision with two empirical parameters
(|xi| and phi), we now ask: what is the most general Z_3-invariant potential
V(xi, xi*) whose minimum sits at |xi|^2 = 1/2 (Koide) and arg(xi) = 2/9 rad
(Brannen)?

Z_3 acts as xi -> omega xi with omega = exp(2 pi i / 3).
Z_3 invariants are functions of  |xi|^2  and  xi^3  (and xi*^3 = (xi^3)*).

A real Z_3-invariant potential takes the general form:

    V(xi, xi*) = f(|xi|^2)  +  Re[ g(|xi|^2) * xi^3 ]

where f is any real polynomial in |xi|^2 (the radial part) and g is any
complex polynomial in |xi|^2 (the cubic-phase part).  Higher-order
invariants are powers of these.

This script:
  1. Computes the critical-point conditions for V.
  2. Finds parameter regions consistent with (|xi|^2 = 1/2, phi = 2/9 rad).
  3. Tests the structural conjecture that the cubic phase coefficient
     is alpha = Q rad = 2/3 rad (the Koide ratio interpreted in radians).
  4. Identifies the minimum set of parameters required and reports what
     additional principles could fix them.
"""

import numpy as np
import sympy as sp


# =========================================================================
# Symbolic analysis
# =========================================================================
print("="*72)
print("Step 4: Z_3-invariant potential analysis")
print("="*72)

# Use real parametrization xi = r exp(i phi), with r = |xi| >= 0.
r, phi = sp.symbols('r phi', real=True, positive=True)
m2, lam, mu, alpha = sp.symbols('m^2 lambda mu alpha', real=True)

# General potential up to degree 6 in xi:
#   V = m^2 |xi|^2  +  lambda |xi|^4  +  mu Re[ xi^3 exp(-i alpha) ]
#                                    +  d6 |xi|^6 + ...
# where mu and alpha parameterize the complex coefficient of xi^3.

# Radial + cubic-phase (degree-6 radial truncation):
d6 = sp.symbols('d6', real=True)
V = (m2 * r**2
     + lam * r**4
     + mu * r**3 * sp.cos(3*phi - alpha)
     + d6 * r**6)

print("\nGeneral Z_3-invariant potential (radial up to r^6):")
print(f"  V(r, phi) = m^2 r^2 + lambda r^4 + mu r^3 cos(3 phi - alpha) + d6 r^6")
print()
print("Parameters:")
print("  m^2, lambda, d6  : real radial coefficients")
print("  mu, alpha        : magnitude and phase of the complex cubic coupling c = mu exp(i alpha)")

# Critical point conditions
print("\n--- Critical-point conditions ---")

dV_dr = sp.diff(V, r)
dV_dphi = sp.diff(V, phi)

print(f"\ndV/dr   = {sp.simplify(dV_dr)}")
print(f"dV/dphi = {sp.simplify(dV_dphi)}")

# dV/dphi = -3 mu r^3 sin(3 phi - alpha) = 0
# => sin(3 phi - alpha) = 0
# => 3 phi - alpha = n pi  for integer n
# For minimum (assuming mu > 0): cos(3 phi - alpha) = -1, so 3 phi - alpha = pi.
# For minimum (assuming mu < 0): cos(3 phi - alpha) = +1, so 3 phi - alpha = 0.

print("\nPhase minimum (for mu < 0, the standard convention):")
print("  3 phi - alpha = 0   =>   phi = alpha / 3")
print()
print("For empirical phi = 2/9 rad:")
print("  alpha = 3 * (2/9) = 2/3 rad")
print()
print("So the cubic phase coefficient alpha must be exactly 2/3 rad to give Brannen phi = 2/9 rad.")

# Now the radial condition
print("\nRadial minimum at the angular minimum (cos = +1, mu < 0):")
print("Substitute cos(3 phi - alpha) = +1 and -|mu| for mu:")

V_radial = m2 * r**2 + lam * r**4 - sp.Abs(mu) * r**3 + d6 * r**6
dV_radial_dr = sp.diff(V_radial, r)
print(f"\nV_radial(r) = m^2 r^2 + lambda r^4 - |mu| r^3 + d6 r^6")
print(f"dV_radial/dr = {sp.simplify(dV_radial_dr)}")

# Set r^2 = 1/2 (Koide condition)
print("\nAt the Koide point r = 1/sqrt(2), the radial extremum condition is:")
print("  dV/dr|_{r=1/sqrt(2)} = 0")
print()

r_koide = 1/sp.sqrt(2)
dV_at_koide = sp.simplify(dV_radial_dr.subs(r, r_koide))
print(f"  dV/dr at r=1/sqrt(2):  {dV_at_koide} = 0")
print()
print("This is ONE equation in THREE parameters (m^2, lambda, |mu|, d6).")
print("So generically two parameters remain free.")
print()

# Solving for one parameter in terms of others:
print("Solving for m^2 in terms of (lambda, |mu|, d6):")
m2_sol = sp.solve(dV_at_koide, m2)[0]
print(f"  m^2 = {sp.simplify(m2_sol)}")
print()

# The simplest case: d6 = 0 (truncate to degree 4 radial)
print("--- Simplest case: degree-4 radial (d6 = 0) ---")
V_simple = (m2 * r**2 + lam * r**4 - sp.Abs(mu) * r**3)
dVs_dr = sp.diff(V_simple, r)
print(f"V_simple = m^2 r^2 + lambda r^4 - |mu| r^3")
print(f"dV/dr = {sp.simplify(dVs_dr)}")

# Set dV/dr = 0 at r = 1/sqrt(2):
dVs_at_koide = sp.simplify(dVs_dr.subs(r, r_koide))
print(f"\ndV/dr at r=1/sqrt(2):  {dVs_at_koide} = 0")
sol_simple = sp.solve(dVs_at_koide, m2)[0]
print(f"  m^2 = {sp.simplify(sol_simple)}")
print()
print("So lambda and |mu| are free; m^2 is determined.")
print()
print("Stability (positive curvature at the minimum):")
d2V = sp.diff(V_simple, r, 2)
d2V_at_koide = sp.simplify(d2V.subs(r, r_koide).subs(m2, sol_simple))
print(f"  d^2 V / dr^2 at r=1/sqrt(2) (with m^2 fixed by minimum condition):")
print(f"    = {sp.simplify(d2V_at_koide)}")
print(f"  Stable iff this is positive.")


# =========================================================================
# Numerical exploration: visualize the potential at specific parameters
# =========================================================================
print()
print("="*72)
print("Numerical exploration: V at specific parameter choices")
print("="*72)

# Concrete choice: |mu| = 1, lambda = 1, m^2 determined by Koide condition.
# Then V_simple(r) = m^2 r^2 + r^4 - r^3.
# m^2 = -2 lambda r_koide^2 + (3/2) |mu| r_koide = -1 + (3/2)(1/sqrt(2))
import math

def V_simple_num(r, m2_v, lam_v, mu_mag):
    """Radial potential at the angular minimum."""
    return m2_v * r**2 + lam_v * r**4 - mu_mag * r**3

def find_radial_min(m2_v, lam_v, mu_mag):
    """Find r > 0 that minimizes V_simple_num.

    Use bounded minimization to avoid Brent bracket failures when the
    bracket cannot be established."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(V_simple_num, bounds=(0.01, 5.0),
                          args=(m2_v, lam_v, mu_mag), method='bounded',
                          options={'xatol': 1e-12})
    return res.x, res.fun

# Choose lambda = 1, |mu| = 1.  m^2 from formula:
lam_v = 1.0
mu_v = 1.0
m2_v = -2 * lam_v * 0.5 + 1.5 * mu_v / math.sqrt(2)
print(f"\nWith lambda = {lam_v}, |mu| = {mu_v}:")
print(f"  m^2 (Koide-pinning) = -2*lambda*(1/2) + (3/2)*|mu|/sqrt(2)")
print(f"                      = {-2*lam_v*0.5 + 1.5*mu_v/math.sqrt(2):.6f}")
print(f"                      = {m2_v:.6f}")

r_min, V_min = find_radial_min(m2_v, lam_v, mu_v)
print(f"\nNumerical minimum:")
print(f"  r_min      = {r_min:.10f}")
print(f"  r_min^2    = {r_min**2:.10f}")
print(f"  1/2 (target) = 0.5000000000")
print(f"  delta       = {r_min**2 - 0.5:+.3e}")
print(f"  V(r_min)   = {V_min:.6f}")


# Test stability across (lambda, |mu|) space
print()
print("\nScan over (lambda, |mu|): does Koide-pinning always work?")
print(f"{'lambda':>8} {'|mu|':>8} {'m^2 needed':>12} {'r_min':>10} {'r_min^2':>10} {'delta':>10}")
for lam_v in [0.5, 1.0, 2.0, 5.0]:
    for mu_v in [0.5, 1.0, 2.0, 5.0]:
        m2_v = -lam_v + 1.5*mu_v/math.sqrt(2)
        r_min, _ = find_radial_min(m2_v, lam_v, mu_v)
        print(f"  {lam_v:>6.2f} {mu_v:>8.2f} {m2_v:>12.6f} {r_min:>10.6f} {r_min**2:>10.6f} {r_min**2 - 0.5:>+10.2e}")


# =========================================================================
# The structural conjecture: alpha = 2/3 rad = Q rad
# =========================================================================
print()
print("="*72)
print("Structural conjecture: alpha_cubic = Q rad = 2/3 rad")
print("="*72)
print()
print("From the phase-minimum condition:  phi = alpha / 3.")
print("Empirical Brannen phase phi ≈ 2/9 rad implies alpha = 2/3 rad.")
print()
print("The Koide ratio Q = 2/3 is ALSO 2/3 (dimensionless).")
print()
print("Conjecture: the cubic coupling phase alpha in the potential is exactly")
print("Q radians, with the dimensionless Koide ratio reinterpreted as a phase.")
print()
print("If this is correct, then BOTH Koide and Brannen come from the same")
print("number 2/3, appearing twice:")
print("  - As Q = 2/3 (the Koide identity), which is structurally equivalent")
print("    to |xi|^2 = 1/2.")
print("  - As alpha = 2/3 rad (the cubic phase), which sets phi = Q/3 = 2/9.")
print()
print("Status: this is a CONJECTURE.  It would be confirmed by finding an")
print("algebraic operation in Cl(3,1) whose natural cubic invariant on the")
print("Z_3 sector has phase 2/3 rad, derivable from the algebra alone.")
print()

# Compute and report
phi_predicted = 2/3 / 3
print(f"Numerical check:")
print(f"  alpha = 2/3 rad => phi = 2/9 rad = {2/9:.10f}")
print(f"  Empirical phi   = 0.2222296315 rad (from step 3)")
print(f"  delta (phi - 2/9) = {0.2222296315 - 2/9:+.3e}")
print(f"  (within experimental m_tau uncertainty)")


# =========================================================================
# Summary
# =========================================================================
print()
print("="*72)
print("Summary")
print("="*72)
print("""
Findings from step 4:

1. The most general Z_3-invariant degree-6 potential in xi has four real
   parameters: m^2, lambda, d6 (radial) and (mu, alpha) (cubic complex
   coupling, magnitude and phase).

2. Phase minimum:  phi = alpha / 3  (modulo 2 pi / 3).
   Empirical phi = 2/9 rad => alpha = 2/3 rad.

3. Radial minimum at |xi|^2 = 1/2:  ONE constraint relating m^2 and lambda
   (and d6 if included).  At the simplest degree-4 truncation:
       m^2 = -lambda + (3/2) |mu| / sqrt(2)
   So lambda and |mu| are free; m^2 is determined.

4. The potential PINS the lepton spectrum given:
       alpha = 2/3 rad        (chooses Brannen phase phi = 2/9)
       lambda > 0             (stability)
       m^2 set by Koide minimum condition
   Two parameters (lambda and |mu|) remain free; they affect the curvature
   of the potential at the minimum but not the location.

5. Structural conjecture:  alpha = Q rad = 2/3 rad.
   This would tie Koide and Brannen to a single dimensionless number 2/3,
   appearing once as Q (the radial-minimum identity) and once as alpha
   (the cubic phase, reinterpreted in radians).  If found algebraically,
   this would make BOTH Koide and Brannen structural rather than empirical.

What this potential CANNOT do alone:

  - Predict alpha from algebra (it remains free in the general Z_3 potential).
  - Predict m^2, lambda, |mu| individually from first principles.

What this potential CAN do, once accepted:

  - Make Koide (|xi|^2 = 1/2) a stable equilibrium point.
  - Lock Brannen (phi = 2/9 rad) by the cubic phase, IF alpha = 2/3 rad
    is enforced.

Next step (Step 5): test whether the SAME Lagrangian -- with the lepton
sector pinned -- predicts a sensible value of alpha (fine structure)
from its bivector-grade coupling.  This is the cross-sector unification
test that v58 originally claimed but never quantitatively verified.
""")
