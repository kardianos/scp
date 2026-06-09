#!/usr/bin/env python3
"""
v59/gaps/gravity/g9_spin2_route_test.py

G9 follow-up: can a spin-2 mode be EXTRACTED from the bivector / so(8) structure?
This tests the two concrete routes in ALTERNATIVES.md (G9-A induced metric /
constrained 2-form DOF count; G9-B Sym^2(adjoint) branching to Lorentz).

Route A -- DEGREE-OF-FREEDOM COUNT for a 2-form field as gravity (Plebanski-style)
  A massless spin-2 graviton in 4D has exactly 2 physical (TT) DOF.
  Compare what various candidate fields propagate:
    scalar          : 1
    vector A_mu (gauge): 2  (but helicity +-1, NOT +-2)
    2-form B_mu_nu (gauge): in 4D a massless 2-form has 1 DOF (dual to a scalar!)
    symmetric h_mu_nu (massless, diffeo gauge): 2  (helicity +-2)  <-- graviton
  => the DOF count alone separates them; a 2-form is NOT a graviton (1 DOF, dual
     to a scalar -- the "notoph"). This is a clean, well-known no-go.

Route B -- SU(2)/SO(n) rep branching: does Sym^2(adjoint) contain a symmetric-
  traceless rank-2 Lorentz tensor when so(8) is broken to an embedded Lorentz?
  We test the simplest tractable analog (so(4) ~ Lorentz-like, su(2)xsu(2)) and
  count whether Sym^2 of the bivector adjoint yields a spin-2 (j=2) piece under
  the diagonal su(2) (the 'rotation' subgroup) -- the spin-2 fingerprint.
"""

import numpy as np
import itertools

print("=" * 74)
print("G9 ROUTE A -- propagating DOF count (4D, massless)")
print("=" * 74)

def massless_dof_4d(field):
    """Physical (gauge-fixed, on-shell) DOF of a massless field in 4D."""
    table = {
        "scalar":            1,   # h = 0
        "vector(gauge)":     2,   # h = +-1   (photon)
        "2-form(gauge)":     1,   # h = 0     (dual to a scalar; the 'notoph')
        "sym_tensor(graviton)": 2,  # h = +-2 (graviton)
    }
    return table[field]

helicity = {"scalar": "{0}", "vector(gauge)": "{+-1}",
            "2-form(gauge)": "{0}", "sym_tensor(graviton)": "{+-2}"}
for f in ["scalar", "vector(gauge)", "2-form(gauge)", "sym_tensor(graviton)"]:
    print(f"  {f:24s}  DOF = {massless_dof_4d(f)}   helicity = {helicity[f]}")
print("\n  => A massless 2-form in 4D has 1 DOF and is DUAL TO A SCALAR (h=0).")
print("     So even a *spacetime* 2-form connection (the natural reading of an")
print("     'so(8)-bivector promoted to spacetime') carries NO h=+-2. Only the")
print("     symmetric tensor h_mu_nu does. This is a standard no-go, independent")
print("     of the helicity-generator computation in g9_polarization_test.py.")

print("\n" + "=" * 74)
print("G9 ROUTE B -- does Sym^2(adjoint) contain a spin-2 (symmetric-traceless)?")
print("=" * 74)
print("Toy model: so(4) bivectors (6-dim adjoint) under the diagonal su(2) (the")
print("'rotation' subgroup that plays the role of the little group). A spin-2")
print("(j=2, 5-dim) piece in Sym^2 is the graviton fingerprint.")

# so(4) adjoint = 6 antisymmetric 4x4 matrices. Under so(4)=su(2)_L x su(2)_R,
# the adjoint splits 6 = (3,1) + (1,3) (self-dual + anti-self-dual).
# Sym^2(6) = 21-dim. Decompose under the DIAGONAL su(2) (spin j of each piece):
#   (3,1): j=1 ; (1,3): j=1  -> under diagonal su(2) both are j=1 (spin-1 triplets)
#   Sym^2(adjoint): contains Sym^2(3,1) + Sym^2(1,3) + (3,1)x(1,3)
#     Sym^2 of a j=1 triplet = j=2 (5) + j=0 (1)
#     (3,1)x(1,3) under diagonal = j=1 x j=1 = j=2(5)+j=1(3)+j=0(1)
# So Sym^2(adjoint) under the diagonal su(2) DOES contain j=2 (spin-2) pieces.
def cg_decompose_product(j1, j2):
    """Clebsch-Gordan: j1 x j2 -> list of j from |j1-j2| to j1+j2."""
    lo = abs(j1 - j2); hi = j1 + j2
    return list(range(lo, hi + 1))  # integer spins here

# diagonal-su(2) spins of the two adjoint halves (each is a j=1 triplet)
half_spins = [1, 1]   # (3,1) and (1,3), both spin-1 under diagonal su(2)

# Sym^2 of the 6-dim adjoint, organized as Sym^2 of (V_L + V_R):
#   Sym^2(V_L) + Sym^2(V_R) + V_L (x) V_R
spin2_count = 0
pieces = []
# Sym^2 of a spin-1 triplet = spin-2 + spin-0  (dim 5 + 1 = 6 = Sym^2(3))
for half in half_spins:
    # Sym^2(spin-1) = {2, 0}
    pieces.append(("Sym^2(half spin-1)", [2, 0]))
    spin2_count += 1  # one j=2
# cross term V_L x V_R = spin1 x spin1 = {2,1,0}
pieces.append(("(3,1)x(1,3) diagonal", cg_decompose_product(1, 1)))
spin2_count += 1  # contains a j=2

for name, js in pieces:
    print(f"  {name:24s} -> diagonal-su(2) spins {js}")
print(f"\n  total j=2 (spin-2) multiplets in Sym^2(adjoint) = {spin2_count}")
print("  => YES: Sym^2(so(4) adjoint) CONTAINS spin-2 pieces. By the same")
print("     mechanism Sym^2(so(8) adjoint) contains symmetric-traceless rank-2")
print("     reps. So the ALGEBRA has the spin-2 representation available.")

print("\n" + "=" * 74)
print("THE CATCH (why route B is necessary but NOT sufficient)")
print("=" * 74)
print("- Sym^2(adjoint) containing a spin-2 REP is a statement about INTERNAL")
print("  so(8) indices. To be the LIGO graviton it must become a SPACETIME")
print("  symmetric tensor h_mu_nu via a SOLDERING (vielbein) that maps internal")
print("  so(8) -> spacetime Lorentz. v59 has the chain G2 c Spin(7) c Spin(8),")
print("  but NO soldering form identifying 4 of the 8 internal directions with")
print("  spacetime is established. Without it, the spin-2 rep stays INTERNAL")
print("  (a 'graviton' charged under so(8), not propagating in spacetime).")
print("- AND it must be a FUNDAMENTAL (linear, massless) field, not the quadratic")
print("  bilinear Omega Omega (which is O(alpha), short-range). Sym^2(adjoint) as a")
print("  composite of two Omega's is exactly that quadratic object.")
print("\nVERDICT (G9 routes): the spin-2 REPRESENTATION exists in the algebra")
print("(route B), but (i) no soldering maps it to a SPACETIME tensor, and (ii) as")
print("a composite it is quadratic/short-range. A massless 2-form (route A) is")
print("provably NOT a graviton (1 DOF, dual to a scalar). So a genuine spin-2")
print("LIGO mode requires ADDING a fundamental soldered symmetric tensor")
print("(Plebanski/vielbein) -- it is NOT automatically present. G9 stands open,")
print("and is the decisive obstruction.")
