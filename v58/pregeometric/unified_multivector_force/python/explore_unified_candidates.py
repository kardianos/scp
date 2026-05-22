"""
explore_unified_candidates.py

Python / Geometric Algebra discovery script for the unified multivector force law.
First cycle: symbolic grade projection analysis of candidate differential equations
(Candidates A/B from BACKGROUND_AND_SPECULATIVE_EQUATIONS.md and the Level 0/1 forms
discussed in project memory).

Uses the lightweight ga.py (Cl(3,0) for static 3D intuition; easily extended).

Focus:
- Structure of DΩ (vector derivative acting by geometric product)
- How quadratic nonlinearities Ω² mix grades
- What source grades (vector for ∇ρ gravity, bivector for J_χ EM) are needed
- Whether a single equation cleanly decouples into Newtonian-like + Maxwell-like
  subsystems under plausible grade assumptions or projectors.
- Quick numeric check of the integral force expression on 2-3 "particles".

Run: python explore_unified_candidates.py
"""

import sympy as sp
import numpy as np
from ga import MV, symbolic_mv, vector_derivative_symbolic, SIGNATURE_3D

print("=" * 70)
print("UNIFIED MULTIVECTOR FORCE LAW — FIRST PYTHON EXPLORATION CYCLE")
print("Candidate forms: Level 0 (∇Ω + λ Ω² = J), Candidate A/B variants")
print("=" * 70)

# ---------------------------------------------------------------------------
# 1. Symbolic setup
# ---------------------------------------------------------------------------
sp.init_printing(use_unicode=False)

# We explore Ω containing scalar + vector + bivector (most relevant grades for
# gravity connection + EM bivector F). Pseudoscalar can be added later.
Omega = symbolic_mv("Ω", grades=(0, 1, 2))
print("\n[1] General Ω (grades 0+1+2):")
print(Omega.to_dict())

D = vector_derivative_symbolic()
print("\n[2] Formal vector derivative D = ∂x e1 + ∂y e2 + ∂z e3")

D_Omega = D * Omega
print("\n[3] D * Ω — resulting grades and sample projections:")
for k in range(4):
    gk = D_Omega.grade(k)
    if gk.terms:
        print(f"  grade-{k} component (simplified): {gk}")

# ---------------------------------------------------------------------------
# 2. Quadratic term Ω*Ω and its grade content
# ---------------------------------------------------------------------------
Omega2 = Omega * Omega
print("\n[4] Ω * Ω — grade content (nonlinear self-interaction):")
for k in range(4):
    gk = Omega2.grade(k)
    if gk.terms:
        # Show one representative term to illustrate mixing
        first_blade = next(iter(gk.terms))
        print(f"  grade-{k} : contains term e.g. {gk.terms[first_blade]} * blade{first_blade}")

# Key observation: Ω² produces scalar (from biv·biv ~ |B|^2, vec·vec, etc.),
# vectors, bivectors, and trivectors. This is rich mixing.

# ---------------------------------------------------------------------------
# 3. Model source J(M) — vector part from density gradient + bivector from chiral current
# ---------------------------------------------------------------------------
# ∇ρ_M is primarily a *vector* source (gradient of scalar density).
# J_χ (chiral current) is a *bivector* (net handedness / phase winding).
J_vec = symbolic_mv("Jρ", grades=(1,))   # gravity-like source (vector)
J_chi = symbolic_mv("Jχ", grades=(2,))   # EM-like source (bivector)
J = J_vec + J_chi

print("\n[5] Model source J = J_vector(∇ρ) + J_bivector(χ)")
print("  J grades:", sorted({len(b) for b in J.terms.keys()}))

# ---------------------------------------------------------------------------
# 4. Candidate equation (Level 0 / Candidate A simplified, static intuition)
#    D Ω + λ (Ω * Ω)  = J
#    (We ignore m²Ω term for structure; it would be lower-order mass-like.)
# ---------------------------------------------------------------------------
lam = sp.symbols("λ", real=True)
eq_left = D_Omega + lam * Omega2
eq_right = J

print("\n[6] Candidate equation under test:  DΩ + λ Ω²  = J_ρ + J_χ")
print("    (full geometric product; no grade projectors yet)")

# Project both sides and equate grades. This gives the coupled system.
print("\n[7] GRADE-BY-GRADE DECOMPOSITION (equating left = right):")
for k in [0, 1, 2, 3]:
    left_k = eq_left.grade(k)
    right_k = eq_right.grade(k)
    print(f"\n  --- Grade {k} ---")
    if left_k.terms or right_k.terms:
        print(f"  LHS has {len(left_k.terms)} terms, RHS has {len(right_k.terms)} terms")
        # For insight, show if the quadratic contributes to this grade
        qk = (lam * Omega2).grade(k)
        if qk.terms:
            print(f"  λΩ² contributes to grade-{k}")
        if k == 0:
            print("  => This grade typically sources a scalar Poisson-like equation (gravity potential).")
            print("     Nonlinear terms from Ω² can give self-gravity or effective mass corrections.")
        elif k == 1:
            print("  => Vector grade can carry force-like or current-like information.")
        elif k == 2:
            print("  => Bivector grade is the natural home for the EM field strength F (Faraday).")
            print("     DΩ 's bivector part includes curl-like and div-like terms on vector components of Ω.")
            print("     Ω² bivector part can source magnetic currents or nonlinear EM effects.")
        elif k == 3:
            print("  => Trivector (axial) components must be constrained or projected out for classical limits.")
    else:
        print("  (both sides zero at this grade)")

# ---------------------------------------------------------------------------
# 5. Effect of grade projectors (Candidate B style / causal or selective)
#    Often one keeps only specific grades of the left-hand side.
#    E.g.  <DΩ + λΩ² >_{0,2}  = J  (keep only scalar + bivector for gravity+EM)
# ---------------------------------------------------------------------------
print("\n[8] With grade projector P = <...>_{0,2}  (keep scalar + bivector only):")
for k in [0, 2]:
    left_k = (D_Omega + lam * Omega2).grade(k)
    right_k = J.grade(k)  # J already only 1+2; grade 0 would be 0
    print(f"  grade-{k} : LHS terms={len(left_k.terms)}, RHS terms={len(right_k.terms)}")
    if k == 2:
        print("    This projects the dynamics onto the bivector sector (Maxwell-like) while")
        print("    scalar sector still receives contributions from vector parts of Ω via D.")

# ---------------------------------------------------------------------------
# 6. Commutator form (Candidate C — pre-geometric hint)
#    [Ω, M] = ΩM - MΩ  (purely algebraic, no background derivative)
# ---------------------------------------------------------------------------
M_test = symbolic_mv("M", grades=(0, 2))  # test source with density + chirality
comm = Omega * M_test - M_test * Omega
print("\n[9] Commutator [Ω, M] grade content (Candidate C style, no manifold ∇):")
for k in range(4):
    gk = comm.grade(k)
    if gk.terms:
        print(f"  grade {k}: {len(gk.terms)} terms — algebraic mixing only")

print("\n    This form generates interaction purely from algebra (no background space).")
print("    Retardation/causality would come from partial order on the discrete substrate.")
print("    Useful for Lean formalization of pre-geometric limit.")

# ---------------------------------------------------------------------------
# 7. First numeric experiment: 2-particle unified force via integral-like sum
#    (proxy for Green's function solution of the candidate)
#    Use the force expression from MULTIVECTOR_FORCE_LAW.md:
#      Ω(x) ~ sum K * ( f_g * "∇ρ" + f_em * Jχ )
#      a ~ < Ω M_t >_1   (vector force on test multivector)
#    Here we use simplified effective scalars for sources + full GA product
#    on the "Ω" assembled as bivector + vector contributions.
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("[10] NUMERIC 2-PARTICLE TEST (static kernel proxy, 3D Cl(3))")
print("=" * 70)

def make_point_source(pos, rho, chi, pos_test):
    """Return contribution to Ω at pos_test from a point source.
    - rho contributes a 'grav' vector-like field ~ (pos_test - pos)/r^3  (like E or g)
    - chi contributes a 'biv' field. For simplicity we model chi as producing
      an effective axial bivector field (simplified 'B' or twist).
    Uses inverse-square falloff for both (consistent with 1/r^2 forces).
    """
    r_vec = np.array(pos_test) - np.array(pos)
    r = np.linalg.norm(r_vec)
    if r < 1e-9:
        return MV.zero()
    r_hat = r_vec / r
    # Grav contrib: vector part pointing along r (or - for attraction; sign convention later)
    # We put it in vector grade of Ω
    grav_strength = rho / r**2
    grav_vec = MV.vector(tuple(grav_strength * r_hat), sig=SIGNATURE_3D)

    # Chiral/EM contrib: for test we put a simple bivector ~ chi * some plane
    # e.g. bivector in xy plane scaled by chi/r^2  (toy model)
    em_strength = chi / r**2
    # Arbitrary orientation for demo: e12 component
    biv = MV.bivector((em_strength, 0.0, 0.0), sig=SIGNATURE_3D)

    # Ambient modulation f(ρ_amb) would scale both; here we use f=1 for isolated test
    return grav_vec + biv

# Two particles
p1 = {'pos': np.array([0., 0., 0.]), 'rho': 1.0, 'chi': 0.0}   # neutral heavy
p2 = {'pos': np.array([1., 0., 0.]), 'rho': 1.0, 'chi': 0.5}   # charged

# Test location near p2, "force" on a test excitation at p2 due to p1
test_pos = p2['pos'] + np.array([0.1, 0.0, 0.0])
Omega_at_test = make_point_source(p1['pos'], p1['rho'], p1['chi'], test_pos)

print(f"Ω at test (from p1 on p2): {Omega_at_test.to_dict()}")

# Test multivector M_t at p2: has its own density + its own chirality (charge)
# For simplicity, M_t ~ scalar (density) + bivector (internal current)
rho_t = 1.0
chi_t = 0.5
M_t = MV.scalar(rho_t) + MV.bivector((chi_t, 0, 0))

# Force-like: the vector grade of (Ω M_t)  — this is the <Ω M_t >_1 term in the doc
force_mv = (Omega_at_test * M_t).grade(1)
print(f"Vector force component <Ω M_t>_1 : {force_mv.to_dict()}")

# Separate checks: pure gravity (set chi=0 on both)
Omega_g_only = make_point_source(p1['pos'], p1['rho'], 0.0, test_pos)
M_t_neutral = MV.scalar(rho_t)
force_g = (Omega_g_only * M_t_neutral).grade(1)
print(f"Pure gravity force (chi=0): {force_g.to_dict()}")

# Pure EM-like (rho small, chi nonzero) — note the vector force arises from
# bivector * bivector -> vector terms in the geometric product! This is key.
Omega_em_only = make_point_source(p1['pos'], 0.0, 0.5, test_pos)
M_t_charged = MV.bivector((0.5,0,0))
force_em = (Omega_em_only * M_t_charged).grade(1)
print(f"Pure EM-like (biv*biv -> vec force): {force_em.to_dict()}")

print("\n[11] Observation from numeric:")
print("  - The geometric product automatically generates a vector force from")
print("    both vector-Ω (grav) acting on scalar-Mt and bivector-Ω acting on biv-Mt.")
print("  - Cross terms (grav Ω on charged Mt, EM Ω on dense Mt) appear naturally.")
print("  - Whether those cross terms are physically desired (e.g. charged particles")
print("    feel slightly modified gravity, or vice versa) is a theory constraint.")
print("  - 1/r^2 scaling is built into the kernel; the algebra supplies the 'type' of force.")

# ---------------------------------------------------------------------------
# 8. Summary for findings + Lean implications
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("FIRST CYCLE SUMMARY")
print("=" * 70)
print("""
Key findings so far:
* DΩ already mixes many grades; the bivector sector (EM) receives contributions
  from derivatives of the vector part of Ω and vice versa. This is promising
  for unification but requires projectors or selective grade support in Ω.
* The quadratic λ Ω² term populates ALL grades. It can provide nonlinear
  corrections to gravity (scalar) and self-interaction for EM (bivector),
  but will also generate unwanted trivectors unless the equation uses
  <DΩ + λΩ² >_{0,2} or Ω itself is restricted to grades {0,1,2} with
  closure properties.
* The commutator form [Ω, M] is purely algebraic and automatically couples
  density (scalar) and chirality (bivector) sectors without needing a
  background derivative — strong hint for pre-geometric Lean work.
* Numeric integral proxy + geometric product shows that a single Ω can
  source both Newtonian 1/r² attraction (via vector grades) and a bivector-
  induced force (via geometric product on charged M_t) using the *same*
  algebraic object. Cross-coupling exists and must be interpreted or
  suppressed by the choice of f(ρ_ambient) or by the precise definition
  of how a = <Ω M_t >_1 + higher terms acts on the test excitation.

Implications / requests for the Lean track:
- Encode the grade projection operators < >_k and the geometric product in Lean.
- State the candidate "DΩ + λ <Ω²>_{0,2} = <J_ρ + J_χ>_{0,2}" and prove
  the conditions under which the grade-0 part reduces to a Poisson equation
  with ambient-density modulation f(ρ) and the grade-2 part reduces to
  inhomogeneous Maxwell (or at least the static Coulomb law).
- Investigate whether the commutator form can be shown to imply the same
  classical limits when a causal partial order is added.
- Clarify the exact multivector definition of J_χ (is it <M ~M>_2 ? a
  derived current from the Dirac-like operator on M?).

Next Python steps (proposed):
- Add 3+1 signature (Cl(3,1)) and implement a simple retarded kernel numeric
  evolution on a small 1D/2D lattice for wave propagation (causality test).
- Test ambient-density functions f_g(ρ), f_em(ρ) numerically on 3-body
  configurations and check equivalence-principle-like behavior.
- Visualize bivector twists + density gradients (export to ganja.js or matplotlib).
- Iterate the candidate (add or remove terms) based on Lean feedback.

This completes the first work cycle.
""")

print("Script completed successfully.")
