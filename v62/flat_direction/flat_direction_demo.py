#!/usr/bin/env python3
"""
v62/flat_direction/flat_direction_demo.py

DEMONSTRATION of the FLAT-DIRECTION LAW (see ../THESIS.md, ./FLAT_DIRECTION.md):

  Representation theory pins the symmetric skeleton of a sector and is
  constitutionally BLIND to the flat (angular) directions it leaves behind.
  Every residual conjecture of v59-v61 sits on such a flat direction:

    * the Brannen PHASE phi      -- a flat direction of every symmetric
                                    (Z3-invariant) mass invariant;
    * the EW DEMOCRACY (R1)      -- a flat direction (S^783) of the
                                    O(784)-symmetric Frobenius hat.

  Same structure both times: a symmetry-invariant functional fixes a
  radial/magnitude part and leaves an angular part flat.  The flat part is
  exactly what the algebra cannot fix.

  Resolution differs by NUMBER-TYPE of the flat coordinate:
    - EW democracy: a DEGENERACY (which point on the sphere) -> closable by
      state-selection [../residual_audit/democracy_selection.py];
    - Brannen phase: a coordinate whose VALUE is transcendental (cos(2/3))
      -> only a dynamical loop-ratio can lift it [../residual_audit/phase_dynamical_home.py].

This script proves both flatness statements: symbolically (Brannen) and
numerically (EW Hessian Goldstone count), with embedded self-checks.
"""
import sympy as sp
import numpy as np
import math

SEP = "=" * 78


def banner(s):
    print(SEP); print(s); print(SEP)


# ---------------------------------------------------------------------------
# Part A -- the Brannen phase is a flat direction of every symmetric invariant
# ---------------------------------------------------------------------------
def part_A_brannen_phase_flatness():
    banner("Part A: phi is a flat direction of all Z3-symmetric mass invariants")
    a, t, phi = sp.symbols('a t phi', positive=True, real=True)

    # Z3-symmetric Brannen amplitudes (the symmetric ansatz the algebra forces):
    s = [a * (1 + 2 * t * sp.cos(2 * sp.pi * k / 3 + phi)) for k in range(3)]
    m = [sp.expand_trig(sp.simplify(sk**2)) for sk in s]

    # Symmetric (elementary/power-sum) invariants:
    P1_s = sp.simplify(sum(s))                       # sum s_k
    e2_s = sp.simplify(s[0]*s[1] + s[1]*s[2] + s[2]*s[0])
    P1_m = sp.simplify(sum(m))                       # sum m_k
    # Koide ratio Q = (sum m) / (sum sqrt m)^2  = P1_m / P1_s^2
    Q = sp.simplify(P1_m / P1_s**2)

    # phi-DEPENDENT invariants (carry cos(3 phi)):
    e3_s = sp.simplify(s[0] * s[1] * s[2])           # product
    P3_s = sp.simplify(sum(sk**3 for sk in s))       # sum s_k^3

    print("  Symmetric invariants and their phi-derivative:")
    results = {}
    for name, expr in [("P1(s)=Sum s_k", P1_s), ("e2(s)", e2_s),
                       ("P1(m)=Sum m_k", P1_m), ("Koide Q", Q)]:
        d = sp.simplify(sp.diff(expr, phi))
        flat = (d == 0)
        results[name] = flat
        print(f"    {name:18s} = {sp.nsimplify(sp.trigsimp(expr))}")
        print(f"        d/dphi = {d}   -> {'FLAT (phi-independent)' if flat else 'depends on phi'}")

    print("\n  phi-DEPENDENT invariants (the only ones that see phi):")
    for name, expr in [("e3(s)=prod s_k", e3_s), ("P3(s)=Sum s_k^3", P3_s)]:
        d = sp.simplify(sp.diff(expr, phi))
        # Show the dependence is through sin(3 phi):
        carries_3phi = sp.simplify(d) != 0
        # confirm proportional to sin(3 phi)
        ratio = sp.simplify(d / sp.sin(3 * phi)) if carries_3phi else sp.Integer(0)
        results[name] = not carries_3phi
        print(f"    {name:18s}: d/dphi = {sp.trigsimp(d)}")
        print(f"        d/dphi / sin(3 phi) = {sp.trigsimp(ratio)}  (phi-dependence is via cos(3 phi))")

    print("""
  LAW (Brannen instance): every SYMMETRIC invariant -- the ones representation
  theory computes (dimensions, Koide Q, e2) -- is phi-INDEPENDENT.  phi enters
  ONLY the cos(3 phi)-class invariants.  So the symmetric structure fixes the
  whole skeleton and leaves phi as a flat direction it cannot see.
""")
    return results


# ---------------------------------------------------------------------------
# Part B -- the EW democracy is a flat direction (Goldstones) of the O(n) hat
# ---------------------------------------------------------------------------
def part_B_ew_democracy_flatness(n=8):
    banner(f"Part B: the O(n) Frobenius hat leaves a flat S^(n-1) (demo n={n}; real n=784)")
    # Frobenius Mexican hat: V(Y) = (lam/4)(||Y||_F^2 - v0)^2 depends ONLY on the
    # radial coordinate r = ||Y||_F.  Any vacuum with ||Y||_F = sqrt(v0) is a minimum;
    # the vacuum manifold is the sphere S^(n-1).  The Hessian at a vacuum has exactly
    # ONE positive (radial/Higgs) mode and n-1 ZERO modes (Goldstones).
    lam, v0 = 1.0, 1.0
    rng = np.random.default_rng(0)

    def V(Y):
        return 0.25 * lam * (Y @ Y - v0) ** 2

    def grad(Y):
        return lam * (Y @ Y - v0) * Y

    def hess(Y):
        # d2V = lam[(||Y||^2 - v0) I + 2 Y Y^T]
        return lam * ((Y @ Y - v0) * np.eye(n) + 2 * np.outer(Y, Y))

    # A vacuum: any unit-norm * sqrt(v0).  Test the DEMOCRATIC one and a random one.
    Y_dem = np.ones(n) / math.sqrt(n) * math.sqrt(v0)      # all components equal
    Y_rand = rng.standard_normal(n); Y_rand *= math.sqrt(v0) / np.linalg.norm(Y_rand)

    for label, Y in [("democratic (all equal)", Y_dem), ("random point on sphere", Y_rand)]:
        g = grad(Y)
        H = hess(Y)
        evals = np.sort(np.linalg.eigvalsh(H))
        nzero = int(np.sum(np.abs(evals) < 1e-9))
        npos = int(np.sum(evals > 1e-9))
        print(f"  {label}:")
        print(f"    ||Y||_F^2 = {Y@Y:.6f} (= v0),  |grad V| = {np.linalg.norm(g):.2e} (vacuum)")
        print(f"    Hessian spectrum: {npos} positive (Higgs), {nzero} zero (Goldstones)")
        print(f"    -> vacuum manifold dimension (flat) = {nzero} = n-1 = {n-1}")

    print(f"""
  LAW (EW instance): the O({n})-symmetric hat fixes only the RADIUS (the Higgs
  mode); the {n-1} angular directions are exactly degenerate (Goldstones).  At
  n=784 this is the S^783 of v61/EwVevHome.lean.  The DEMOCRATIC point is just
  ONE point on this flat sphere -- not selected by the symmetric potential.
  Same shape as the phase: magnitude fixed, angle flat.
""")
    return {"n": n, "goldstones": n - 1}


# ---------------------------------------------------------------------------
# Embedded self-checks
# ---------------------------------------------------------------------------
def run_checks():
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    a, t, phi = sp.symbols('a t phi', positive=True, real=True)
    s = [a * (1 + 2 * t * sp.cos(2 * sp.pi * k / 3 + phi)) for k in range(3)]
    P1_s = sp.simplify(sum(s))
    P1_m = sp.simplify(sum(sk**2 for sk in s))
    e2_s = sp.simplify(s[0]*s[1] + s[1]*s[2] + s[2]*s[0])
    Q = sp.simplify(P1_m / P1_s**2)
    e3_s = sp.simplify(s[0]*s[1]*s[2])

    check("Koide Q is phi-independent (d/dphi = 0)", sp.simplify(sp.diff(Q, phi)) == 0)
    check("e2(s) is phi-independent", sp.simplify(sp.diff(e2_s, phi)) == 0)
    check("P1(s) is phi-independent", sp.simplify(sp.diff(P1_s, phi)) == 0)
    check("e3(s) DOES depend on phi (the flat coordinate is visible only here)",
          sp.simplify(sp.diff(e3_s, phi)) != 0)
    # Koide value at the lepton symmetric point t^2 = 1/2 is exactly 2/3:
    Q_lepton = sp.simplify(Q.subs(t, sp.sqrt(sp.Rational(1, 2))))
    check("Koide Q = 2/3 at t^2 = 1/2 (symmetric structure fixes it)",
          sp.simplify(Q_lepton - sp.Rational(2, 3)) == 0)

    # EW: Hessian Goldstone count = n-1 for the Frobenius hat.
    n = 6
    Y = np.ones(n) / math.sqrt(n)
    H = (n * 0 + 1) * ((Y @ Y - 1) * np.eye(n) + 2 * np.outer(Y, Y))  # lam=1, v0=1
    evals = np.sort(np.linalg.eigvalsh(H))
    nzero = int(np.sum(np.abs(evals) < 1e-9))
    check(f"O({n}) Frobenius hat has n-1={n-1} Goldstones at a vacuum", nzero == n - 1)

    print(SEP)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    part_A_brannen_phase_flatness()
    print()
    part_B_ew_democracy_flatness()
    print()
    ok = run_checks()
    import sys
    sys.exit(0 if ok else 1)
