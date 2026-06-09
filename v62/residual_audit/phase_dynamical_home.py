#!/usr/bin/env python3
"""
v62/residual_audit/phase_dynamical_home.py

THE ONLY VIABLE CATEGORY FOR THE BRANNEN PHASE.  (See ../THESIS.md, ../no_go/NOGO.md.)

The no-go shows the phase invariant cos(3 phi) = cos(2/3) is transcendental, so
NO algebraic construction (Casimir / character / J-eigenphase) can produce it.
The ONLY category of object that CAN equal a transcendental is a ratio of
dynamical (loop) coefficients.  This script makes that placement concrete and
states the honest open status -- it does NOT claim to derive phi = 2/9.

Established here:
  (1) The phi-potential is a Fourier series in 3 phi:  V_eff(phi) = c0 +
      c3 cos(3 phi) + c6 cos(6 phi) + ...  (a Z3-invariant function), so its
      critical points obey  cos(3 phi) = -c3 / (4 c6)  -- a RATIO of loop
      coefficients.  We extract c3, c6 from a self-contained Coleman-Weinberg
      potential and confirm they are nonzero.
  (2) Loop coefficients carry logarithms, so the ratio is a CONTINUOUS,
      transcendental function of the loop data (scale mu, counting N_X).  To
      reproduce the observed phi = 2/9 the ratio must equal 4 cos(2/3) =
      3.14355... -- a transcendental value, reachable by a loop ratio and
      UNreachable by any algebraic eigenphase.

Honest open status (unchanged from v59 Step-6): the simple 1-loop CW at the
physical lepton point is non-perturbative (one amplitude ~ vanishes -- the real
electron/tau hierarchy) and its leading term gives only the COMMENSURATE minima
sin(3 phi)=0.  Producing the small nonzero phi = 2/9 needs the full loop
structure (higher orders + the sector counting N_X derived from the G2-equivariant
projection).  That calculation is NOT done.  The DEFINITIVE result of v62 is only
the category placement: the phase is a P3 (dynamical) object, provably not a
P1/P2 (representation-theory) one.
"""
import numpy as np
import math

SEP = "=" * 78
COS_2_3 = math.cos(2.0 / 3.0)          # transcendental target ~ 0.785887
REQUIRED_RATIO = 4 * COS_2_3           # the coefficient ratio that yields phi=2/9


def banner(s):
    print(SEP); print(s); print(SEP)


def brannen_amps(phi, t, a=1.0):
    """Brannen amplitudes s_k = a(1 + 2 t cos(2 pi k/3 + phi))  (sqrt-masses)."""
    return np.array([a * (1 + 2 * t * math.cos(phi + 2 * math.pi * k / 3))
                     for k in range(3)])


def V_eff(phi, t, mu_sq, a=1.0, floor=1e-9):
    """Self-contained 1-loop Coleman-Weinberg potential for the Z3 kernel:
       V = sum_k m_k^2 (log(m_k^2/mu^2) - 3/2),  m_k = s_k^2  (lepton masses).
       abs+floor regularize the isolated zero-crossings (electron ~ 0)."""
    s = brannen_amps(phi, t, a)
    m = s**2
    return float(np.sum(m * (np.log(np.maximum(m, floor) / mu_sq) - 1.5)))


def fourier_cos(V_vals, n, phi_grid):
    dphi = phi_grid[1] - phi_grid[0]
    return (1.0 / math.pi) * float(np.sum(V_vals * np.cos(n * phi_grid)) * dphi)


def coeffs(t, mu_sq):
    g = np.linspace(0, 2 * math.pi, 6000, endpoint=False)
    V = np.array([V_eff(p, t, mu_sq) for p in g])
    c0 = (1.0 / (2 * math.pi)) * float(np.sum(V) * (g[1] - g[0]))
    return c0, fourier_cos(V, 3, g), fourier_cos(V, 6, g)


# ---------------------------------------------------------------------------
# Part 1 -- the phi-potential is a cos(3 phi) Fourier series
# ---------------------------------------------------------------------------
def part1_fourier_structure():
    banner("Part 1: V_eff(phi) is a Fourier series in 3 phi (perturbative t)")
    t = 0.3   # perturbative amplitude where the CW expansion is clean
    c0, c3, c6 = coeffs(t, 1.0)
    print(f"  t = {t} (perturbative), mu^2 = 1")
    print(f"    V(phi) = {c0:+.5f} + ({c3:+.5f}) cos(3 phi) + ({c6:+.5f}) cos(6 phi) + ...")
    print(f"    c3 nonzero: {abs(c3) > 1e-6};  c6 nonzero: {abs(c6) > 1e-8}")
    print("""
  V_eff is Z3-invariant, hence a function of cos(3 phi).  Its critical points
  (dV/dphi = 0) solve  sin(3 phi) (3 c3 + 12 c6 cos(3 phi) + ...) = 0, i.e.
        cos(3 phi) = -c3 / (4 c6)   [the nontrivial branch].
  The phase is therefore a RATIO of loop coefficients -- not an algebraic
  invariant.
""")
    return {"t": t, "c3": c3, "c6": c6}


# ---------------------------------------------------------------------------
# Part 2 -- the ratio is a continuous transcendental dial
# ---------------------------------------------------------------------------
def part2_transcendental_dial():
    banner("Part 2: the coefficient ratio is a continuous transcendental function")
    # NOTE: in this minimal CW potential the scale mu multiplies Sum m_k (which
    # is phi-flat), so mu only shifts c0 -- it does NOT dial the phase here.  The
    # genuine continuous input is the amplitude t (per-sector loop data, fixed
    # structurally by (1-t^2) D = 14, used here as a stand-in dial).
    print("  -c3/(4 c6) as the loop amplitude t varies (mu^2 = 1):")
    ts = np.linspace(0.15, 0.65, 11)
    vals = []
    for t in ts:
        _, c3, c6 = coeffs(t, 1.0)
        r = -c3 / (4 * c6)
        vals.append(r)
        print(f"    t = {t:4.2f}:  c3 = {c3:+.5f}, c6 = {c6:+.6f},  -c3/(4 c6) = {r:+.4f}")
    vals = np.array(vals)
    print(f"\n  The ratio is a CONTINUOUS, transcendental-valued function of the loop")
    print(f"  data (it sweeps {vals.min():.2f} .. {vals.max():.2f}); it is not a frozen")
    print(f"  algebraic number.")
    print(f"\n  At leading 1-loop order |c3| >> 4|c6|, so -c3/(4 c6) is far outside")
    print(f"  [-1, 1]: the minimal CW gives only the COMMENSURATE minima sin(3 phi)=0")
    print(f"  (exactly v59's finding).  An in-range nonzero phi needs the further")
    print(f"  sector-tilt term (N_X * alpha) -- the unfinished part.")
    print(f"\n  To reproduce the observed phi = 2/9 the relevant ratio must equal")
    print(f"    4 cos(2/3) = {REQUIRED_RATIO:.6f}  (a TRANSCENDENTAL value): reachable")
    print(f"  by a loop ratio, UNreachable by any algebraic eigenphase.")
    return {"ratio_range": [float(vals.min()), float(vals.max())],
            "required_ratio": REQUIRED_RATIO}


# ---------------------------------------------------------------------------
# Part 3 -- honest status + the category contrast that IS definitive
# ---------------------------------------------------------------------------
def part3_status_and_contrast():
    banner("Part 3: honest status + the definitive category contrast")
    print("""  OPEN (not done here, = v59 Step-6): the explicit loop calculation that fixes
  mu and derives the sector counting N_X, so that -c3/(4 c6) = cos(2/3) comes
  out rather than being matched.  The naive 1-loop CW at the physical point is
  non-perturbative (electron amplitude ~ 0) and gives only commensurate minima.

  DEFINITIVE (the v62 contribution):

    number-type      producible by                 can equal cos(2/3)?
    -------------    --------------------------    -------------------
    algebraic        P1 restriction, P2 Hessian,   NO   <- the no-go
                     character, J-eigenphase            (transcendental
                                                         target)
    transcendental   P3 loop-coefficient ratio     YES  (ratio of CW
                     -c3/(4 c6)  (mu, N_X dial)          coefficients is
                                                         transcendental)

  So: stop trying to make the algebra emit the phase (barred); the only path
  that can work is the loop Lagrangian.  Effort allocation follows from the
  number-type, not from taste.
""")


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

    t = 0.3
    _, c3, c6 = coeffs(t, 1.0)
    check("c3 (cos 3phi coefficient) nonzero", abs(c3) > 1e-6)
    check("c6 (cos 6phi coefficient) nonzero", abs(c6) > 1e-8)
    # the ratio is a continuous function of the loop data: it changes with t.
    _, c3a, c6a = coeffs(0.25, 1.0)
    _, c3b, c6b = coeffs(0.55, 1.0)
    ra, rb = -c3a / (4 * c6a), -c3b / (4 * c6b)
    check("ratio -c3/(4 c6) varies continuously with loop amplitude t",
          abs(ra - rb) > 1e-3)
    # the matching value is the transcendental 4 cos(2/3), NOT a simple rational.
    # (it differs from every low-height rational p/q with q<=50 by > 1e-4).
    best = min(abs(REQUIRED_RATIO - p / q)
               for q in range(1, 51) for p in range(0, 5 * q + 1))
    check("required ratio 4cos(2/3) is not a low-height rational", best > 1e-4)
    # masses finite (regularized) at the physical point.
    check("regularized V_eff finite at phi=2/9, t^2=1/2",
          math.isfinite(V_eff(2 / 9, math.sqrt(0.5), 1.0)))

    print(SEP)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    part1_fourier_structure()
    print()
    part2_transcendental_dial()
    print()
    part3_status_and_contrast()
    print()
    ok = run_checks()
    import sys
    sys.exit(0 if ok else 1)
