#!/usr/bin/env python3
"""
v63/dh_localization.py  --  READING B: Duistermaat-Heckman / character localization
on the canonical internal G2 geometry.

The G2 7-rep character is  chi_7(X) = sum_{weights mu} e^{i <mu, X>}  (a DH/Weyl
localization: a finite sum over the torus weights = the 6 short roots + 0).  We ask
the gating question: does a CANONICAL (parameter-free) evaluation point give the
phase 2/3?

Result: NO, and it was already PROVEN in v59 (`PhaseExclusions.lean`):
  * At a CANONICAL torsion element (e.g. the triality Z3 point), every weight
    angle <mu,X> is a rational multiple of 2*pi (COMMENSURATE), so each
    e^{i<mu,X>} is a root of unity and chi_7 is an ALGEBRAIC number (here = 2).
  * To get a contribution e^{i*2/3} (cos(2/3) transcendental) you need
    <mu,X> = 2/3, i.e. 2/3 in 2*pi*Q -- but 2/3 is NOT a rational multiple of pi
    (`phase_not_pi_rational`/`koide_not_pi_rational`).  So 2/3 can only appear at
    a GENERIC, non-commensurate X -- a free continuous parameter, not canonical.

So both the canonical-measure route (here) and the v59 holonomy analysis give the
same dichotomy:  TRANSCENDENTAL phase  <=>  OFF-LATTICE (non-pi-rational) point
<=>  a free continuous parameter.  Geometry/characters hand you pi-rational
(commensurate) angles; 2/3 is a BARE radian ratio -- the signature of a DYNAMICAL
loop-ratio, not a geometric/algebraic angle.
"""
import math
import cmath
import numpy as np
from fractions import Fraction

SEP = "=" * 78
TWO_PI = 2 * math.pi

# G2 fundamental 7-rep weights: the zero weight + the 6 short roots (unit hexagon).
WEIGHTS = [
    (0.0, 0.0),
    (1.0, 0.0), (-1.0, 0.0),
    (0.5, math.sqrt(3) / 2), (-0.5, -math.sqrt(3) / 2),
    (0.5, -math.sqrt(3) / 2), (-0.5, math.sqrt(3) / 2),
]


def banner(s):
    print(SEP); print(s); print(SEP)


def character(X):
    """chi_7(X) = sum_mu e^{i <mu, X>}  (real, since weights come in +/- pairs)."""
    return sum(cmath.exp(1j * (mu[0] * X[0] + mu[1] * X[1])) for mu in WEIGHTS)


# ---------------------------------------------------------------------------
# B1 -- canonical (commensurate) evaluation gives an ALGEBRAIC character
# ---------------------------------------------------------------------------
def B1_canonical_is_algebraic():
    banner("B1: canonical (triality) evaluation -> commensurate angles -> algebraic")
    Xc = (TWO_PI / 3, 0.0)   # a canonical torsion element
    chi = character(Xc)
    print(f"  X_canonical = (2*pi/3, 0)")
    angles = [(mu[0] * Xc[0] + mu[1] * Xc[1]) for mu in WEIGHTS]
    fracs = [Fraction(a / TWO_PI).limit_denominator(12) for a in angles]
    print(f"  weight angles / (2*pi): {[str(f) for f in fracs]}  (all RATIONAL = commensurate)")
    print(f"  chi_7(X_canonical) = {chi.real:+.6f} + {chi.imag:+.6f}i  -> = 2 (ALGEBRAIC/rational)")
    # every angle is a rational multiple of 2*pi -> e^{i*6*angle} = 1 (6th roots of unity)
    all_roots_of_unity = all(abs(cmath.exp(1j * 6 * a) - 1) < 1e-9 for a in angles)
    print(f"  every e^(i*angle) is a root of unity (e^(i*6*angle)=1): {all_roots_of_unity}")
    print("""
  A canonical torsion element gives commensurate (rational-multiple-of-2*pi)
  angles, so chi_7 is a sum of roots of unity -- an algebraic number.  No
  transcendental (cos(2/3)) can appear here.
""")
    return {"chi": chi.real, "all_roots_of_unity": all_roots_of_unity}


# ---------------------------------------------------------------------------
# B2 -- the phase 2/3 needs a GENERIC, non-commensurate point
# ---------------------------------------------------------------------------
def B2_phase_needs_offlattice():
    banner("B2: the phase 2/3 lives only at a non-commensurate (off-lattice) point")
    Xp = (2.0 / 3.0, 0.0)    # chosen so <(1,0),X> = 2/3
    chi = character(Xp)
    print(f"  X_phase = (2/3, 0): <(1,0),X> = 2/3  -> contributes 2*cos(2/3) (transcendental)")
    print(f"  chi_7(X_phase) = {chi.real:+.6f}  = 1 + 2cos(2/3) + 4cos(1/3) (TRANSCENDENTAL)")
    # Is 2/3 commensurate (a rational multiple of 2*pi)?  Search low-height q.
    best = min((abs(2.0 / 3.0 - TWO_PI * (p / r)), p, r)
               for r in range(1, 51) for p in range(0, 3 * r + 1))
    print(f"\n  Is 2/3 a rational multiple of 2*pi?  closest 2*pi*(p/r), r<=50:")
    print(f"    min |2/3 - 2*pi*(p/r)| = {best[0]:.4e}  (p/r = {best[1]}/{best[2]})  -- bounded away from 0")
    print(f"  -> 2/3 is NOT commensurate (proven exactly: phase_not_pi_rational, v59).")
    print(f"  -> e^(i*2/3) is never a root of unity; X_phase=(2/3,0) is a GENERIC point")
    print(f"     chosen BY HAND -- not a canonical/torsion element.  No parameter-free")
    print(f"     structure puts you there.")
    return {"chi_phase": chi.real, "min_commensurate_gap": best[0]}


# ---------------------------------------------------------------------------
# B3 -- the character is a continuous dial; canonical points are isolated
# ---------------------------------------------------------------------------
def B3_continuous_dial():
    banner("B3: chi is a continuous dial; the phase sits at no special point")
    xs = np.linspace(0, TWO_PI, 13)
    print("  X=(x,0):   x        chi_7(x)")
    for x in xs[::2]:
        print(f"            {x:6.3f}   {character((x, 0.0)).real:+.5f}")
    # demonstrate continuity / dial: chi changes smoothly with x
    c1, c2 = character((0.6, 0)).real, character((0.7, 0)).real
    print(f"\n  chi varies continuously (chi(0.6)={c1:.4f} != chi(0.7)={c2:.4f}).")
    print("  2/3=0.667 is an unremarkable interior point -- no canonical structure")
    print("  selects it.  Reaching 2/3 = choosing a continuous parameter = a dial.")
    return {"dial": abs(c1 - c2) > 1e-3}


# ---------------------------------------------------------------------------
def verdict():
    banner("VERDICT (Reading B) + the unifying dichotomy")
    print("""  Reading B reproduces, concretely, what v59 already PROVED
  (PhaseExclusions.lean E1/E1'):

    * canonical (commensurate / torsion) evaluation of the G2 character gives
      pi-rational angles -> roots of unity -> ALGEBRAIC values (chi_7 = 2);
    * the transcendental phase 2/3 is NOT pi-rational, so it appears only at a
      GENERIC, off-lattice point = a free continuous parameter, unforced by the
      algebra ("a free Aharonov-Bohm phase", koide_not_pi_rational).

  UNIFYING DICHOTOMY (both readings, and the whole v62->v63 arc, converge here):

      TRANSCENDENTAL phase 2/3  <=>  OFF-LATTICE (non-pi-rational) point
                                <=>  a free CONTINUOUS parameter.

  Geometry, characters, holonomies, and the uniform/canonical measure all hand
  you pi-rational (commensurate) angles -> algebraic.  The phase 2/3 is a BARE
  radian ratio (NOT pi-times-rational): the signature of a DYNAMICAL loop-ratio
  (v62 P3), not a geometric/algebraic angle.  This is exactly v59's E1 reading:
  "the phase must be a RATIO (like Q/3) inserted as a radian value, not an angle
  read off a rotation."

  So: integrating according to the CANONICAL measure (Reading B) cannot give a
  parameter-free 2/3 -- not because integration fails, but because 2/3 is not a
  geometric angle.  The only home left is the continuous, dynamically-pinned
  loop-ratio (v62/residual_audit/phase_dynamical_home.py).
""")


def run_checks(b1, b2, b3):
    banner("EMBEDDED SELF-CHECKS (Reading B)")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # canonical character is the algebraic value 2.
    check("chi_7 at the canonical triality element = 2 (algebraic/rational)",
          abs(b1["chi"] - 2.0) < 1e-9)
    check("all canonical weight phases are roots of unity (commensurate)",
          b1["all_roots_of_unity"])
    # 2/3 is not commensurate with 2*pi (off-lattice).
    check("2/3 is not a low-height rational multiple of 2*pi (off-lattice)",
          b2["min_commensurate_gap"] > 1e-3)
    # the phase character value is transcendental-flavored (!= the rational 2).
    check("chi_7 at the phase point != the canonical algebraic value",
          abs(b2["chi_phase"] - 2.0) > 1e-3)
    # character is a continuous dial.
    check("chi_7 is a continuous dial in X", b3["dial"])

    print(SEP)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    b1 = B1_canonical_is_algebraic()
    print()
    b2 = B2_phase_needs_offlattice()
    print()
    b3 = B3_continuous_dial()
    print()
    verdict()
    print()
    ok = run_checks(b1, b2, b3)
    import sys
    sys.exit(0 if ok else 1)
