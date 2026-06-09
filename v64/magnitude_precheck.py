#!/usr/bin/env python3
"""
v64 Guardrail §3 — Magnitude pre-check for the density->proper-time (lapse) mechanism.

Question: does the lapse  alpha(P) = 1/sqrt(1 + kappa*P^n)  produce a *gravitational*
time-dilation, or the same nuclear-scale, too-strong, short-range effect every prior
BLV / effective-metric computation in this project found?

This is analytic (no simulation). It computes, for a plausible range of soliton core
densities P:
  - the lapse alpha and the effective potential  Phi/c^2 = ln(alpha)
  - the implied force range (Yukawa, set by the soliton core scale)
  - the ratio of Phi/c^2 to the *gravitational* Phi/c^2 at the proton scale
and reports which side of the "nuclear wall" the mechanism sits on.

Robust to the kP vs kP^2 ambiguity (CONCEPT posited 1+kP; kernel dynamics carry 1+kP^2)
-- both forms are evaluated. Embedded self-checks at the bottom (exit 1 on failure).

Refs: v60/README.md (scalar g00 is LIGO-fatal), MEMORY 'BLV Effective Metric'
(Phi_min/c^2=-0.275, Yukawa ~0.8 fm, 10^30-10^35x too strong),
CLAUDE.md 'Parameter Fitting' (code-unit -> physical conversions).
"""
import math

# ---------------------------------------------------------------------------
# Fixed theory / unit anchors
# ---------------------------------------------------------------------------
KAPPA = 50.0                      # CLAUDE.md standard nonlinearity

# Code-unit -> physical conversions (CLAUDE.md 'Parameter Fitting', e=1, rho0=1)
CODE_L_FM = 0.5624                # 1 code length in fm
CODE_E_MEV = 9.098               # 1 code energy in MeV

# Physical constants (SI)
G    = 6.674e-11                  # m^3 kg^-1 s^-2
C    = 2.998e8                    # m/s
HBAR = 1.0546e-34                # J s
M_PROTON = 1.6726e-27            # kg
R_PROTON_FM = 0.84               # charge radius, fm
FM = 1e-15                        # m

# ---------------------------------------------------------------------------
# The gravitational yardstick: Phi/c^2 = G m / (r c^2) at the proton scale
# ---------------------------------------------------------------------------
def newtonian_phi_over_c2(m_kg, r_m):
    return G * m_kg / (r_m * C * C)

PHI_GRAV_PROTON = newtonian_phi_over_c2(M_PROTON, R_PROTON_FM * FM)

# ---------------------------------------------------------------------------
# The lapse mechanism
# ---------------------------------------------------------------------------
def lapse(P, n):
    """alpha = 1/sqrt(1 + kappa P^n).  n=1 posited (CONCEPT), n=2 kernel dynamics."""
    return 1.0 / math.sqrt(1.0 + KAPPA * P**n)

def phi_over_c2(P, n):
    """Effective gravitational potential of the lapse: Phi/c^2 = ln(alpha)."""
    return math.log(lapse(P, n))

# ---------------------------------------------------------------------------
# Plausible core densities P = phi0*phi1*phi2 (local product at the soliton core)
# Bracket: field amplitudes O(0.5..2) at the core -> P ~ 0.1 .. 8.
# Also the integrated-density reading (P_int ~320, CLAUDE.md) as an upper bracket.
# ---------------------------------------------------------------------------
CORE_P_SCAN = [0.1, 0.5, 1.0, 2.0, 3.4, 5.0, 8.0]
P_INT_UPPER = 320.0              # integrated-density reading (loose upper bound)

# Yukawa range: set by the soliton core size ~ 1 code length (memory: ~0.5-0.8 fm)
RANGE_FM = CODE_L_FM             # ~0.56 fm half-width

def report():
    print("=" * 74)
    print("v64 MAGNITUDE PRE-CHECK — lapse density->time mechanism vs gravity")
    print("=" * 74)
    print(f"kappa = {KAPPA}   range (Yukawa half-width) ~ {RANGE_FM:.2f} fm")
    print(f"Gravitational yardstick: Phi/c^2 at proton scale "
          f"= G m_p/(r_p c^2) = {PHI_GRAV_PROTON:.3e}")
    print()
    header = f"{'P (core)':>9} | {'alpha (n=1)':>11} {'Phi/c2 n=1':>11} | " \
             f"{'alpha (n=2)':>11} {'Phi/c2 n=2':>11} | {'ratio/grav (n=2)':>16}"
    print(header)
    print("-" * len(header))
    rows = []
    for P in CORE_P_SCAN + [P_INT_UPPER]:
        a1, a2 = lapse(P, 1), lapse(P, 2)
        p1, p2 = phi_over_c2(P, 1), phi_over_c2(P, 2)
        ratio2 = abs(p2) / PHI_GRAV_PROTON
        tag = "  <-- P_int reading" if P == P_INT_UPPER else ""
        print(f"{P:>9.2f} | {a1:>11.4f} {p1:>11.3f} | "
              f"{a2:>11.4f} {p2:>11.3f} | {ratio2:>16.2e}{tag}")
        rows.append((P, a1, p1, a2, p2, ratio2))
    print()

    # Headline comparison at a representative core density P=1.0
    a2 = lapse(1.0, 2); p2 = phi_over_c2(1.0, 2)
    depth_mev = abs(p2) * CODE_E_MEV    # rough potential depth in MeV (Phi/c^2 * E_scale)
    print("Representative core (P=1.0, n=2, kernel form):")
    print(f"  alpha          = {a2:.4f}   (proper time runs {1/a2:.1f}x slower at core)")
    print(f"  Phi/c^2        = {p2:.3f}   (cf. gravity {PHI_GRAV_PROTON:.1e})")
    print(f"  ratio to grav  = {abs(p2)/PHI_GRAV_PROTON:.2e} x too strong")
    print(f"  depth          ~ {depth_mev:.0f} MeV   (nuclear scale, cf. BLV ~260 MeV)")
    print(f"  range          ~ {RANGE_FM:.2f} fm Yukawa (NOT 1/r; cf. memory ~0.8 fm)")
    print()

    # Verdict
    min_ratio = min(r[5] for r in rows)
    print("VERDICT")
    print("-" * 7)
    print(f"  Smallest |Phi/c^2|/grav over the whole scan = {min_ratio:.2e}")
    if min_ratio > 1e20:
        print("  => NUCLEAR WALL confirmed: the lapse mechanism is >1e20x too strong")
        print("     across every plausible core density, for BOTH lapse forms.")
        print("     This reproduces the documented BLV result (10^30-10^40x, Yukawa")
        print("     ~0.8 fm). The density->time relation is a NUCLEAR effect, not")
        print("     gravity. Decision gate (PLAN.md §5) -> characterization-and-exit:")
        print("     pin the density->time LAW (Thread A) and document why the magic")
        print("     sector is short-range, rather than launching GPU gravity runs.")
    else:
        print("  => Some regime approaches the gravitational scale: escalate to the")
        print("     full three-thread simulation program (PLAN.md §5).")
    print("=" * 74)
    return rows, min_ratio


def self_checks():
    """Embedded self-checks. Exit 1 on any failure."""
    ok = True
    def check(name, cond):
        nonlocal ok
        status = "PASS" if cond else "FAIL"
        if not cond: ok = False
        print(f"  [{status}] {name}")

    print("\nSELF-CHECKS")
    print("-" * 11)
    # 1. lapse monotone decreasing in P (higher density -> slower time)
    check("alpha decreases with P (n=2)", lapse(0.5, 2) > lapse(2.0, 2))
    # 2. Phi/c^2 negative (potential well) and deepening with P
    check("Phi/c^2 < 0 and deepens with P (n=2)",
          phi_over_c2(1.0, 2) < 0 and phi_over_c2(2.0, 2) < phi_over_c2(1.0, 2))
    # 3. correct sign: -d ln alpha / dP > 0  =>  force up the density gradient (attractive)
    dP = 1e-6
    dlnalpha = (math.log(lapse(1.0 + dP, 2)) - math.log(lapse(1.0, 2))) / dP
    check("-d(ln alpha)/dP > 0 (drift up-gradient = attractive)", -dlnalpha > 0)
    # 4. gravitational yardstick is ~1e-39 (proton)
    check("Phi_grav(proton) ~ 1e-39", 1e-40 < PHI_GRAV_PROTON < 1e-38)
    # 5. nuclear wall: ratio >> 1e20 at representative core
    ratio = abs(phi_over_c2(1.0, 2)) / PHI_GRAV_PROTON
    check("nuclear wall ratio > 1e20 at P=1", ratio > 1e20)
    # 6. conclusion robust to lapse form: both n=1 and n=2 give O(1) |Phi/c^2| at P=1
    check("|Phi/c^2| ~ O(1) for both forms at P=1",
          0.5 < abs(phi_over_c2(1.0, 1)) < 5 and 0.5 < abs(phi_over_c2(1.0, 2)) < 5)
    return ok


if __name__ == "__main__":
    import sys
    report()
    if not self_checks():
        print("\nSELF-CHECKS FAILED")
        sys.exit(1)
    print("\nAll self-checks passed.")
