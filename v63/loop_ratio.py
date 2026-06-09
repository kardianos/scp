#!/usr/bin/env python3
"""
v63/loop_ratio.py  --  PROCEEDING on the one surviving route:
the phase as a dynamical loop-ratio  cos(3φ) = −c₃/(4c₆).

(The gating checks closed c_eff, holonomy, and canonical localization; v62 left
the loop-ratio as the only number-type that can equal the transcendental cos(2/3).)

This script does the first concrete step honestly: compute c₃, c₆ from the
1-loop Coleman-Weinberg V_eff(φ) of the lepton Brannen spectrum (structural
t²=1/2), and test whether the parameter-free 1-loop ratio reproduces
cos(3φ)=cos(2/3)=0.7859 (i.e. φ=2/9).

FINDING (negative, sharp):
  * For every plausible loop FORM, −c₃/(4c₆) is far OUTSIDE [−1,1] (−7.70, +38.1,
    ...) and even sign-unstable across forms.  So at 1-loop the only critical
    points are the COMMENSURATE sin(3φ)=0 minima (φ=0, π/3, ...), which are
    π-rational ≠ 2/9 (proven: v59 PhaseExclusions z3_potential_does_not_select_2_9,
    cos6_potential_does_not_select_2_9).
  * Getting the offset 2/9 requires the cos(6φ) coefficient boosted to
    |c₆| ≈ |c₃|/(4cos(2/3)) ≈ |c₃|/3.14 -- a HIGHER-loop / sub-leading quantity
    that 1-loop suppresses, and that is currently a FREE coefficient (a dial).

So "proceeding" reduces to ONE well-posed, hard question (v59 Step-6):
  does the G₂-equivariant loop (a FINITE sum over the sector subspace of
  Cl(7)_even, NOT the continuum CW) produce c₃/c₆ = −4cos(2/3) with NO free
  coefficient?  Until that is computed, the loop-ratio is a real home for the
  phase but NOT yet a parameter-free prediction.
"""
import math
import numpy as np

SEP = "=" * 78
COS_2_3 = math.cos(2.0 / 3.0)          # 0.785887  -- target for -c3/(4c6)
NEEDED_RATIO = -4 * COS_2_3            # c3/c6 must equal this for phi=2/9
T = math.sqrt(0.5)                     # lepton structural point t^2 = 1/2


def banner(s):
    print(SEP); print(s); print(SEP)


def s_amps(phi, t=T):
    return np.array([1 + 2 * t * math.cos(phi + 2 * math.pi * k / 3) for k in range(3)])


def mlogm(x):
    """x*log x with the finite limit x*log x -> 0 as x -> 0 (the real CW integrand)."""
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    nz = x > 1e-300
    out[nz] = x[nz] * np.log(x[nz])
    return out


def Veff(phi, form):
    s = s_amps(phi)
    m = s * s                                   # lepton mass m = s^2
    if form == "M2logM2":                       # Sum m (log m - 3/2)   (v_eff_loop form)
        return float(np.sum(mlogm(m) - 1.5 * m))
    if form == "M4logM2":                       # Sum m^2 (log m^2 - 3/2)  (standard fermion CW)
        return float(np.sum(mlogm(m * m) - 1.5 * m * m))
    if form == "s4logs4":                       # treat s as the mass: s^4 log s^4
        return float(np.sum(mlogm(s ** 4) - 1.5 * s ** 4))
    raise ValueError(form)


def fourier(form, n, N=20000):
    g = np.linspace(0, 2 * math.pi, N, endpoint=False)
    V = np.array([Veff(p, form) for p in g])
    d = g[1] - g[0]
    return (1.0 / math.pi) * float(np.sum(V * np.cos(n * g)) * d)


# ---------------------------------------------------------------------------
def part1_target():
    banner("Part 1: the target")
    print(f"  critical-point equation:   cos(3φ) = −c₃ / (4 c₆)")
    print(f"  observed φ = 2/9  ⇒  cos(3·2/9) = cos(2/3) = {COS_2_3:.6f}")
    print(f"  ⇒ the parameter-free loop must give  c₃/c₆ = −4 cos(2/3) = {NEEDED_RATIO:.4f}")


def part2_oneloop_fails():
    banner("Part 2: the 1-loop ratio (parameter-free) — does it land on 0.786?")
    rows = []
    print(f"  {'form':10s} {'c3':>12s} {'c6':>12s} {'-c3/(4c6)':>12s}  {'in[-1,1]?':>10s}")
    for form in ["M2logM2", "M4logM2", "s4logs4"]:
        c3, c6 = fourier(form, 3), fourier(form, 6)
        r = -c3 / (4 * c6) if abs(c6) > 1e-12 else float("nan")
        rows.append((form, c3, c6, r))
        print(f"  {form:10s} {c3:12.5f} {c6:12.6f} {r:12.4f}  {str(abs(r) <= 1):>10s}")
    ratios = [r for _, _, _, r in rows]
    print(f"""
  Every form gives −c₃/(4c₆) FAR outside [−1,1] (target {COS_2_3:.3f}), and the
  values are sign-unstable across forms ({ratios[0]:+.2f}, {ratios[1]:+.2f}, ...).
  So 1-loop has NO in-range nontrivial critical point: only the commensurate
  minima sin(3φ)=0 → φ = 0, π/3, 2π/3, ...  These are π-rational, hence ≠ 2/9
  (v59 PhaseExclusions: z3_/cos6_potential_does_not_select_2_9).
  The parameter-free 1-loop ratio does NOT reproduce the phase.
""")
    return rows


def part3_whats_needed(rows):
    banner("Part 3: what an in-range offset requires (and why it's a dial now)")
    c3 = rows[0][1]
    c6_needed = c3 / NEEDED_RATIO
    c6_actual = rows[0][2]
    print(f"  using form M2logM2:  c₃ = {c3:.4f}")
    print(f"  to get cos(3φ)=cos(2/3) need  c₆ = c₃ / (−4cos(2/3)) = {c6_needed:.4f}")
    print(f"  1-loop actually gives          c₆ = {c6_actual:.4f}  "
          f"(off by ×{c6_needed / c6_actual:.1f})")
    print(f"""
  The cos(6φ) coefficient must be ≈ {abs(c6_needed):.3f} (comparable to c₃), but
  1-loop suppresses the higher harmonic (c₆ ≪ c₃).  Boosting it needs a 2-loop /
  sub-leading term whose coefficient is, right now, a FREE fit -- i.e. a dial.
  (This is the empirical 'tilt' ε in v59/algebra/v_eff_loop.py.)
""")
    return {"c6_needed": c6_needed, "c6_actual": c6_actual}


def part4_frontier():
    banner("Part 4: the make-or-break frontier (v59 Step-6)")
    print("""  v59's 'loop' is NOT the continuum Coleman-Weinberg used above; it is a FINITE
  SUM over the G₂-equivariant subspace of Cl(7)_even (v_eff_loop.py docstring).
  A finite, structured sum can have DIFFERENT harmonic content -- c₆ need not be
  suppressed -- so its ratio c₃/c₆ is fixed by the G₂ representation, with no free
  loop scale.

  THE one well-posed question that decides parameter-freeness:

     Does the G₂-equivariant finite-sum loop give  c₃/c₆ = −4 cos(2/3)  with NO
     free coefficient (only the structural counting N_X from the G₂ projection)?

  If YES -> the phase is a genuine parameter-free prediction (the win).
  If it needs a fitted ε / loop order -> the phase is a DIAL (or α-downstream).

  Requirements to compute it (unbuilt):
    1. Build the G₂ ⊂ Spin(7) action on the sector ambient (Λ⁴(R⁷) for d-quark;
       the L⊕F structure for u-quark) and identify the N_X-dim invariant subspace.
    2. Evaluate the finite-sum V_eff(φ) restricted to that subspace.
    3. Extract c₃, c₆ and test c₃/c₆ = −4cos(2/3) with zero free coefficients.

  Honest risk: the 1-loop form-sensitivity above (−7.70 vs +38.1) is a warning
  that the ratio may stay form/scheme dependent even in the finite sum.  This is
  a real frontier, not a sure thing.
""")


def run_checks(rows, need):
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    ratios = [r for _, _, _, r in rows]
    check("1-loop ratio is OUT of [-1,1] for all forms (no in-range crit pt)",
          all(abs(r) > 1 for r in ratios))
    check("ratio is form-sensitive (not a robust prediction)",
          max(ratios) - min(ratios) > 1.0)
    check("ratio is sign-unstable across forms (M2 vs M4)",
          ratios[0] * ratios[1] < 0)
    check("required ratio = -4cos(2/3) ~ -3.144", abs(NEEDED_RATIO + 3.1436) < 1e-3)
    check("needed c6 differs from 1-loop c6 (offset needs extra structure)",
          abs(need["c6_needed"] - need["c6_actual"]) > 1e-2)
    # the commensurate minima sin(3φ)=0 are pi-rational, hence != 2/9 (sanity).
    check("commensurate minimum φ=π/3 != 2/9", abs(math.pi / 3 - 2 / 9) > 0.1)

    print(SEP)
    print("ALL CHECKS PASSED (1-loop loop-ratio does NOT give the phase; frontier scoped)"
          if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    part1_target()
    print()
    rows = part2_oneloop_fails()
    print()
    need = part3_whats_needed(rows)
    print()
    part4_frontier()
    print()
    ok = run_checks(rows, need)
    import sys
    sys.exit(0 if ok else 1)
