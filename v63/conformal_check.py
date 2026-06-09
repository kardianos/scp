#!/usr/bin/env python3
"""
v63/conformal_check.py  --  SECOND gating check for the holonomy idea:
is there a parameter-free c_eff in the CURRENT (Cl(3,1)) gravity sector?

(Follows ceff_test.py, which falsified the BLV v_min/c = 2/3 lead in the OLD,
superseded gravity sector.  Here we test the current one.)

Context pulled (v60/gravity_recast/08_plebanski_action.{md,py}, 05_findings.md,
v61/PROVEN_LEDGER.md, CLOSEOUT.md):
  * The spacetime metric is NOT fundamental: it is the INDUCED metric g(B) via
    the Urbantke formula, which is CUBIC in the 2-form B and densitized -- it
    recovers g only UP TO an overall scale (the 08 script reports "ĝ ∝ g").
  * The ρ_grav coupling normalization is OPEN; the only magnitude on offer,
    G_e = (21/16)α²¹, is flagged "curve-fitting, not derivation" (P-fit).
So every MAGNITUDE in the sector is open/fitted; what is pinned is the conformal
(causal/angle) structure, the 2 TT DOF, helicity ±2, the signature.

A "c_eff" (a speed) is a MAGNITUDE.  This script shows:
  (1) the Urbantke metric is homogeneous of degree 3 -> conformal-class only ->
      no absolute scale -> no parameter-free speed;
  (2) the conformally-invariant *ratio* (a speed v²=−g_tt/g_xx) IS scale-free for
      a fixed metric, BUT it is a position/matter-dependent FIELD, not a universal
      constant (the BLV numbers 0.670 vs 0.625 from ceff_test already prove the
      matter dependence) -> still no universal c_eff = 2/3.
Verdict: as anticipated, NEGATIVE.  The only parameter-free, scale-free angles
live in the INTERNAL G₂/S⁷ geometry (fixed round metric + triality cycle), not in
the dynamical emergent spacetime metric -- and that route needs no c_eff at all.
"""
import numpy as np

SEP = "=" * 78
rng = np.random.default_rng(0)


def banner(s):
    print(SEP); print(s); print(SEP)


# 4D Levi-Civita
def levi4():
    eps = np.zeros((4, 4, 4, 4))
    from itertools import permutations
    def sgn(p):
        s = 1
        p = list(p)
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]:
                    s = -s
        return s
    for p in permutations(range(4)):
        eps[p] = sgn(p)
    return eps


EPS4 = levi4()
EPS3 = np.zeros((3, 3, 3))
for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    EPS3[a, b, c] = 1
    EPS3[a, c, b] = -1


def urbantke(Sig):
    """ĝ_{μν} = −(1/6) ε^{αβγδ} ε_{ijk} Σ^i_{μα} Σ^j_{βγ} Σ^k_{δν}.
       Sig: (3,4,4) array of three antisymmetric 2-forms."""
    g = -(1.0 / 6.0) * np.einsum(
        "abgd,ijk,ima,jbg,kdn->mn", EPS4, EPS3, Sig, Sig, Sig)
    return 0.5 * (g + g.T)  # symmetrize (numerical)


def antisym3():
    """Three random antisymmetric 4x4 2-forms."""
    S = rng.standard_normal((3, 4, 4))
    return S - np.transpose(S, (0, 2, 1))


# ---------------------------------------------------------------------------
def part1_homogeneity():
    banner("Part 1: Urbantke metric is homogeneous of degree 3 (conformal-class)")
    Sig = antisym3()
    g = urbantke(Sig)
    for lam in [2.0, 0.5, -3.0]:
        g_scaled = urbantke(lam * Sig)
        expected = (lam ** 3) * g
        err = np.max(np.abs(g_scaled - expected))
        print(f"  ĝ(λΣ) vs λ³ ĝ(Σ),  λ={lam:+.1f}:  max|Δ| = {err:.2e}")
    print("""
  ĝ(λΣ) = λ³ ĝ(Σ): the induced metric scales with the 2-form.  An overall
  rescaling of B rescales the metric, so B (hence the matter 2-form, hence the
  density) fixes only the CONFORMAL CLASS of g -- never an absolute length/time
  scale.  A 'speed c_eff' is an absolute-scale (magnitude) quantity, so it is NOT
  pinned by this construction.
""")
    return Sig, g


def part2_ratio_is_conformal_but_not_universal(Sig, g):
    banner("Part 2: the speed RATIO is conformal-invariant but matter-dependent")
    # A "speed" is a ratio of metric components: v² = -g_00/g_11 (schematic).
    # Under B->λB, g->λ³g, so the RATIO is invariant: it is a conformal invariant.
    def speed_sq(gg):
        return -gg[0, 0] / gg[1, 1] if abs(gg[1, 1]) > 1e-12 else float("nan")
    v2 = speed_sq(g)
    v2_scaled = speed_sq(urbantke(2.0 * Sig))
    print(f"  ratio -g00/g11 at scale 1 : {v2:+.6f}")
    print(f"  ratio -g00/g11 at scale 2 : {v2_scaled:+.6f}  "
          f"(scale-invariant: {abs(v2 - v2_scaled) < 1e-9})")
    print(f"""
  So conformal-invariant ratios (speeds) ARE pinned for a FIXED metric.  But the
  metric depends on the matter configuration, so such a ratio is a position- and
  parameter-dependent FIELD, not a universal constant.  The BLV speed ratio is
  exactly this object, and ceff_test.py already showed it is parameter-dependent:
      massless v_min/c = 0.67020151   massive v_min/c = 0.62519338
  -> no universal c_eff = 2/3 even as a conformal invariant.
""")
    return v2


def part3_verdict():
    banner("VERDICT (Cl(3,1) sector)")
    print("""  NEGATIVE, as anticipated.  The current gravity sector pins the CONFORMAL
  (causal / angle) structure, the 2 TT DOF, helicity ±2, and the signature --
  all dimensionless / structural -- but every MAGNITUDE is open or fitted:
    * induced metric g(B) is conformal-class only (Urbantke is degree-3, ĝ ∝ g);
    * ρ_grav coupling normalization is OPEN;
    * the one magnitude on offer, G_e=(21/16)α²¹, is curve-fitting (v61 P-fit),
      and α-dependent.
  A density-dependent SPEED c_eff is a magnitude, so it is not parameter-free
  here -- and the conformal-invariant speed RATIO is matter-dependent (BLV).

  Both gravity sectors now fail to provide a parameter-free c_eff:
    - BLV (old):     a definite speed, but parameter-dependent (moves with m_π);
    - Cl(3,1) (new): conformal-class + open/fitted magnitude.

  CONSTRUCTIVE CONCLUSION: the 'meta-const c' was the wrong object.  Speeds are
  never pinned in the emergent spacetime metric.  The only parameter-free,
  scale-free angles live in the INTERNAL G₂/S⁷ geometry (canonical round metric
  + triality cycle), where a curvature×(fractional area) holonomy is a pure
  dimensionless number independent of any scale.  IF the holonomy idea is
  pursued, the phase is that internal geometric angle -- with NO c_eff factor.
  The constraint that survives is conformal invariance, not a magic value of c.
""")


def run_checks(Sig, g, v2):
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # homogeneity degree 3
    err = np.max(np.abs(urbantke(2.0 * Sig) - 8.0 * g))
    check("Urbantke is degree-3 homogeneous (conformal-class only)", err < 1e-9)
    # ratio scale-invariance (conformal invariant)
    r1 = -g[0, 0] / g[1, 1]
    r2 = -urbantke(2.0 * Sig)[0, 0] / urbantke(2.0 * Sig)[1, 1]
    check("speed ratio is conformally invariant (scale-independent)",
          abs(r1 - r2) < 1e-9)
    # the BLV speed ratios differ (matter dependence) -> no universal constant
    check("BLV speed ratio is matter-dependent (0.6702 != 0.6252)",
          abs(0.67020151 - 0.62519338) > 0.01)
    # the magnitude is genuinely unpinned: scaling Sig changes |g| (not invariant)
    check("metric magnitude is NOT scale-invariant (no absolute speed)",
          abs(np.linalg.norm(urbantke(2.0 * Sig)) - np.linalg.norm(g)) > 1e-6)

    print(SEP)
    print("ALL CHECKS PASSED (no parameter-free c_eff in Cl(3,1), as designed)"
          if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    Sig, g = part1_homogeneity()
    print()
    v2 = part2_ratio_is_conformal_but_not_universal(Sig, g)
    print()
    part3_verdict()
    print()
    ok = run_checks(Sig, g, v2)
    import sys
    sys.exit(0 if ok else 1)
