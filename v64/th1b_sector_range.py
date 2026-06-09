#!/usr/bin/env python3
"""
v64 T-H1b (analytic precheck) — does the long-range 1/r tail live ONLY in the massless
sector, and what does it require?

Hypothesis H1 (PLAN.md §1c): gravity = the MASSLESS sector (m_theta^2=0); strong = the
MASSIVE phi sector (m^2=2.25). Prediction: a localized source produces a long-range 1/r
tail ONLY when the carrier is massless AND the source carries a nonzero monopole. A mass
screens the tail (Yukawa); a total-divergence source (zero monopole) has NO 1/r tail
even when massless (the v6 Gauss obstruction).

This precheck establishes the falsifiable law and the measurement protocol BEFORE any
GPU run (PLAN.md guardrail: prechecks first). It computes the linear screened response
phi of (-nabla^2 + m^2) phi = rho on a radial grid and measures the effective falloff
power  p(r) = d ln|phi| / d ln r  in a far-field window, for:
  (A) carrier mass scan  m^2 = 0, small, ..., 2.25   (massless vs massive)
  (B) source-monopole structure: rho_B = Gaussian (nonzero monopole, the v6 |omega|^2
      channel) vs rho_A = nabla^2(Gaussian) (ZERO monopole, the v6 algebraic/Gauss
      obstruction). Massless response: B -> 1/r tail; A -> NO 1/r tail.

Analytic, no simulation. Embedded self-checks (exit 1 on fail).

Refs: PLAN.md §1c (H1, T-H1b), MEMORY 'V6 — Conserved Torsion Field'
(massless wave: box delta-rho = -1/2|omega|^2 -> 1/r; algebraic: -> 1/r^6, zero
monopole; G_eff/G_N ~ 3.6e40), CLAUDE.md (m^2=2.25 phi, m_theta^2=0).
"""
import math

CODE_L_FM = 0.5624

# ---------------------------------------------------------------------------
# Screened linear response of (-nabla^2 + m^2) phi = rho for a spherically symmetric
# source rho(r), via the radial Green's function G(r,r') of the modified Helmholtz op.
# For m>0:  phi(r) = (1/(m r)) \int_0^inf r' rho(r') [e^{-m|r-r'|} - e^{-m(r+r')}]/2 dr'
# For m=0 (Poisson): phi(r) = (1/r)\int_0^r r'^2 rho dr' + \int_r^inf r' rho dr'
# We use a unified small-m form; m=0 handled by the Poisson branch.
# ---------------------------------------------------------------------------
def response(r, m, rho, hi, n=6000):
    dr = hi/n
    if m < 1e-9:  # Poisson
        inner = 0.0; outer = 0.0
        for i in range(1, n+1):
            rp = i*dr
            if rp <= r: inner += rp*rp*rho(rp)*dr
            else:       outer += rp*rho(rp)*dr
        return inner/max(r,1e-9) + outer
    s = 0.0
    for i in range(1, n+1):
        rp = i*dr
        s += rp*rho(rp)*(math.exp(-m*abs(r-rp)) - math.exp(-m*(r+rp)))*dr
    return s/(2.0*m*max(r,1e-9))

def falloff_power(m, rho, hi, r1, r2, npts=10):
    """effective power p = d ln|phi| / d ln r over [r1,r2]."""
    rs = [r1*(r2/r1)**(i/(npts-1)) for i in range(npts)]
    lr = [math.log(r) for r in rs]
    lp = [math.log(abs(response(r, m, rho, hi))) for r in rs]
    mlr = sum(lr)/npts; mlp = sum(lp)/npts
    num = sum((lr[i]-mlr)*(lp[i]-mlp) for i in range(npts))
    den = sum((lr[i]-mlr)**2 for i in range(npts))
    return num/den

# Sources
W = 0.6
def rho_gauss(r):            # nonzero monopole: Q = \int rho dV != 0  (v6 |omega|^2 channel)
    return math.exp(-(r/W)**2)
def rho_divergence(r):       # rho = nabla^2(Gaussian): ZERO monopole (total divergence)
    # nabla^2 e^{-(r/W)^2} = (4 r^2/W^4 - 6/W^2) e^{-(r/W)^2}
    g = math.exp(-(r/W)**2)
    return (4.0*r*r/W**4 - 6.0/W**2)*g

def monopole(rho, hi, n=8000):
    dr = hi/n
    return sum((i*dr)**2 * rho(i*dr) * dr for i in range(1, n+1))  # \int r^2 rho dr (~Q/4pi)

def report():
    print("=" * 72)
    print("v64 T-H1b (precheck) — long-range tail lives only in the massless sector?")
    print("=" * 72)
    HI = 14.0
    # window in the far field, a few code lengths out
    R1, R2 = 4.0, 8.0
    print(f"Far-field window r in [{R1}, {R2}] code lengths "
          f"({R1*CODE_L_FM:.1f}-{R2*CODE_L_FM:.1f} fm); source width W={W}")
    print()

    # ---- (A) carrier-mass scan, monopole source ----
    print("(A) Falloff power p = d ln|phi|/d ln r vs carrier mass (monopole source):")
    print(f"    {'m^2':>7} {'m':>7} | {'power p':>9}   interpretation")
    print("    " + "-" * 56)
    m2_list = [0.0, 0.01, 0.04, 0.25, 1.0, 2.25]
    powers = {}
    for m2 in m2_list:
        m = math.sqrt(m2)
        p = falloff_power(m, rho_gauss, HI, R1, R2)
        powers[m2] = p
        if m2 == 0.0:
            interp = "1/r LONG-RANGE (gravity carrier)"
        elif p < -2.5:
            interp = "screened / Yukawa (short-range)"
        else:
            interp = "near-massless, tail still ~1/r in window"
        print(f"    {m2:>7.2f} {m:>7.3f} | {p:>9.3f}   {interp}")
    print()
    print(f"  massless p = {powers[0.0]:+.3f} (-> -1, pure 1/r);  "
          f"standard m^2=2.25 p = {powers[2.25]:+.3f} (steep/screened)")
    print()

    # ---- (B) source-monopole structure (the v6 Gauss bypass), massless carrier ----
    # The 1/r tail COEFFICIENT is the monopole Q; existence of a tail is an AMPLITUDE
    # question, not a power-law one. (A zero-monopole source has no 1/r term; its tiny
    # residual response is numerical roundoff that can spuriously refit to power -1.)
    print("(B) Massless (m=0) response vs SOURCE monopole structure (tail amplitude):")
    qB = monopole(rho_gauss, HI); qA = monopole(rho_divergence, HI)
    rmid = 6.0
    aB = abs(response(rmid, 0.0, rho_gauss, HI))
    aA = abs(response(rmid, 0.0, rho_divergence, HI))
    print(f"    rho_B = Gaussian (v6 |omega|^2):    monopole Q ~ {qB:+.4f}  "
          f"|phi(r={rmid:.0f})| = {aB:.3e}")
    print(f"    rho_A = nabla^2(Gaussian):          monopole Q ~ {qA:+.3e}  "
          f"|phi(r={rmid:.0f})| = {aA:.3e}")
    print(f"    amplitude ratio A/B in far field = {aA/aB:.2e}  (<<1 => A has NO tail)")
    print("    => a 1/r tail requires BOTH masslessness AND a nonzero-monopole source.")
    print("       A total-divergence source (Q=0) has NO 1/r tail even when massless")
    print("       (the v6 Gauss obstruction: algebraic channel -> 1/r^6, not gravity).")
    print()

    print("VERDICT")
    print("-" * 7)
    massless_is_1r   = abs(powers[0.0] + 1.0) < 0.15
    massive_screened = powers[2.25] < -2.5
    tail_needs_monopole = (abs(qA) < 1e-2*abs(qB)) and (aA/aB < 1e-3)
    if massless_is_1r and massive_screened and tail_needs_monopole:
        print("  => CONSISTENT with H1: the long-range 1/r tail lives ONLY in the")
        print("     massless sector, and ONLY for a monopole-carrying source. The")
        print("     massive phi sector (m^2=2.25) is Yukawa-screened (strong/nuclear);")
        print("     the massless sector with a |omega|^2-type monopole source carries")
        print("     the 1/r (gravity-range) tail. This is the v6 result, recovered as")
        print("     the falsifiable PROTOCOL for the GPU run:")
        print("       * measure radial power of the theta-sector response -> expect ~ -1")
        print("       * measure radial power of the phi-sector response   -> expect << -1")
        print("       * confirm the theta source has a nonzero monopole in real dynamics")
        print("     OPEN (needs the sim, not analytics): does the REAL nonlinear theta")
        print("     dynamics actually generate a nonzero-monopole source? (v6 claims yes")
        print("     via Q_eff=E2/rho0.) And magnitude stays ~10^40 too strong.")
    else:
        print("  => prediction not cleanly recovered; revisit before the sim.")
    print("=" * 72)
    return powers, (qA, qB, aA, aB)

def self_checks():
    ok = True
    def check(name, cond):
        nonlocal ok
        if not cond: ok = False
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    print("\nSELF-CHECKS")
    print("-" * 11)
    HI = 14.0; R1, R2 = 4.0, 8.0
    # 1. massless monopole source -> 1/r (power ~ -1)
    p0 = falloff_power(0.0, rho_gauss, HI, R1, R2)
    check("massless + monopole source -> power ~ -1", abs(p0 + 1.0) < 0.15)
    # 2. standard-mass response is screened (much steeper than -1)
    pm = falloff_power(math.sqrt(2.25), rho_gauss, HI, R1, R2)
    check("massive (m^2=2.25) response screened (power < -2.5)", pm < -2.5)
    # 3. monotone: heavier carrier -> steeper falloff in window
    p_lo = falloff_power(0.1, rho_gauss, HI, R1, R2)
    p_hi = falloff_power(1.0, rho_gauss, HI, R1, R2)
    check("heavier carrier -> steeper falloff", p_hi < p_lo)
    # 4. divergence source has ~zero monopole; Gaussian has nonzero
    qA = monopole(rho_divergence, HI); qB = monopole(rho_gauss, HI)
    check("nabla^2(Gaussian) has ~zero monopole", abs(qA) < 1e-2*abs(qB))
    check("Gaussian has nonzero monopole", abs(qB) > 0.05)
    # 5. massless response to zero-monopole source has NO 1/r tail: its far-field
    #    amplitude is negligible vs the monopole source (the tail coefficient is Q).
    aA = abs(response(6.0, 0.0, rho_divergence, HI))
    aB = abs(response(6.0, 0.0, rho_gauss, HI))
    check("massless + zero-monopole source -> negligible tail amplitude (A/B < 1e-3)",
          aA/aB < 1e-3)
    return ok

if __name__ == "__main__":
    import sys
    report()
    if not self_checks():
        print("\nSELF-CHECKS FAILED"); sys.exit(1)
    print("\nAll self-checks passed.")
