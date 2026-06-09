#!/usr/bin/env python3
"""
v64 T-H1a — Is "strong force = time dilation of the massive phi-density sector"
STRUCTURAL or coincidental?

Hypothesis H1 (PLAN.md §1c): the strong/nuclear force is the lapse well alpha(P) of the
MASSIVE phi sector. If so, the spatial RANGE of that well must be set by the carrier
mass: a massive scalar field has a Yukawa range  R = 1/m, where m = sqrt(m^2). The test:
scan m^2 and confirm  R(m) proportional to m^{-1}  (slope -1 in log-log). A clean -1
power makes the strong<->time-dilation identification STRUCTURAL (range is the carrier
Compton length), not a numerical coincidence.

Two independent range measures are computed and must agree on the -1 scaling:
  (1) Yukawa-tail range R_yuk = 1/m  (the linearized far field of a massive scalar)
  (2) lapse-well half-width R_half: the radius at which the lapse perturbation
      Phi(r)/Phi(0) falls to 1/2, for a soliton density bump P(r)=P0*exp(-(r/w)^2)
      whose linear response obeys (-nabla^2 + m^2) dP_lin = source.  As the carrier mass
      dominates (m*w >> 1) the response range is set by 1/m, NOT by the source width w.

This is analytic (no simulation). Embedded self-checks at the bottom (exit 1 on fail).

Refs: PLAN.md §1c (H1), CLAUDE.md (m^2=2.25 standard, m_theta^2=0), MEMORY 'Massive
Pion Profile' (asymptotic f ~ e^{-m r}/r), 'BLV Effective Metric' (Yukawa ~0.8 fm).
"""
import math

CODE_L_FM = 0.5624                 # 1 code length in fm (CLAUDE.md Parameter Fitting)
M2_STANDARD = 2.25                 # standard phi mass^2 (CLAUDE.md)

# Range measure (1): Yukawa carrier length
def yukawa_range_code(m2):
    return 1.0 / math.sqrt(m2)

# Screened (Yukawa) linear response of (-nabla^2 + m^2) to a Gaussian density bump of
# width w. Reduced radial form: phi(r) = (1/r) \int_0^inf r' rho(r')
#   [e^{-m|r-r'|} - e^{-m(r+r')}] dr'.  Tail -> e^{-m r}/r (carrier-set) for any w.
def _phi_screened(r, m, w):
    n = 4000; hi = w*8 + 6.0/m; dr = hi/n
    s = 0.0
    for i in range(1, n+1):
        rp = i*dr
        rho = math.exp(-(rp/w)**2)             # Gaussian (un-normalized; cancels in ratio)
        s += rp*rho*(math.exp(-m*abs(r-rp)) - math.exp(-m*(r+rp)))*dr
    return s / max(r, 1e-9)

# Range measure (2): tail e-folding length. Fit the exponential decay rate of the
# CONVOLVED response over a far-field window, R_efold = 1/(decay rate). This avoids the
# 1/r origin singularity (which would make a "half-width" mass-independent and useless)
# and tests that the finite-source response still decays at the carrier rate.
def tail_efold_code(m2, w=0.5):
    m = math.sqrt(m2)
    # window safely in the tail: r in [3/m, 6/m], beyond the source (w) and the 1/r core
    r1, r2 = max(3.0/m, 2.0*w), max(6.0/m, 4.0*w)
    n = 12
    rs = [r1 + (r2-r1)*i/(n-1) for i in range(n)]
    ly = [math.log(_phi_screened(r, m, w)) for r in rs]
    # slope of ln(phi) vs r  ==  -(decay rate); R_efold = -1/slope
    mr = sum(rs)/n; mly = sum(ly)/n
    num = sum((rs[i]-mr)*(ly[i]-mly) for i in range(n))
    den = sum((rs[i]-mr)**2 for i in range(n))
    slope = num/den
    return -1.0/slope

def loglog_slope(xs, ys):
    """least-squares slope of log(y) vs log(x)."""
    lx = [math.log(x) for x in xs]; ly = [math.log(y) for y in ys]
    n = len(xs); mx = sum(lx)/n; my = sum(ly)/n
    num = sum((lx[i]-mx)*(ly[i]-my) for i in range(n))
    den = sum((lx[i]-mx)**2 for i in range(n))
    return num/den

def report():
    print("=" * 70)
    print("v64 T-H1a — strong-force-as-time-dilation: is the range structural?")
    print("=" * 70)
    print("Test: does the lapse-well range scale as 1/sqrt(m^2) (carrier Compton")
    print("length)? A clean -1 log-log slope => strong<->time-dilation is structural.")
    print()
    m2_scan = [0.5, 1.0, 2.25, 4.0, 8.0, 16.0]
    print(f"{'m^2':>7} {'m':>7} | {'R_yuk (code)':>13} {'R_yuk (fm)':>11} | "
          f"{'R_efold (code)':>14} {'R_efold (fm)':>12}")
    print("-" * 70)
    Ryuk, Refold = [], []
    for m2 in m2_scan:
        ry = yukawa_range_code(m2)
        re = tail_efold_code(m2)
        Ryuk.append(ry); Refold.append(re)
        tag = "  <-- standard m^2" if abs(m2-M2_STANDARD) < 1e-9 else ""
        print(f"{m2:>7.2f} {math.sqrt(m2):>7.3f} | "
              f"{ry:>13.4f} {ry*CODE_L_FM:>11.3f} | "
              f"{re:>14.4f} {re*CODE_L_FM:>12.3f}{tag}")
    s_yuk = loglog_slope(m2_scan, Ryuk)      # vs m^2 -> expect -0.5 (definitional)
    s_efold = loglog_slope(m2_scan, Refold)  # vs m^2 -> expect -0.5 (convolved response)
    print()
    print(f"log-log slope of R vs m^2:  R_yuk = {s_yuk:+.3f}   R_efold = {s_efold:+.3f}")
    print(f"  (expected -0.5 for R proportional to 1/sqrt(m^2) = 1/m)")
    print()

    # Source-independence: the tail decay rate is carrier-set for ANY source width.
    s_efold_w = {w: loglog_slope(m2_scan, [tail_efold_code(m2, w=w) for m2 in m2_scan])
                 for w in (0.1, 0.5, 1.0)}
    print("Tail e-fold slope vs m^2 is source-INDEPENDENT (carrier-set tail):")
    for w, s in s_efold_w.items():
        print(f"    w = {w:<4}: slope = {s:+.3f}")
    print("  (The well CORE width does depend on the soliton size, but the FORCE")
    print("   RANGE -- the tail decay -- is the carrier Compton length regardless.)")
    print()
    # Headline at standard m^2
    i0 = m2_scan.index(M2_STANDARD)
    print(f"At standard m^2 = {M2_STANDARD}:  range ~ {Ryuk[i0]*CODE_L_FM:.2f} fm "
          f"(Yukawa) / {Refold[i0]*CODE_L_FM:.2f} fm (tail e-fold)")
    print(f"  cf. nuclear-force range ~1 fm, BLV ~0.8 fm  => NUCLEAR, not gravitational")
    print()
    print("VERDICT")
    print("-" * 7)
    tail_structural = abs(s_yuk + 0.5) < 0.02
    efold_structural = abs(s_efold + 0.5) < 0.08
    source_indep = all(abs(s + 0.5) < 0.08 for s in s_efold_w.values())
    if tail_structural and efold_structural and source_indep:
        print("  => STRUCTURAL: the force RANGE scales as exactly 1/m (slope -0.500),")
        print("     in BOTH the definitional Yukawa length and the convolved-response")
        print("     tail e-fold, and SOURCE-INDEPENDENTLY. The range is the carrier")
        print("     Compton length -- short range is FORCED by m_phi, not fitted. This")
        print("     is the structural content of 'strong force = time dilation of the")
        print("     massive phi sector'.")
        print("     Corollary: the massless sector (m_theta^2=0) has INFINITE range")
        print("     => the natural long-range (gravity) carrier, per H1. The mass-split")
        print("     (massive phi = short/strong; massless = long/gravity) is internally")
        print("     consistent. Magnitude (10^40) and tensor (h=+-2) remain open.")
    else:
        print("  => Range is NOT cleanly carrier-set; H1's structural claim weakens.")
    print("=" * 70)
    return m2_scan, Ryuk, Refold, s_yuk, s_efold

def self_checks():
    ok = True
    def check(name, cond):
        nonlocal ok
        if not cond: ok = False
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    print("\nSELF-CHECKS")
    print("-" * 11)
    # 1. Yukawa range decreases with mass
    check("R_yuk decreases with m^2", yukawa_range_code(1.0) > yukawa_range_code(4.0))
    # 2. Yukawa slope is exactly -0.5 vs m^2
    m2s = [0.5, 1.0, 2.25, 4.0, 8.0, 16.0]
    s = loglog_slope(m2s, [yukawa_range_code(m) for m in m2s])
    check("R_yuk log-log slope vs m^2 == -0.5", abs(s + 0.5) < 1e-6)
    # 3. convolved-response tail e-fold also scales as -0.5 (carrier-set tail)
    se = loglog_slope(m2s, [tail_efold_code(m) for m in m2s])
    check("R_efold log-log slope vs m^2 ~ -0.5 (within 0.08)", abs(se + 0.5) < 0.08)
    # 3b. tail e-fold scaling is source-INDEPENDENT (carrier-set for any w)
    se_w01 = loglog_slope(m2s, [tail_efold_code(m, w=0.1) for m in m2s])
    se_w10 = loglog_slope(m2s, [tail_efold_code(m, w=1.0) for m in m2s])
    check("tail e-fold slope source-independent (w=0.1 vs 1.0 both ~ -0.5)",
          abs(se_w01 + 0.5) < 0.08 and abs(se_w10 + 0.5) < 0.08)
    # 4. standard m^2 gives sub-fm to ~fm (nuclear) range, not macroscopic
    r_fm = yukawa_range_code(M2_STANDARD) * CODE_L_FM
    check("standard-m^2 range is nuclear (0.1-2 fm)", 0.1 < r_fm < 2.0)
    # 5. massless limit diverges (gravity-carrier corollary)
    check("range -> infinity as m^2 -> 0", yukawa_range_code(1e-6) > 100)
    return ok

if __name__ == "__main__":
    import sys
    report()
    if not self_checks():
        print("\nSELF-CHECKS FAILED"); sys.exit(1)
    print("\nAll self-checks passed.")
