#!/usr/bin/env python3
"""
rg_fixedpoint_test.py  —  Avenue C (ALTERNATIVES.md): can RG running explain the
g_W^2 = 5 sqrt(alpha) coincidence WITHOUT it being a law?

Key honest point from README.md §4: the SM-definitional law is g_W^2 = 18 pi alpha
(LINEAR, holds at every scale once sin^2θ_W=2/9 is imposed). The "5 sqrt(alpha)"
form only touches this line at alpha(M_Z). If RG running of the LINEAR law passes
through the v59 value alpha(M_Z)=25/(324pi^2) at the Z scale, then the sqrt-form is
explained as "the value the correct (linear) running happens to cross" -- no new law.

This script does a TRANSPARENT 1-loop QED-only running of alpha(0) -> alpha(M_Z) to
see how close standard running gets to the v59 EW value, and reports the gap. It is
deliberately minimal (QED 1-loop with the standard light-fermion content), NOT a full
SM 2-loop GUT analysis -- it tests the DIRECTION and rough magnitude, and flags what a
real test (avenue C) would need.

NOTE: full alpha(M_Z) running is dominated by hadronic vacuum polarization (non-pert);
the perturbative leptonic piece is computed here for orientation only.
"""
import math
PI = math.pi

ALPHA0     = 1/137.035999084
ALPHA_MZ   = 1/127.951            # PDG MSbar (the empirical target)
ALPHA_V59  = 25/(324*PI**2)       # = 1/127.910 (the v59 conjecture)
MZ         = 91.1876              # GeV
ME, MMU, MTAU = 0.000511, 0.10566, 1.77686   # GeV

def reldev(a,b): return abs(a-b)/abs(b)

print("="*78)
print("Avenue C: does standard running explain the EW-scale value?")
print("="*78)
print(f"  alpha(0)        = 1/{1/ALPHA0:.4f}")
print(f"  alpha(M_Z) PDG  = 1/{1/ALPHA_MZ:.4f}")
print(f"  alpha(M_Z) v59  = 25/(324 pi^2) = 1/{1/ALPHA_V59:.4f}  "
      f"(v59 vs PDG: {reldev(ALPHA_V59,ALPHA_MZ)*100:.4f}%)")

# ---------------------------------------------------------------------------
# Leptonic 1-loop QED running of 1/alpha from 0 to M_Z.
#   d(1/alpha)/d ln(mu) = -(1/(2 pi)) * sum_f Q_f^2 * (2/3) * theta(mu>m_f) ...
# Standard result: Delta(1/alpha)_lep(M_Z) ~ -0.03142*... ; the known leptonic
# contribution to alpha(M_Z)^-1 shift is about -3.1 (1/alpha drops ~3.1 from leptons),
# hadronic ~ -3.6 more, top/W small. We compute the leptonic piece explicitly.
# ---------------------------------------------------------------------------
# Delta(alpha)_lep = (alpha/(3 pi)) * sum_lep [ ln(MZ^2/m_l^2) - 5/3 ], Q_l^2=1
def delta_alpha_lep(MZ):
    s = 0.0
    for m in (ME, MMU, MTAU):
        s += math.log(MZ**2/m**2) - 5.0/3.0
    return (ALPHA0/(3*PI)) * s
da_lep = delta_alpha_lep(MZ)
# alpha(MZ) = alpha(0) / (1 - Delta_alpha)
alpha_mz_lep_only = ALPHA0 / (1 - da_lep)
print("\n[leptonic 1-loop running, orientation only]")
print(f"  Delta(alpha)_lep            = {da_lep:.5f}")
print(f"  alpha(M_Z) leptonic-only    = 1/{1/alpha_mz_lep_only:.4f}")
print(f"  (1/alpha shift from leptons)= {1/ALPHA0 - 1/alpha_mz_lep_only:+.3f}")
print(f"  remaining gap to PDG 127.951: {1/alpha_mz_lep_only - 1/ALPHA_MZ:+.3f} "
      f"(the rest is the hadronic+top+W vacuum polarization, non-perturbative)")

# Orientation: full SM running takes 1/alpha 137.036 -> ~128 at M_Z; the leptonic
# 1-loop piece alone (this script) accounts for ~4 of the ~9-unit drop, the
# hadronic+heavy piece for the rest. PDG quotes MSbar 1/127.951.
print("\n[interpretation]")
print("""  - Standard running of alpha(0)->alpha(M_Z) does land in the 1/128 ballpark
    (leptonic 1-loop alone moves 1/alpha by ~4; full SM reaches ~1/128).
  - So the v59 value 1/127.91 is CONSISTENT with being 'the value standard running
    passes through at M_Z' -- supporting the README §4 thesis that 5 sqrt(alpha) is
    a single-scale coincidence, NOT an independent law.
  - But this does NOT derive alpha: the running takes alpha(0) as INPUT. To DERIVE
    alpha via avenue C one needs a FIXED POINT (beta=0) or a forced boundary value
    from the 5:5:2 unification pattern -- neither is supplied by SM running.""")

# ---------------------------------------------------------------------------
# Fixed-point check: does the v59 Spin(7) embedding give a beta-function zero?
# 1-loop SM beta for the U(1) that becomes EM is POSITIVE (alpha grows in UV) -> no
# IR/UV fixed point in pure SM. Report this as the honest null for the naive version.
# ---------------------------------------------------------------------------
print("\n[fixed-point check, 1-loop SM]")
print("""  SM 1-loop: b0(U(1)_em-like) > 0  =>  alpha grows monotonically in the UV,
  NO Banks-Zaks / asymptotically-safe zero in the pure-SM gauge sector.
  => A genuine fixed point requires NEW matter/structure from the Spin(7) embedding
     (the exceptional-group content). That computation is NOT done here and is the
     real open work for avenue C. Status: DIRECTION supported, DERIVATION open.""")

print("\n" + "="*78)
print(f"VERDICT (avenue C): v59 alpha(M_Z) is consistent with standard running")
print(f"(supports 'sqrt-form = single-scale coincidence'); a true derivation still")
print(f"needs a forced boundary value or a fixed point -- OPEN, highest-leverage.")
print("="*78)
