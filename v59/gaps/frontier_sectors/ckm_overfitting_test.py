#!/usr/bin/env python3
"""
ckm_overfitting_test.py  —  G10 (CKM) rigor / overfitting analysis.

Goal: stress-test the v59 CKM conjectures against PDG, and HONESTLY quantify
how easy it is to hit the observed values with small-integer×α (or simple
ratio) forms.  This cluster is the most prone to spurious matches; the most
valuable output is a bound on what is and isn't a real signal.

PDG 2024 central values (CKM):
  V_us = sin θ_C = 0.22500   (Cabibbo)
  V_cb            = 0.04182
  V_ub            = 0.00369
  V_td            = 0.00857
  δ_CP            = 1.20 rad   (CKM phase, in the standard parametrization γ≈65.5°→δ via Wolfenstein)
  Wolfenstein:  λ=0.22500, A=0.826, ρ̄=0.159, η̄=0.348
"""
import math
import itertools

# ---------------------------------------------------------------------------
# Empirical inputs
# ---------------------------------------------------------------------------
ALPHA0  = 1.0 / 137.035999084     # α(0), IR / Thomson value
ALPHAMZ = 1.0 / 127.951           # α(M_Z)

V_us = 0.22500
V_cb = 0.04182
V_ub = 0.00369
V_td = 0.00857
delta_CP = 1.20                   # rad
A_wolf  = V_cb / V_us**2          # Wolfenstein A
rhobar, etabar = 0.159, 0.348

def pct(pred, emp):
    return abs(pred - emp) / abs(emp) * 100.0

print("="*74)
print("G10 — CKM:  v59 conjecture tests + overfitting analysis")
print("="*74)

# ---------------------------------------------------------------------------
# 1.  The headline conjecture: sin²θ_C = 7·α(0),  i.e.  sinθ_C = √(7α)
# ---------------------------------------------------------------------------
print("\n--- 1. Cabibbo:  sin²θ_C = 7·α(0)  (7 = dim Im O) ---")
sin2_emp = V_us**2
for label, a in [("α(0)", ALPHA0), ("α(M_Z)", ALPHAMZ)]:
    pred_sin2 = 7 * a
    pred_sin  = math.sqrt(7 * a)
    print(f"  with {label:7s}: 7α = {pred_sin2:.6f}  vs sin²θ_C={sin2_emp:.6f}  "
          f"[{pct(pred_sin2, sin2_emp):.3f}%];  √(7α)={pred_sin:.5f} vs {V_us} "
          f"[{pct(pred_sin, V_us):.3f}%]")

# ---------------------------------------------------------------------------
# 2.  OVERFITTING TEST.  How many (small integer)×α or simple-ratio forms
#     land within X% of sin²θ_C = 0.0506?  If "lots", 7α is not special.
# ---------------------------------------------------------------------------
print("\n--- 2. OVERFITTING: how many simple ×α forms hit sin²θ_C=0.0506? ---")
target = sin2_emp
# v59 structural integers (the documented roster) and small integers
v59_ints = [5, 7, 8, 9, 14, 16, 18, 21, 27, 28, 35, 63, 72]
small    = list(range(1, 16))
# candidate multipliers: n, n/m for small n,m, and v59 ratios
cands = {}
for n in sorted(set(small + v59_ints)):
    cands[f"{n}·α"] = n * ALPHA0
for n in sorted(set(small + v59_ints)):
    for m in range(1, 16):
        if math.gcd(n, m) == 1 and n != m:
            cands[f"({n}/{m})·α"] = (n/m) * ALPHA0
# also n·α(M_Z)
for n in sorted(set(small + v59_ints)):
    cands[f"{n}·α(MZ)"] = n * ALPHAMZ

for thresh in (0.5, 1.0, 2.0, 5.0):
    hits = [(k, v) for k, v in cands.items() if pct(v, target) < thresh]
    print(f"  within {thresh:>3}% of sin²θ_C: {len(hits):3d}/{len(cands)} forms"
          f"   e.g. {sorted(hits, key=lambda kv: pct(kv[1],target))[:4]}")
# Rank where 7α sits
ranked = sorted(cands.items(), key=lambda kv: pct(kv[1], target))
rank7 = next(i for i,(k,v) in enumerate(ranked) if k == "7·α")
print(f"  → '7·α' ranks #{rank7+1} of {len(cands)} by closeness "
      f"({pct(cands['7·α'],target):.3f}%). Best form: {ranked[0]}")
print("  INTERPRETATION: 7·α is one of MANY ×α forms near the target. Its claim to")
print("  significance is STRUCTURAL (7=dimImO, recurrent in v59), not statistical.")

# ---------------------------------------------------------------------------
# 3.  V_cb via the 7/9 (Wolfenstein A = t²_u) conjecture
# ---------------------------------------------------------------------------
print("\n--- 3. V_cb / V_ub / V_td: structural candidates ---")
print(f"  Wolfenstein A (=V_cb/V_us²) empirical = {A_wolf:.4f}")
for label, val in [("7/9 = t²_u = cos²θ_W", 7/9),
                   ("4/5", 4/5),
                   ("(1-2/9)(1-1/14)", (1-2/9)*(1-1/14)),
                   ("28/35 = 4/5", 28/35)]:
    print(f"    A ?= {label:24s} = {val:.4f}   [{pct(val, A_wolf):.2f}%]")
Vcb_pred = (7/9) * 7 * ALPHA0       # A·λ² = (7/9)·(7α)
print(f"  V_cb pred = (7/9)·(7α) = {Vcb_pred:.5f}  vs {V_cb}  [{pct(Vcb_pred,V_cb):.2f}%]")

# V_ub, V_td: in Wolfenstein these depend on (ρ̄,η̄) — i.e. on δ_CP.
# Test whether |ρ̄+iη̄| or the phase have v59 readings.
rho_mag = math.hypot(rhobar, etabar)
print(f"  |ρ̄+iη̄| = {rho_mag:.4f}   candidates:")
for label, val in [("3/8", 3/8), ("sin θ_C·1.7", V_us*1.7), ("√(2/9)·0.8", math.sqrt(2/9)*0.8)]:
    print(f"    ?= {label:14s} = {val:.4f}  [{pct(val, rho_mag):.2f}%]")
Vub_pred = (7/9) * V_us**3 * rho_mag
print(f"  V_ub (Wolfenstein, using empirical |ρ̄+iη̄|) = {Vub_pred:.5f} vs {V_ub} [{pct(Vub_pred,V_ub):.1f}%]")
print("  → V_ub, V_td have NO v59-structural reading independent of empirical (ρ̄,η̄).")

# ---------------------------------------------------------------------------
# 4.  δ_CP — search for any clean form (expect NULL)
# ---------------------------------------------------------------------------
print("\n--- 4. δ_CP CP phase: search for clean form (expect NULL) ---")
print(f"  δ_CP = {delta_CP} rad = {math.degrees(delta_CP):.1f}°")
delta_cands = {
    "2/9·π·? ": None,
    "π·(7/18)": math.pi*7/18,
    "Q=2/3 rad": 2/3,
    "3·(2/9)=2/3 rad": 2/3,
    "π/e": math.pi/math.e,
    "2/9·... no": None,
    "1 rad": 1.0,
    "π/√7": math.pi/math.sqrt(7),
    "atan(η̄/ρ̄)·...": math.atan2(etabar, rhobar),
}
for k, v in delta_cands.items():
    if v is None: continue
    print(f"    {k:18s} = {v:.4f}  [{pct(v, delta_CP):.2f}%]")
# the standard CKM δ in the std parametrization ≈ 65.5° = 1.144 rad (=γ)
print("  γ (unitarity-triangle angle) = atan2(η̄,ρ̄) = "
      f"{math.degrees(math.atan2(etabar,rhobar)):.1f}° = {math.atan2(etabar,rhobar):.4f} rad")
print("  → NO clean v59-structural form for δ_CP. Genuinely FREE.")

# ---------------------------------------------------------------------------
# 5.  Sector-overlap of Brannen eigenvectors → CKM structure (qualitative)
# ---------------------------------------------------------------------------
print("\n--- 5. Brannen-eigenvector overlap (recap of 09/10_*.py result) ---")
print("  ℂ-circulant kernels share the DFT basis → V_CKM = permutation (WRONG).")
print("  ℍ-quaternion ξ in sector-specific slices → mixing ~0.4-0.6 (TOO BIG).")
print("  No principled slice choice yet reproduces the near-diagonal CKM.")

print("\n" + "="*74)
print("VERDICT (G10): one soft signal (sin²θ_C≈7α, 0.45%, structural-7), one weak")
print("(A≈7/9, 6%); V_ub/V_td/δ_CP have NO structural reading — genuinely free.")
print("Overfitting test shows 7α is NOT statistically unique among ×α forms.")
print("="*74)
