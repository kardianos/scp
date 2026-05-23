"""
v59 SCALE BRIDGE INVESTIGATION
==============================

Question (Frontier 4 / Q4-1 of ROADMAP):  What structural factor connects the
Brannen lepton scale a to the SM Higgs VEV v?

Key insight reframed:
  - The ROADMAP states "(v/a)^2 ≈ 1.92e8" — but a has units √mass and v has
    units mass, so the dimensionless ratio is v / a^2 (or equivalently (v/a)^2 / v).
  - v / a^2 ≈ 246,220 MeV / 313.84 MeV ≈ 784.6
  - This is suspicious: 28^2 = 784.

If the conjecture v_Higgs = D_lepton^2 · a^2 holds, then v59 derives the
Higgs VEV from the lepton sector (Brannen scale a) and the ambient
dimension (28).

We also test downstream implications:
  - m_W from g_W^2 = 5·√α (v59 conjecture) + this scale bridge
  - Weak-mixing-angle candidates
  - Higgs mass candidates
"""

import json
import math
import os

# ---------------- empirical inputs (PDG 2024 / CODATA 2022) ----------------
# Charged lepton pole masses (MeV)
M_E   = 0.51099895069
M_MU  = 105.6583755
M_TAU = 1776.86

# Electroweak
V_HIGGS_GEV = 246.2196          # v = (sqrt(2) G_F)^(-1/2); CODATA/PDG
V_HIGGS_MEV = V_HIGGS_GEV * 1000.0
M_W_GEV     = 80.3692            # PDG 2024 average (was 80.379 pre-CDF)
M_Z_GEV     = 91.1876
M_H_GEV     = 125.20             # Higgs mass (PDG 2024)
M_TOP_GEV   = 172.57             # top quark pole mass

ALPHA_INV   = 137.035999084      # 1/alpha at q=0
ALPHA       = 1.0 / ALPHA_INV
ALPHA_S_MZ  = 0.1180             # alpha_s at M_Z

# Weinberg angle (PDG 2024)
SIN2_THW_ONSHELL = 0.22305       # sin^2(theta_W) on-shell = 1 - m_W^2/m_Z^2
SIN2_THW_MSBAR   = 0.23121       # sin^2(theta_W) MS-bar

# Brannen lepton scale (compute from lepton masses)
A_LEPTON = (math.sqrt(M_E) + math.sqrt(M_MU) + math.sqrt(M_TAU)) / 3.0  # √MeV
A_SQ     = A_LEPTON * A_LEPTON                                          # MeV

# ---------------- v59 structural numbers ----------------
# All integers identified in v59 as structural (not empirical):
STRUCT = {
    "7":          7,    # dim ImO = dim S^7
    "8":          8,    # dim O
    "14":        14,    # dim G_2
    "16":        16,    # dim Cl(3,1)
    "21":        21,    # dim Spin(7) = (7 choose 2)
    "28":        28,    # D_lepton = dim(Lambda^2 + Lambda^6) of R^7 = dim Spin(8)
    "35":        35,    # D_d_quark = dim Lambda^3(R^7) = dim Lambda^4(R^7)
    "63":        63,    # D_u_quark = D_lepton + D_d_quark = Lambda^2+Lambda^4+Lambda^6
    "64":        64,    # dim Cl(7)_even = dim(C tensor O)
    "5":          5,    # Killing index so(3) in so(7), also 21-16
    "23":        23,    # Koide Q_u numerator = 23
    "11":        11,    # Koide Q_d numerator = 11
    "2":          2,
    "3":          3,
    "4":          4,
}

print("=" * 70)
print("v59 SCALE BRIDGE INVESTIGATION")
print("=" * 70)
print(f"Brannen lepton scale a   = {A_LEPTON:.6f} √MeV")
print(f"Brannen lepton scale a^2 = {A_SQ:.4f} MeV")
print(f"Higgs VEV v               = {V_HIGGS_MEV:.1f} MeV = {V_HIGGS_GEV:.4f} GeV")
print(f"Dimensionless ratio v/a^2 = {V_HIGGS_MEV/A_SQ:.4f}")
print()

# ---------------- TEST 1: v / a^2 vs structural integer squares ----------------
print("-" * 70)
print("TEST 1: scan v/a^2 against structural-integer squares N^2")
print("-" * 70)
TARGET_1 = V_HIGGS_MEV / A_SQ
print(f"Target v/a^2 = {TARGET_1:.4f}")
print(f"{'N':>5} {'N^2':>8} {'gap %':>10}")
candidates_1 = []
for name, N in STRUCT.items():
    Nsq = N*N
    gap = abs(Nsq - TARGET_1) / TARGET_1 * 100
    candidates_1.append((name, N, Nsq, gap))
candidates_1.sort(key=lambda x: x[3])
for name, N, Nsq, gap in candidates_1[:8]:
    print(f"{name:>5} {Nsq:>8} {gap:>9.3f}%")
print()
print(f"BEST: N = {candidates_1[0][1]} ({candidates_1[0][0]}), N^2 = {candidates_1[0][2]}, gap = {candidates_1[0][3]:.4f}%")
print()

# ---------------- TEST 2: scan multiplicative combinations ----------------
print("-" * 70)
print("TEST 2: scan v/a^2 against products A * B of structural integers")
print("-" * 70)
TARGET_2 = V_HIGGS_MEV / A_SQ
prods = []
for nameA, A in STRUCT.items():
    for nameB, B in STRUCT.items():
        P = A * B
        if P == 0:
            continue
        gap = abs(P - TARGET_2) / TARGET_2 * 100
        prods.append((nameA, A, nameB, B, P, gap))
prods.sort(key=lambda x: x[5])
print(f"{'A':>5} {'B':>5} {'A*B':>8} {'gap %':>10}")
seen = set()
shown = 0
for nameA, A, nameB, B, P, gap in prods:
    key = tuple(sorted([A, B]))
    if key in seen:
        continue
    seen.add(key)
    print(f"{nameA:>5} {nameB:>5} {P:>8} {gap:>9.4f}%")
    shown += 1
    if shown >= 6:
        break
print()
print("→ 28*28 = 784 is the cleanest hit; 16*49 = 7^2 * 16 is the same number.")
print()

# ---------------- TEST 3: predict m_W from scale bridge + g_W^2 = 5√α ----------------
print("-" * 70)
print("TEST 3: predict m_W from (v = 28^2 · a^2) and (g_W^2 = 5 · √α)")
print("-" * 70)
# v59 conjectures (so far)
D_LEPTON = 28
v_v59 = D_LEPTON * D_LEPTON * A_SQ            # predicted Higgs VEV (MeV)
g_W_sq_v59 = 5.0 * math.sqrt(ALPHA)            # v59 conjecture
g_W_v59 = math.sqrt(g_W_sq_v59)

# m_W = g_W · v / 2  (standard SM tree-level)
m_W_v59 = g_W_v59 * v_v59 / 2.0                # MeV
m_W_v59_gev = m_W_v59 / 1000.0

print(f"Predicted v   = {D_LEPTON}^2 · a^2     = {v_v59:.2f} MeV = {v_v59/1000:.4f} GeV")
print(f"Observed  v                            = {V_HIGGS_MEV:.2f} MeV = {V_HIGGS_GEV:.4f} GeV")
print(f"  gap (v)                              = {abs(v_v59-V_HIGGS_MEV)/V_HIGGS_MEV*100:.4f}%")
print()
print(f"Predicted g_W = √(5·√α)               = {g_W_v59:.5f}")
print(f"Observed  g_W(M_Z)                    ≈ 0.6517")
print(f"  gap (g_W)                            = {abs(g_W_v59-0.6517)/0.6517*100:.3f}%")
print()
print(f"Predicted m_W = (1/2)·g_W·(28^2·a^2)   = {m_W_v59_gev:.4f} GeV")
print(f"Observed  m_W                          = {M_W_GEV:.4f} GeV")
print(f"  gap (m_W)                            = {abs(m_W_v59_gev-M_W_GEV)/M_W_GEV*100:.4f}%")
print()

# ---------------- TEST 4: weak-mixing-angle candidates ----------------
print("-" * 70)
print("TEST 4: weak-mixing-angle candidates from v59 Brannen-t² values")
print("-" * 70)
print(f"sin^2(θ_W)|_onshell  = 1 - m_W^2/m_Z^2 = {SIN2_THW_ONSHELL:.5f}")
print(f"sin^2(θ_W)|_MSbar                     = {SIN2_THW_MSBAR:.5f}")
print()
v59_t_sq = {
    "lepton t²    = 1/2 = 0.5":               1/2,
    "d-quark t²   = 3/5 = 0.6":               3/5,
    "u-quark t²   = 7/9 ≈ 0.7778":             7/9,
    "Brannen phase φ_e/3 = 2/9 ≈ 0.2222":      2/9,
    "1 - t²_u      = 1 - 7/9 = 2/9":           1 - 7/9,
    "1 - t²_d      = 1 - 3/5 = 2/5":           1 - 3/5,
    "(1-t²_u)/(t²_u) = (2/9)/(7/9) = 2/7":     2/7,
    "5·√α":                                    5 * math.sqrt(ALPHA),
}
for name, val in v59_t_sq.items():
    gap_on = abs(val - SIN2_THW_ONSHELL) / SIN2_THW_ONSHELL * 100
    gap_ms = abs(val - SIN2_THW_MSBAR) / SIN2_THW_MSBAR * 100
    print(f"  {name:40s} → {val:.5f}  | gap_onshell {gap_on:6.3f}%  gap_MSbar {gap_ms:6.3f}%")
print()

# cos²θ_W check
print(f"cos²(θ_W)|_onshell = m_W^2/m_Z^2 = {1-SIN2_THW_ONSHELL:.5f}")
print("v59 candidates for cos²θ_W:")
v59_costhw = {
    "t²_u = 7/9":  7/9,
    "1 - 2/9":     1 - 2/9,    # = 7/9 (same)
    "1 - 1/(2·5)": 1 - 0.1,
}
for name, val in v59_costhw.items():
    gap = abs(val - (1-SIN2_THW_ONSHELL)) / (1-SIN2_THW_ONSHELL) * 100
    print(f"  {name:30s} → {val:.5f}  | gap {gap:.4f}%")
print()

# ---------------- TEST 5: derived ratios and consistency checks ----------------
print("-" * 70)
print("TEST 5: combined v59 derivation of (m_W, m_Z, m_H)")
print("-" * 70)
# Conjecture set:
#   v       = 28^2 · a²
#   g_W²    = 5 · √α
#   sin²θ_W = 2/9
# Then:
#   cos²θ_W = 7/9
#   m_W^2   = g_W^2 · v^2 / 4
#   m_Z^2   = m_W^2 / cos²θ_W = (9/7)·m_W^2
#   g'^2    = g_W^2 · tan²θ_W = (2/7)·g_W^2

sin2_thw_v59 = 2/9
cos2_thw_v59 = 7/9

m_W_v59 = math.sqrt(g_W_sq_v59) * v_v59 / 2.0
m_Z_v59_sq = m_W_v59**2 / cos2_thw_v59
m_Z_v59 = math.sqrt(m_Z_v59_sq)

print(f"Inputs:  a (lepton) = {A_LEPTON:.4f} √MeV  (empirical)")
print(f"         α          = {ALPHA:.4e}      (empirical, conjecture -ln α+2α=π²/2)")
print()
print(f"Predictions (m_W and m_Z purely from v59 + α + a):")
print(f"  m_W (v59) = {m_W_v59/1000:.4f} GeV  vs  {M_W_GEV:.4f} GeV  | gap {abs(m_W_v59/1000-M_W_GEV)/M_W_GEV*100:.4f}%")
print(f"  m_Z (v59) = {m_Z_v59/1000:.4f} GeV  vs  {M_Z_GEV:.4f} GeV  | gap {abs(m_Z_v59/1000-M_Z_GEV)/M_Z_GEV*100:.4f}%")
print()

# ---------------- TEST 6: Higgs mass candidates ----------------
print("-" * 70)
print("TEST 6: Higgs mass m_H candidates")
print("-" * 70)
# In the |ξ|² potential V = (λ/4)(|ξ|² - 1/2)², the radial mode mass is m_H² = λ.
# This requires identifying λ in v59 structural terms.
# m_H = 125.20 GeV.  m_H² ≈ 15,675 GeV².
# v² = 60,624 GeV².  m_H² / v² = 0.2586.
# Candidates for m_H² / v²:
print(f"m_H / v       = {M_H_GEV/V_HIGGS_GEV:.5f}")
print(f"m_H^2 / v^2   = {(M_H_GEV/V_HIGGS_GEV)**2:.5f}")
print()
mh_over_v_sq = (M_H_GEV/V_HIGGS_GEV)**2
candidates_mh = {
    "1/4":                   0.25,
    "5√α":                   5 * math.sqrt(ALPHA),
    "5α":                    5 * ALPHA,
    "(2/9)·(7/9·8)·1":       (2/9)*(56/9),
    "Mh²/v² = 1 - cos²β?":   None,
    "g_W²/φ_e/2 = 5√α/(2/9)/2": 5*math.sqrt(ALPHA)/(2/9)/2,
}
for name, val in candidates_mh.items():
    if val is None: continue
    gap = abs(val - mh_over_v_sq) / mh_over_v_sq * 100
    print(f"  {name:35s} → {val:.5f}  | gap {gap:.3f}%")
print()
# Also check: m_H² / a^4
mh_sq_mev = (M_H_GEV*1000.0)**2
print(f"m_H² / a⁴     = {mh_sq_mev / A_SQ**2:.4e}")
print(f"  cf. 28^4 ·5√α/4 = {28**4 * 5*math.sqrt(ALPHA)/4:.4e}")
print()

# ---------------- TEST 7: Yukawa couplings ----------------
print("-" * 70)
print("TEST 7: Yukawa couplings y_l = m_l · √2 / v")
print("-" * 70)
# y_e v / √2 = m_e
# In v59 with v = 28²·a², y_l = m_l √2 / (28²·a²) = (m_l/a²) · √2/784
print(f"v59: y_l = (m_l/a²) · √2 / D_L²  where D_L = 28")
for name, m in [("e", M_E), ("μ", M_MU), ("τ", M_TAU)]:
    y_l_sm = m * math.sqrt(2) / V_HIGGS_MEV
    y_l_v59 = (m / A_SQ) * math.sqrt(2) / (28**2)
    s_sq = m / A_SQ
    s = math.sqrt(s_sq)
    print(f"  y_{name:s} = {y_l_sm:.4e} (SM)  vs {y_l_v59:.4e} (v59 with v=28²a²)  | s²=m_l/a²={s_sq:.4f}")
print()

# ---------------- Save results ----------------
results = {
    "inputs": {
        "M_e_MeV": M_E, "M_mu_MeV": M_MU, "M_tau_MeV": M_TAU,
        "v_Higgs_GeV": V_HIGGS_GEV, "m_W_GeV": M_W_GEV,
        "m_Z_GeV": M_Z_GEV, "m_H_GeV": M_H_GEV,
        "alpha_inv": ALPHA_INV,
        "sin2_thw_onshell": SIN2_THW_ONSHELL,
        "sin2_thw_MSbar": SIN2_THW_MSBAR,
    },
    "derived": {
        "a_lepton_sqrtMeV": A_LEPTON,
        "a_sq_MeV": A_SQ,
        "v_over_a_sq": V_HIGGS_MEV / A_SQ,
    },
    "test1_v_over_a_sq_vs_N_squared": {
        "target": V_HIGGS_MEV / A_SQ,
        "candidates": [{"N": c[1], "N_sq": c[2], "gap_pct": c[3]} for c in candidates_1[:5]],
    },
    "test3_m_W_prediction": {
        "v_v59_GeV":   v_v59/1000.0,
        "v_obs_GeV":   V_HIGGS_GEV,
        "v_gap_pct":   abs(v_v59-V_HIGGS_MEV)/V_HIGGS_MEV*100,
        "g_W_v59":     g_W_v59,
        "m_W_v59_GeV": m_W_v59/1000.0,
        "m_W_obs_GeV": M_W_GEV,
        "m_W_gap_pct": abs(m_W_v59/1000-M_W_GEV)/M_W_GEV*100,
    },
    "test4_sin2_thw_2_over_9": {
        "v59_sin2_thw":      2.0/9.0,
        "obs_sin2_thw_onshell": SIN2_THW_ONSHELL,
        "obs_sin2_thw_MSbar":   SIN2_THW_MSBAR,
        "gap_pct_onshell":   abs(2.0/9.0 - SIN2_THW_ONSHELL)/SIN2_THW_ONSHELL*100,
        "gap_pct_MSbar":     abs(2.0/9.0 - SIN2_THW_MSBAR)/SIN2_THW_MSBAR*100,
    },
    "test5_m_Z_prediction": {
        "m_Z_v59_GeV": m_Z_v59/1000.0,
        "m_Z_obs_GeV": M_Z_GEV,
        "m_Z_gap_pct": abs(m_Z_v59/1000-M_Z_GEV)/M_Z_GEV*100,
    },
}

out_path = os.path.join(os.path.dirname(__file__), "scale_bridge.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"Results saved to {out_path}")
