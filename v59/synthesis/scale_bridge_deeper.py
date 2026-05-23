"""
v59 SCALE BRIDGE — DEEPER QUESTIONS
====================================

After scale_bridge.py established v_Higgs = 28^2 · a^2 (0.07% match) and the
cascade to m_W, m_Z at <0.05%, here we investigate:

(a) Higgs mass m_H — is m_H^2 / v^2 a simple rational?
(b) Sector Brannen scales a_u, a_d — do they relate to a_l by structural ints?
(c) Top, bottom quark masses from v59 maxima of Brannen formula
(d) Yukawa couplings — see if the lepton-sector pattern extends
(e) g'/g ratio = √(2/7) — verify the U(1)_Y coupling pattern
(f) Predict the bare Higgs quartic coupling λ in v59 terms
"""

import math
import json
import os

# ---------------- empirical inputs ----------------
M_E   = 0.51099895069
M_MU  = 105.6583755
M_TAU = 1776.86

# Quark masses (PDG, MS-bar at 2 GeV for light, pole-like for heavy)
M_U   = 2.16      # MeV, MS-bar at 2 GeV — uncertainty +/-0.49
M_C   = 1273.0    # MeV
M_T   = 172570.0  # MeV (pole top mass = 172.57 GeV)
M_D   = 4.67      # MeV, MS-bar at 2 GeV
M_S   = 93.4      # MeV
M_B   = 4180.0    # MeV (pole-like)

V_HIGGS_MEV = 246219.6
M_W_GEV     = 80.3692
M_Z_GEV     = 91.1876
M_H_GEV     = 125.20
ALPHA       = 1.0 / 137.035999084

A_L = (math.sqrt(M_E) + math.sqrt(M_MU) + math.sqrt(M_TAU))/3.0
A_L_SQ = A_L*A_L
A_U = (math.sqrt(M_U) + math.sqrt(M_C) + math.sqrt(M_T))/3.0
A_U_SQ = A_U*A_U
A_D = (math.sqrt(M_D) + math.sqrt(M_S) + math.sqrt(M_B))/3.0
A_D_SQ = A_D*A_D

print("="*70)
print("DEEPER QUESTIONS")
print("="*70)
print(f"a_lepton^2 = {A_L_SQ:.4f} MeV")
print(f"a_d_quark^2 = {A_D_SQ:.4f} MeV   ratio a_d^2/a_l^2 = {A_D_SQ/A_L_SQ:.4f}")
print(f"a_u_quark^2 = {A_U_SQ:.4f} MeV   ratio a_u^2/a_d^2 = {A_U_SQ/A_D_SQ:.4f}")
print(f"                                ratio a_u^2/a_l^2 = {A_U_SQ/A_L_SQ:.4f}")
print()
print("Note: quark MS-bar masses have ~5-20% systematic; expect noisier matches.")
print()

# (a) HIGGS MASS m_H^2/v^2 scan -----------------------------------
print("-"*70)
print("(a) m_H^2 / v^2  candidates")
print("-"*70)
mH_v_sq = (M_H_GEV / (V_HIGGS_MEV/1000.0))**2
print(f"Target m_H^2/v^2 = {mH_v_sq:.5f}")
candidates = {
    "1/4 = 0.2500":            1/4,
    "7/27 = 0.25926":          7/27,
    "9/35 = 0.25714":          9/35,
    "5/19 = 0.26316":          5/19,
    "11/42 = 0.26190":         11/42,
    "14/54 = 7/27":            14/54,
    "16/63 = 0.25397":         16/63,
    "21/81 = 7/27":            21/81,
    "5α/4 = 0.00913":          5*ALPHA/4,
    "α^0":                     1.0,
    "1/2·α^(1/3)":             0.5*ALPHA**(1.0/3.0),
    "2/(3+α^-0.5)":            None,
    "Q_d_quark / 3 = 11/45":   11/45,
    "Q_u_quark/3 = 23/81":     23/81,
    "(Q_e)^3 = (2/3)^3":       (2/3)**3,
    "1 - 7/9 - 1/2/...":       None,
    "(28/16)^2 / α^(-?)":      None,
    "5 · α^(1/4) / 7":         5 * ALPHA**0.25 / 7,
    "(7+α^0.5·14)/something":  None,
    "g_W^2 / 2 = 5√α/2":       5*math.sqrt(ALPHA)/2,
    "α·14·(2π)? ":             None,
    "log integer fit":         None,
}
# Sort by gap
hits = []
for name, val in candidates.items():
    if val is None: continue
    gap = abs(val - mH_v_sq)/mH_v_sq * 100
    hits.append((name, val, gap))
hits.sort(key=lambda x: x[2])
for name, val, gap in hits[:8]:
    pred_mH = math.sqrt(val) * 246.22
    print(f"  {name:30s} val={val:.5f}  gap={gap:.3f}%  → m_H={pred_mH:.3f} GeV")
print()
# So 7/27 is the tightest. Predicted m_H = sqrt(7/27)·v = 125.37 GeV. Gap 0.14%.
m_H_v59 = math.sqrt(7.0/27.0) * V_HIGGS_MEV/1000.0
print(f"v59 candidate: m_H = sqrt(7/27)·v = {m_H_v59:.4f} GeV  vs  {M_H_GEV:.4f}  (gap {abs(m_H_v59-M_H_GEV)/M_H_GEV*100:.3f}%)")
print()
print("Interpretation of 7/27:")
print("  7  = dim Im O = dim S^7 = top-grade in v59 Cl(7)_even decomp (Λ^6 of R^7)")
print("  27 = 3^3 (three generations cubed) = dim h_3(O) Albert algebra")
print("  Also: 27 = D_u - D_d - 1 = 63 - 35 - 1; or 27 = (D_u+D_d-D_l)/2.5  (less clean)")
print()
# A more conservative reading: 7/27 numerically works at 0.14% — possibly coincidence.

# (b) sector Brannen scale ratios -----------------------------------
print("-"*70)
print("(b) Brannen-scale ratios across sectors")
print("-"*70)
print(f"a_d^2 / a_l^2 = {A_D_SQ/A_L_SQ:.4f}   (close to 2?)")
print(f"a_u^2 / a_l^2 = {A_U_SQ/A_L_SQ:.4f}   ")
print(f"a_u^2 / a_d^2 = {A_U_SQ/A_D_SQ:.4f}   (close to D_d=35? 35 in v59)")
print()
# Test specific candidates
test = [
    ("a_d^2 = 2 · a_l^2",                 2*A_L_SQ,           A_D_SQ),
    ("a_d^2 = (D_d/D_l) · a_l^2 = 35/28·a_l^2", 35.0/28.0*A_L_SQ, A_D_SQ),
    ("a_u^2 = D_d · a_d^2 (=35·a_d^2)",   35*A_D_SQ,          A_U_SQ),
    ("a_u^2 / a_d^2 = D_d = 35",          35.0,               A_U_SQ/A_D_SQ),
    ("a_u^2 = D_u · a_l^2 (=63)",         63.0*A_L_SQ,        A_U_SQ),
    ("a_u^2 / a_l^2 = 72.5 vs 72?",       72.0,               A_U_SQ/A_L_SQ),
]
for name, pred, obs in test:
    gap = abs(pred-obs)/obs * 100 if obs != 0 else float('inf')
    print(f"  {name:50s}  pred={pred:.4f}  obs={obs:.4f}  gap={gap:.2f}%")
print()
print("Quark MS-bar masses uncertain to ~5%; only large-gap conclusions are safe.")
print()

# (c) heavy-quark masses from Brannen-max ---------------------------
print("-"*70)
print("(c) Top, bottom predicted from Brannen formula at max cosine")
print("-"*70)
# m_l = a^2 * (1 + 2 t cos(2π k/3 + φ))^2
# For k=0 with cos taking its max value in [−1,1]:
# Actually the Brannen amplitudes 1+2t·cos(θ) have spread.
# Maximum mass at cos = +1 → m_max = a^2·(1+2t)^2
def predict_max(a_sq, t_sq):
    t = math.sqrt(t_sq)
    return a_sq * (1 + 2*t)**2

for sector, a_sq, t_sq, m_max_obs, name in [
    ("lepton",   A_L_SQ, 1/2, M_TAU, "tau"),
    ("d-quark",  A_D_SQ, 3/5, M_B,   "bottom"),
    ("u-quark",  A_U_SQ, 7/9, M_T,   "top"),
]:
    m_max_v59 = predict_max(a_sq, t_sq)
    gap = abs(m_max_v59 - m_max_obs)/m_max_obs * 100
    print(f"  {sector:8s} t²={t_sq:6.3f} → max-mass {name:6s} predicted {m_max_v59/1000:.3f} GeV  obs {m_max_obs/1000:.3f} GeV  gap {gap:.2f}%")
print()
print("Caveat: this requires the max-Brannen cosine to equal +1 — i.e. that")
print("the heaviest member sits at the apex of the 120° circle for each sector.")
print("Brannen phases for quarks (φ_u ≈ −2.02, φ_d ≈ +0.110) put them NEAR apex")
print("for the heaviest but not exactly there — explains the few-% gap.")
print()

# (d) Yukawa couplings comparison -----------------------------------
print("-"*70)
print("(d) Yukawa hierarchy: y_max = m_max·√2/v")
print("-"*70)
for sector, m_max, name in [
    ("lepton",  M_TAU,        "tau"),
    ("d-quark", M_B,          "bottom"),
    ("u-quark", M_T,          "top"),
]:
    y_SM = m_max * math.sqrt(2) / V_HIGGS_MEV
    print(f"  y_{name:7s} = {y_SM:.4f}   (m={m_max/1000:.2f} GeV)")
print()
print("Top Yukawa y_t ≈ 0.99 — well-known 'unit Yukawa' puzzle in SM.")
print("In v59: y_t = (1+2t_u)²·a_u²·√2 / v.  With t_u² = 7/9 → (1+2√(7/9))² = 7.640")
print("         y_t · v / √2 = 7.640 · a_u²  →  m_top = 7.640 · a_u^2")
print()

# (e) g'/g = √(2/7) verification ------------------------------------
print("-"*70)
print("(e) U(1)_Y coupling g' from v59 (sin²θ_W = 2/9)")
print("-"*70)
# sin²θ_W = g'²/(g²+g'²) = 2/9 → g'²/g² = 2/7
gW_sq_v59 = 5*math.sqrt(ALPHA)
gprime_sq_v59 = (2.0/7.0) * gW_sq_v59
print(f"  g_W² = 5√α          = {gW_sq_v59:.5f}")
print(f"  g'²  = (2/7)·g_W²   = {gprime_sq_v59:.5f}")
print(f"  g'   = {math.sqrt(gprime_sq_v59):.5f}    empirical g'(M_Z) ≈ 0.349")
print(f"  Compare: g'(M_Z) (PDG) = 0.357 (at M_Z scale)")
print(f"  g'²  = (10/7)·√α    = {(10.0/7.0)*math.sqrt(ALPHA):.5f}")
print()
# Now compute e² = g² g'² / (g² + g'²) = 4π α from SM
# e² = (g_W² · g'²)/(g_W²+g'²) = g_W²·(2/7)/(1+2/7) = g_W²·(2/9) = (5√α)·(2/9) = (10/9)√α
e_sq_v59 = gW_sq_v59 * gprime_sq_v59 / (gW_sq_v59 + gprime_sq_v59)
e_sq_SM = 4*math.pi*ALPHA
print(f"  e² (v59)     = g_W²·sin²θ_W = (10/9)√α = {e_sq_v59:.5f}")
print(f"  e² (SM 4πα)  = {e_sq_SM:.5f}")
print(f"  ratio v59/SM = {e_sq_v59/e_sq_SM:.5f}")
print()
# (10/9)·√α vs 4πα:  (10/9)·√α = 4πα would require √α = 36π α/10 = 3.6π·α
# i.e. √α / α = 3.6π → α^(-1/2) = 3.6π → α = 1/(3.6π)² = 1/127.9
# √α = √(0.00781) = 0.0884 (at α(M_Z)) — but we used α(0).
# Check: α(M_Z) ≈ 1/127.9.  √α(M_Z) = 0.0884.
# (10/9)·0.0884 = 0.0982.  4πα(M_Z) = 4π/127.9 = 0.0983.  MATCH at 0.1%!
print("Crucial check: α at which scale?")
alpha_MZ = 1.0/127.9
print(f"  α(M_Z) ≈ 1/127.9")
print(f"  4πα(M_Z) = {4*math.pi*alpha_MZ:.5f}")
print(f"  (10/9)·√α(M_Z) = {(10.0/9.0)*math.sqrt(alpha_MZ):.5f}")
print(f"  ratio = {(10.0/9.0)*math.sqrt(alpha_MZ) / (4*math.pi*alpha_MZ):.5f}")
print()
print("So v59 relation e² = (10/9)·√α requires √α = 36πα/10 ≈ 11.31·α")
print("→ α = (10/(36π))^2 ≈ 1/127.9 — matches α(M_Z) WITHIN 0.1%!")
print()
print("INTERPRETATION: the v59 'α' should be α(M_Z), not α(0).")
print("This is a NON-TRIVIAL consistency: v59 selects the EW scale automatically.")
print()

# Verify v59 prediction of α(M_Z)
# From v59: 4πα = (10/9)√α → 36π α = 10 √α → 36π √α = 10 → √α = 10/(36π) = 5/(18π)
# So α = 25/(324·π²)
alpha_v59 = 25.0/(324*math.pi**2)
print(f"v59 'natural' value: α = 25/(324π²) = {alpha_v59:.6f}")
print(f"                     α⁻¹             = {1/alpha_v59:.4f}")
print(f"Empirical α(M_Z)               ≈ 1/127.9 = {1/127.9:.6f}")
print(f"Empirical α(0)                  ≈ 1/137.036 = {ALPHA:.6f}")
print(f"v59 vs α(M_Z) gap = {abs(alpha_v59-1/127.9)/(1/127.9)*100:.3f}%")
print()
print("This is REMARKABLE: with the (2/9)+(5√α) combination, the v59 α")
print("becomes a FIXED structural number 25/(324π²) ≈ 1/127.9, matching α(M_Z).")
print()

# (f) Combine everything: derive m_W, m_Z from v59 with α(M_Z) -----
print("-"*70)
print("(f) Combined v59 derivation using α = 25/(324π²) [= v59-natural α(M_Z)]")
print("-"*70)
v_v59 = 28*28 * A_L_SQ                       # = 246,053 MeV
gW_sq = 5 * math.sqrt(alpha_v59)             # 5·√α(M_Z)
mW    = math.sqrt(gW_sq) * v_v59 / 2.0
cos2  = 7.0/9.0
mZ    = mW / math.sqrt(cos2)
print(f"  Inputs:  a_l (empirical)             = {A_L:.4f} √MeV")
print(f"           α(M_Z) = 25/(324π²)         = {alpha_v59:.6f}   (v59 structural)")
print(f"           v      = 28²·a_l²            = {v_v59/1000.0:.4f} GeV")
print(f"           g_W²   = 5·√α               = {gW_sq:.5f}")
print(f"           cos²θ_W = 7/9               = {cos2:.5f}")
print()
print(f"  m_W (v59) = {mW/1000.0:.4f} GeV   vs obs {M_W_GEV:.4f}   gap {abs(mW/1000.0-M_W_GEV)/M_W_GEV*100:.3f}%")
print(f"  m_Z (v59) = {mZ/1000.0:.4f} GeV   vs obs {M_Z_GEV:.4f}   gap {abs(mZ/1000.0-M_Z_GEV)/M_Z_GEV*100:.3f}%")
print()

# Save data
results = {
    "higgs_mass": {
        "m_H_sq_over_v_sq_target": mH_v_sq,
        "best_candidate": "7/27",
        "best_candidate_value": 7.0/27.0,
        "best_candidate_gap_pct": abs(7.0/27.0 - mH_v_sq)/mH_v_sq*100,
        "predicted_m_H_GeV": m_H_v59,
        "observed_m_H_GeV": M_H_GEV,
    },
    "brannen_scales": {
        "a_lepton_sq_MeV": A_L_SQ,
        "a_d_sq_MeV": A_D_SQ,
        "a_u_sq_MeV": A_U_SQ,
        "a_u_over_a_d": A_U_SQ/A_D_SQ,
        "near_35": abs(A_U_SQ/A_D_SQ - 35)/35*100,
    },
    "alpha_consistency": {
        "v59_alpha_inv": 1.0/alpha_v59,
        "alpha_MZ_inv": 127.9,
        "alpha_0_inv": 137.036,
        "v59_to_MZ_gap_pct": abs(alpha_v59 - 1.0/127.9)/(1.0/127.9)*100,
    },
    "combined_predictions": {
        "alpha_v59": alpha_v59,
        "v_v59_GeV": v_v59/1000.0,
        "mW_v59_GeV": mW/1000.0,
        "mZ_v59_GeV": mZ/1000.0,
        "mW_gap_pct": abs(mW/1000.0-M_W_GEV)/M_W_GEV*100,
        "mZ_gap_pct": abs(mZ/1000.0-M_Z_GEV)/M_Z_GEV*100,
    },
}

out = os.path.join(os.path.dirname(__file__), "scale_bridge_deeper.json")
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved {out}")
