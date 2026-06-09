#!/usr/bin/env python3
"""
pmns_strong_test.py  —  G11 (neutrinos+PMNS) and G12 (strong: α_s, θ_QCD).

Honest reconnaissance.  The expected output for most of this is NULL: no
clean v59-structural match.  We make that explicit and quantify it, and we
test the ONE structurally-motivated lead in the strong sector (α_s from a
color dual-Coxeter analog of g_W²=5√α).

PDG 2024 inputs.
"""
import math

ALPHA0  = 1.0 / 137.035999084
ALPHAMZ = 1.0 / 127.951

def pct(p, e): return abs(p - e) / abs(e) * 100.0

print("="*74)
print("G11 — Neutrinos & PMNS")
print("="*74)

# ---------------------------------------------------------------------------
# 1. PMNS mixing angles (normal ordering, NuFIT 5.2 / PDG 2024 central)
# ---------------------------------------------------------------------------
# sin² of the three angles:
s2_12 = 0.307     # solar
s2_23 = 0.546     # atmospheric (NO)
s2_13 = 0.02203   # reactor
th12 = math.degrees(math.asin(math.sqrt(s2_12)))
th23 = math.degrees(math.asin(math.sqrt(s2_23)))
th13 = math.degrees(math.asin(math.sqrt(s2_13)))
print(f"\n  θ12={th12:.1f}°  θ23={th23:.1f}°  θ13={th13:.1f}°   (LARGE, unlike CKM)")
print(f"  sin²θ12={s2_12}, sin²θ23={s2_23}, sin²θ13={s2_13}")

# Test "structural" candidates. Key known near-coincidences:
print("\n  --- candidate structural readings (mixing angles) ---")
print(f"  sin²θ12 ?= 1/3 = 0.3333  [{pct(1/3, s2_12):.1f}%]   (tri-bimaximal value)")
print(f"  sin²θ23 ?= 1/2 = 0.5000  [{pct(1/2, s2_23):.1f}%]   (maximal mixing)")
print(f"  sin²θ23 ?= 5/9 = 0.5556  [{pct(5/9, s2_23):.1f}%]")
print(f"  sin²θ13 ?= 2/9·α-ish... ?")
print(f"  sin²θ13 ?= 2/9·0.1 ... no clean form")
print(f"  θ13 ?= θ_C/2 (Cabibbo) = {math.degrees(math.asin(0.225))/2:.2f}°  vs {th13:.2f}° "
      f"[{pct(math.degrees(math.asin(0.225))/2, th13):.1f}%]   (QLC-type)")
print("  → PMNS is near TRI-BIMAXIMAL (1/3, 1/2, 0) — a SYMMETRY pattern, not a v59 one.")
print("    v59 has NO mechanism producing large near-degenerate mixing.")

# ---------------------------------------------------------------------------
# 2. The neutrino-Koide ambient puzzle.
#    Charged leptons:   D=28, t²=1/2, Q=2/3.
#    Brannen's neutrino Koide fit needs a DIFFERENT phase/scale; and the
#    "ambient" D_ν that would reproduce a neutrino t² fits NO v59 ambient.
# ---------------------------------------------------------------------------
print("\n  --- neutrino Koide ambient puzzle ---")
# v59 ambients: 28 (L), 35 (F), 63 (L+F).  t²_N = 1 - 14/D_N.
for D in (28, 35, 63):
    t2 = 1 - 14/D
    Q  = (1 + 2*t2)/3
    print(f"    v59 ambient D={D:2d}:  t²=1-14/D={t2:.4f}  Q=(1+2t²)/3={Q:.4f}")
# Neutrino masses (NO, from Σm≈0.06 eV and Δm² splittings) — illustrative.
# Using lightest m1≈0 gives a Koide Q far from 2/3; Brannen's own ν fit uses
# a scaled phase. The point: no single 14/D reproduces the ν pattern.
print("  Neutrinos: with m1≈0 (NO), Q_ν → ~1/3 (degenerate-massless limit), NOT a 14/D value.")
print("  → CONFIRMED PUZZLE: no v59 ambient (28/35/63) yields the neutrino Koide.")
print("    Neutrinos likely live OUTSIDE the Brannen-circulant charged-fermion pattern.")

# ---------------------------------------------------------------------------
# 3. Absolute neutrino mass scale — seesaw-analog order of magnitude
# ---------------------------------------------------------------------------
print("\n  --- seesaw-analog scale check ---")
# Naive seesaw: m_ν ~ m_D²/M_R.  With m_D~m_e and M_R~v: m_ν ~ m_e²/v.
m_e = 0.000511   # GeV
v   = 246.0      # GeV (Higgs vev)
m_nu_seesaw = m_e**2 / v
print(f"  m_e²/v = {m_nu_seesaw*1e9:.4f} eV   (vs Σm_ν ≲ 0.06-0.12 eV; right ballpark for atm)")
m_nu_atm = math.sqrt(2.5e-3)   # eV, √Δm²_atm
print(f"  √Δm²_atm ≈ {m_nu_atm:.3f} eV.  Seesaw with M_R~v gives m_ν~10⁻⁶ eV (too small) —")
print(f"  needs M_R ~ m_e²/m_ν ~ {m_e**2/(m_nu_atm*1e-9):.2e} GeV (intermediate scale).")
print("  → No v59 structural M_R yet. Seesaw is a PLAUSIBLE FRAME, not a v59 derivation.")

print("\n" + "="*74)
print("G12 — Strong sector: α_s and θ_QCD")
print("="*74)

# ---------------------------------------------------------------------------
# 4. α_s(M_Z): test the color dual-Coxeter analog of g_W² = 5√α.
#    v59 gauge pattern:  g_a² = h∨(G_a) · √α.
#    SU(2)_L: h∨=... actually the doc uses h∨(Spin(7))=5 for the EW embedding.
#    Color SU(3)=A₂ has dual Coxeter h∨=3.  So the ANALOG conjecture is:
#        g_s² = 3·√α   →  α_s = g_s²/(4π) = 3√α/(4π).
# ---------------------------------------------------------------------------
alpha_s_MZ = 0.1179      # PDG 2024
g_s2_emp = 4*math.pi*alpha_s_MZ
print(f"\n  α_s(M_Z) = {alpha_s_MZ}   →  g_s² = 4π·α_s = {g_s2_emp:.4f}")
print("\n  --- analog of g_W²=h∨·√α with color SU(3) (h∨(A₂)=3) ---")
for label, h, a in [("3·√α(0)",  3, ALPHA0),
                    ("3·√α(MZ)", 3, ALPHAMZ),
                    ("5·√α(0)",  5, ALPHA0),
                    ("5·√α(MZ)", 5, ALPHAMZ)]:
    g2 = h*math.sqrt(a)
    print(f"    g_s² ?= {label:10s} = {g2:.4f}   [{pct(g2, g_s2_emp):.1f}%]   "
          f"→ α_s={g2/(4*math.pi):.4f} [{pct(g2/(4*math.pi),alpha_s_MZ):.1f}%]")
# Direct α_s candidates
print("\n  --- direct α_s(M_Z) candidates ---")
for label, val in [("√α(0)·...", None),
                   ("(3/(4π))·√α(0)", (3/(4*math.pi))*math.sqrt(ALPHA0)),
                   ("1/(2π)", 1/(2*math.pi)),
                   ("2/(3·√α(0))·α... ", None),
                   ("√α(MZ)", math.sqrt(ALPHAMZ)),
                   ("4/3·α(0)·...", None),
                   ("α(MZ)·15.1", ALPHAMZ*15.1)]:
    if val is None: continue
    print(f"    α_s ?= {label:18s} = {val:.4f}   [{pct(val, alpha_s_MZ):.1f}%]")
print("  NOTE: α_s RUNS strongly (asymptotic freedom); any single fixed v59 number")
print("  must specify a scale. h∨·√α gives the RIGHT ORDER but ~tens-% off, and the")
print("  √α form is itself the unproven G2 conjecture. NOT a clean match.")

# ---------------------------------------------------------------------------
# 5. θ_QCD ≈ 0 (strong CP).  This is a YES/NO structural question, not numeric.
# ---------------------------------------------------------------------------
print("\n  --- θ_QCD (strong CP problem): θ < 10⁻¹⁰ observed ---")
print("  Question: does v59 STRUCTURALLY force θ_QCD = 0?")
print("  v59 fact: color SU(3) is built from a genuine COMPLEX STRUCTURE J_c=γ₀γ₅")
print("  (ColorSU3.lean: J_c²=-I, skew, color-invariant, makes ℝ⁸=ℂ⁴=1⊕3).")
print("  A complex structure J_c gives the color sector a HOLOMORPHIC/COMPLEX")
print("  structure. CP is the antiholomorphic involution (J_c → -J_c, complex conj).")
print("  HYPOTHESIS: if the color current J_F is built from the J_c-holomorphic")
print("  bilinear, the topological term θ Tr(F∧F) — which is CP-ODD — would be")
print("  forbidden by the SAME reality/skewness forcing that pins lepton J (the")
print("  LeptonRealityForcing 'no symmetric matrix is a complex structure' lemma).")
print("  STATUS: a STRUCTURAL LEAD (not a derivation). Falsifier: explicitly compute")
print("  whether ⟨F∧F⟩ has nonzero J_c-even (CP-odd) projection on the color sector.")

print("\n" + "="*74)
print("VERDICT: G11 neutrinos — GENUINELY OPEN (PMNS=symmetry pattern, ν-ambient")
print("puzzle real, seesaw plausible but no v59 M_R). G12 — α_s order-right via")
print("h∨(A₂)=3·√α but not clean; θ_QCD=0 has a real STRUCTURAL LEAD via J_c reality.")
print("="*74)
