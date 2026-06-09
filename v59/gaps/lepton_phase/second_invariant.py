#!/usr/bin/env python3
"""
second_invariant.py — the "second Koide" route: is cos(3φ)=cos(2/3) fixed by a
second symmetric-function constraint on the lepton masses, analogous to how Q=2/3
is fixed by the first?

Logic.  The Brannen kernel has 3 real shape data beyond the scale a:
    p1 = Σ√m = 3a                       (sets a)
    p2 = Σm  = 3a²(1+2t²)               (Koide:  Q=p2/p1² = (1+2t²)/3 fixes t²)
    p3 = Σm^{3/2} = 3a³(1+6t²+2t³cos3φ) (the ONLY phase-dependent moment; fixes cos3φ)
So the phase invariant cos(3φ) is, exactly, the third standardized moment (skewness)
of the √m distribution. A NON-geometric origin for cos(2/3) would be: a structural
relation that fixes this third moment to the value cos(2/3) ≈ 0.78589.

This script:
 (1) confirms cos(3φ) = (normalized third cumulant of √m), expressing the phase
     invariant as a pure data quantity (no angle);
 (2) tests whether the data's third moment equals any clean structural value;
 (3) tests the 'self-similar' hypothesis: that the SECOND invariant is the SAME
     ratio applied one level up — i.e. cos(3φ) = f(Q) for some structural f, and
     whether the observed f is non-geometric/clean.

Run: python3 second_invariant.py
"""
import numpy as np

PDG_LEPTON = [0.51099895e-3, 105.6583755e-3, 1776.86e-3]  # GeV
Q_T, PHI_T = 2/3, 2/9
COS_INV = np.cos(2/3)

def hr(t): print("\n"+"="*74+f"\n{t}\n"+"="*74)

# ---- recover Brannen shape from data ---------------------------------------
sm = np.sqrt(np.array(PDG_LEPTON)); a = sm.mean()
Q  = (sm**2).sum()/sm.sum()**2
t2 = (3*Q-1)/2; t = np.sqrt(t2)
# principal-branch phase from the residuals
k = np.arange(3); ang = 2*np.pi*k/3
C = np.sum((sm/a - 1)*np.cos(ang)); Ss = np.sum((sm/a - 1)*np.sin(ang))
phi = (-np.arctan2(Ss, C)) % (2*np.pi/3)

hr("1. cos(3φ) IS the standardized third moment of the √m distribution")
# p3 = 3a³(1+6t²+2t³cos3φ); solve for cos3φ from the raw √m data:
p1 = sm.sum(); p2 = (sm**2).sum(); p3 = (sm**3).sum()
a_d = p1/3
cos3phi_data = (p3/(3*a_d**3) - 1 - 6*t2) / (2*t2*t)
print(f"  from raw √m moments:  cos(3φ) = {cos3phi_data:.8f}")
print(f"  direct cos(3·φ_fit):           = {np.cos(3*phi):.8f}")
print(f"  target cos(2/3):               = {COS_INV:.8f}")
print(f"  |cos3φ_data − cos(2/3)|        = {abs(cos3phi_data-COS_INV):.2e}")
print("  -> The phase invariant is, literally, the (rescaled) THIRD CUMULANT / skewness")
print("     of the √m values. So 'derive cos(2/3)' = 'fix the skewness of √m structurally'.")
print("     This is the cleanest NON-GEOMETRIC framing: no angle, a pure moment.")

hr("2. The 'two-Koide' picture: which standardized moments are structural?")
# Standardize x_k = (√m_k - a)/a = 2t cos(φ+2πk/3). Its moments:
x = sm/a_d - 1
m2 = np.mean(x**2); m3 = np.mean(x**3)
print(f"  x_k = (√m_k−a)/a:   mean={np.mean(x):.2e} (≈0 by construction)")
print(f"  2nd moment  ⟨x²⟩ = 2t²        : data {m2:.6f}, model {2*t2:.6f}  (Koide: =>Q)")
print(f"  3rd moment  ⟨x³⟩ = (3/2)·2t³cos3φ? check:")
print(f"       data ⟨x³⟩ = {m3:.6f}")
print(f"       model (3/2)·... : per BrannenPhase, Σcos³ = (3/4)cos3φ, so ⟨x³⟩ = (2t)³·(1/4)cos3φ")
model_m3 = (2*t)**3 * 0.25 * np.cos(3*phi)
print(f"       model ⟨x³⟩ = (2t)³·(1/4)cos3φ = {model_m3:.6f}")
print("  -> Koide fixes the 2nd moment (variance) of √m via Q. A 'second Koide' would")
print("     fix the 3rd moment (skewness) -> cos3φ. The DATA says that skewness lands at")
print("     exactly the value making cos3φ = cos(Q) [Q in radians]. THAT coincidence —")
print("     skewness-determined-cos equals cosine-of-the-variance-ratio — is the mystery.")

hr("3. Is cos(3φ)=cos(2/3) explained by a clean f(Q) other than the radian insert?")
# The law is 3φ=Q i.e. cos3φ=cos(Q). Test alternative clean closed forms for cos3φ
# in terms of Q,t² that DON'T require inserting Q as a radian.
candidates = {
    "cos(Q)            [the φ=Q/3 law]": np.cos(Q_T),
    "Q                 (cos3φ=Q?)"      : Q_T,
    "1-Q/3             ":                 1-Q_T/3,
    "(1+2t²)/... no"    :                np.nan,
    "7/9 = cos²θ_W      [ruled out]"    : 7/9,
    "11/14             (B1 numerology)" : 11/14,
    "1 - t²/... "       :                1 - t2/2.3,
    "√(2/3)            ":                np.sqrt(2/3),
}
print(f"  observed cos(3φ) = {np.cos(3*phi):.8f}")
for name,val in candidates.items():
    if np.isnan(val): continue
    print(f"     {name:36s} = {val:.8f}   |Δ|={abs(val-np.cos(3*phi)):.2e}")
print("  -> Only cos(Q) (the radian insert) hits it to 1e-5. Every NON-radian closed form")
print("     (Q itself, 1-Q/3, 7/9, √(2/3), …) MISSES. So the data genuinely prefers the")
print("     transcendental 'cos of Q-as-radian' over any algebraic f(Q). The φ=Q/3 law is")
print("     not disguising a simpler algebraic relation — the radian insert is real.")

hr("VERDICT")
print("""  * cos(3φ) = cos(2/3) is EXACTLY the standardized third moment (skewness) of the
    √m distribution — a pure data invariant, no geometry. [reframing, clean]
  * A 'second Koide' fixing this skewness structurally would close G7. None known.
  * Crucially, the data does NOT prefer any clean ALGEBRAIC f(Q) for cos3φ — only the
    transcendental cos(Q-as-radian) fits at 1e-5. So the magnitude problem cannot be
    dodged by finding cos3φ = (some rational/Lie ratio); it must produce cos(2/3) itself.
  * This SHARPENS G7: the open target is not '2/3' but the specific transcendental
    skewness value cos(2/3), and there is no non-geometric algebraic shortcut to it.""")
