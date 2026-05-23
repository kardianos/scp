"""
HIGGS QUARTIC IDENTITY TEST
============================

Reformulation of m_H = √(7/27)·v:

  m_H^2 / v^2 = 7/27 = (7/9) / 3 = cos²θ_W / generations

Equivalently:

  m_H^2 · m_Z^2 = m_W^2 · v^2 / 3
  λ_Higgs (in V = -μ²|Φ|² + λ|Φ|⁴) = cos²θ_W / (2·generations) = (7/9)/6 = 7/54

These are CONSEQUENCES of (m_H/v)² = 7/27 plus the v59 relations
v = 28²·a², m_W = (g_W/2)·v, m_Z/m_W = 3/√7 (i.e., cos²θ_W = 7/9).

Test all four equivalent forms numerically.
"""
import math

M_W = 80.3692
M_Z = 91.1876
M_H = 125.20
V = 246.2196

# Form 1: (m_H/v)² = 7/27
lhs1 = (M_H/V)**2
rhs1 = 7.0/27.0
print(f"Form 1:  m_H²/v² = {lhs1:.5f}    7/27 = {rhs1:.5f}    gap {abs(lhs1-rhs1)/rhs1*100:.3f}%")

# Form 2: m_H² · m_Z² = m_W² · v² / 3
lhs2 = M_H**2 * M_Z**2
rhs2 = M_W**2 * V**2 / 3.0
print(f"Form 2:  m_H²·m_Z² = {lhs2:.4e}    m_W²·v²/3 = {rhs2:.4e}    gap {abs(lhs2-rhs2)/rhs2*100:.3f}%")

# Form 3: λ_Higgs = m_H²/(2v²) = cos²θ_W / (2·gen)
# With cos²θ_W = m_W²/m_Z² (on-shell) and gen=3:
lambda_lhs = M_H**2 / (2 * V**2)
cos2_thw_obs = M_W**2 / M_Z**2
lambda_v59 = cos2_thw_obs / (2*3)
lambda_pure_v59 = (7.0/9.0) / 6.0    # using cos²θ_W = 7/9
print(f"Form 3:  λ_SM = {lambda_lhs:.5f}    cos²θ_W/6 (obs) = {lambda_v59:.5f}    pure v59 (7/54) = {lambda_pure_v59:.5f}")
print(f"         gap (SM vs pure v59) = {abs(lambda_lhs-lambda_pure_v59)/lambda_lhs*100:.3f}%")

# Form 4: m_H² · m_Z² · gen = m_W² · v²
print()
print("Equivalent SM-quadratic identities:")
print(f"  m_H · m_Z · √gen = m_W · v ?    {M_H*M_Z*math.sqrt(3):.4f}  vs  {M_W*V:.4f}    gap {abs(M_H*M_Z*math.sqrt(3) - M_W*V)/M_W/V*100:.3f}%")
print(f"  m_H²·m_Z²·3 = m_W²·v² ?         {M_H**2*M_Z**2*3:.4e}  vs  {M_W**2*V**2:.4e}    gap {abs(M_H**2*M_Z**2*3 - M_W**2*V**2)/(M_W**2*V**2)*100:.3f}%")
print()

print("Cleanest STRUCTURAL form:")
print("  m_H^2 · m_Z^2 · (n_gen) = m_W^2 · v^2 ")
print(f"  where n_gen = 3 (number of fermion generations)")
print(f"  Match: {abs(M_H**2*M_Z**2*3 - M_W**2*V**2)/(M_W**2*V**2)*100:.3f}%")
print()
print("Or in terms of the SM Higgs quartic:")
print(f"  λ_Higgs = cos²θ_W / (2·n_gen) = 7/54 ≈ {7/54:.5f}")
print(f"  Observed: λ = m_H²/(2v²) = {lambda_lhs:.5f}")
print(f"  Gap: {abs(7/54 - lambda_lhs)/lambda_lhs*100:.3f}%")
print()

# Now do the same with sin²θ_W instead
print("Alternative: sin²θ_W variants")
# m_H²/v² = 1/3 - 5sin⁴θ_W ? Let's check: 1/3 - (1/3)·sin²θ_W·something
sin2 = 2.0/9.0
print(f"  sin²θ_W = 2/9, sin⁴ = 4/81 = {sin2**2:.5f}")
print(f"  m_H²/v² = (1-sin²θ_W)/3 = (1-2/9)/3 = 7/27 = {(1-sin2)/3:.5f}")
print(f"  Observed: {lhs1:.5f}")
print()
print("So: m_H² = (1 - sin²θ_W)·v² / n_gen")
print("    With sin²θ_W = 2/9 (Brannen phase) and n_gen = 3:")
print(f"    m_H² = (7/9)·v²/3 = (7/27)·v²  →  m_H = √(7/27)·v = {math.sqrt(7/27)*V:.4f} GeV")
print(f"    vs observed m_H = {M_H} GeV.  Gap {abs(math.sqrt(7/27)*V - M_H)/M_H*100:.3f}%")
