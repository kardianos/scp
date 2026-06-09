#!/usr/bin/env python3
"""
v60/gaps/rank_tension/01_rank_tension.py

G1 — the rank tension of the electroweak bridge v_Higgs = dim(L)²·a_ℓ² = 784 a².

v59 (gaps/ew_scale_bridge/FINDINGS.md §3) sharpened G1 to: ONE vacuum bilinear
Y ∈ End(L) = M₂₈(ℝ) is required to be simultaneously
  (count/R1) full-rank "democratic"  → ‖Y‖²_F = 784 a²  (the EW scale), and
  (spectrum) rank-3 (3 generations)  → ‖Y‖²_F = Σm = 9Q a² = 6 a²  (gravity charge),
which is impossible for one matrix (784 ≠ 6, full-rank ≠ rank-3), and no single-step
so(8)→H breaking leaves "3 light + 25 eaten."

This script establishes, by computation:
  (1) the baseline numbers (bridge 0.068%, Σm=9Qa², ratio 6/784);
  (2) DEFLATION: the advertised "Σm = (6/784)·v ties gravity to EW at 0.07%" is NOT
      independent evidence — since Σm=9Qa² is DEFINITIONAL, that match is exactly the
      EW bridge match v=784a² re-expressed (cf. v59's g_W²=5√α deflation);
  (3) the single-Y rank tension is real: a "democratic" vacuum gives 28 comparable
      singular values (NOT the rank-3 hierarchy), while a hand-built rank-3+heavy Y
      reproduces both norms but is NOT what an so(8)-invariant potential selects and
      its 25 heavy directions have no symmetry home;
  (4) the subalgebra obstruction (max proper subalgebra of so(8) = 21 < 25);
  (5) verdict: the only consistent reading is TWO objects on DIFFERENT spaces
      (End(L)=784 operator algebra vs the 3-dim generation space), which leaves
      v_Higgs a conditional second scale (tied to a_ℓ only via R1+R2, no dynamical home).

Run:  python 01_rank_tension.py
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

# ---- empirical anchors (PDG, MeV) ----
m_lepton = np.array([0.51099895, 105.6583755, 1776.86])
v_obs = 246220.0
N_gen = 3
dimL = 28


def brannen(m):
    s = np.sqrt(m)
    a = s.mean()                 # a = (Σ√m)/N_gen
    Q = m.sum() / s.sum() ** 2   # Koide
    return a, Q


def section(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main():
    a, Q = brannen(m_lepton)
    Sigma_m = m_lepton.sum()

    section("(1) Baseline: the two Frobenius² readings")
    v_bridge = dimL ** 2 * a ** 2
    print(f"  a_ℓ = (Σ√m)/3            = {a:.6f} MeV^½")
    print(f"  Koide Q                  = {Q:.6f}   (≈ 2/3 = {2/3:.6f})")
    print(f"  Σm  = 9 Q a²  (gravity charge, Frobenius² over 3-gen space) = {Sigma_m:.2f} MeV")
    print(f"        9 Q a²             = {9*Q*a**2:.2f} MeV   (identity check)")
    print(f"  v   = 784 a² (EW scale, Frobenius² over End(L)=M₂₈)         = {v_bridge:.1f} MeV")
    print(f"  v_obs                                                       = {v_obs:.1f} MeV")
    print(f"  bridge match  v_bridge/v_obs = {v_bridge/v_obs:.5f}  (gap {abs(v_bridge-v_obs)/v_obs*100:.3f}%)")
    assert np.isclose(9*Q*a**2, Sigma_m, rtol=1e-9)  # Σm = 9Qa² is exact (definitional)

    section("(2) DEFLATION: 'Σm = (6/784) v at 0.07%' is the bridge match in disguise")
    ratio_struct = 9 * Q / dimL**2          # = 6/784 (uses Q)
    ratio_obs = Sigma_m / v_obs
    print(f"  9Q/dim(L)²  = {ratio_struct:.6f}   (= 6/784 = {6/784:.6f} for Q=2/3)")
    print(f"  Σm / v_obs  = {ratio_obs:.6f}")
    print(f"  'match' rel error = {abs(ratio_struct-ratio_obs)/ratio_struct:.2e}")
    print("  But Σm = 9Q a² is DEFINITIONAL (a=(Σ√m)/3, Q=Koide). Therefore")
    print(f"    (Σm/v_obs) / (9Q/784) = 784 a² / v_obs = {dimL**2*a**2/v_obs:.5f}")
    print(f"    = exactly the EW bridge match v_bridge/v_obs = {v_bridge/v_obs:.5f}.")
    same = np.isclose(ratio_obs/ratio_struct, v_bridge/v_obs, rtol=1e-9)
    print(f"  ⇒ the 'gravity↔EW 6/784 bonus' carries NO information beyond v=784a². [{same}]")
    print("    (Analogous to v59's g_W²=5√α = α(M_Z)-in-disguise deflation.)")
    assert same

    section("(3) The single-Y rank tension is real")
    rng = np.random.default_rng(0)
    # (a) 'democratic' (maximally-mixed) vacuum: random ±a entries
    Yd = rng.choice([-1.0, 1.0], size=(28, 28)) * a
    sv_d = np.linalg.svd(Yd, compute_uv=False)
    print("  (a) democratic vacuum Y (random ±a, 'all 784 components comparable'):")
    print(f"      ‖Y‖²_F = {np.sum(Yd**2):.1f} ≈ 784a² = {784*a**2:.1f}  ✓ (sets EW scale)")
    print(f"      singular values: full rank {np.sum(sv_d > 1e-9*sv_d.max())}, all O(√28·a),"
          f" range [{sv_d.min():.1f}, {sv_d.max():.1f}] (no hierarchy)")
    print("      → 28 comparable 'light' directions, NOT the rank-3 hierarchy. Wrong spectrum.")

    # (b) physical Brannen kernel: rank-3 on the 3-gen space
    print("  (b) physical 3×3 Brannen kernel: singular values = √m_k =",
          np.round(np.sqrt(m_lepton), 3))
    print(f"      ‖M‖²_F = Σm = {Sigma_m:.1f} = 6a² = {6*a**2:.1f}  (hierarchical, rank 3)")

    # (c) hand-built two-piece Y: 3 Brannen + 25 heavy, topped to 784a²
    heavy2 = (784 * a**2 - Sigma_m) / 25.0
    Mheavy = np.sqrt(heavy2)
    sv_2piece = np.concatenate([np.sqrt(m_lepton), np.full(25, Mheavy)])
    print(f"  (c) hand-built two-piece Y: 3 light (√m_k) + 25 heavy (σ={Mheavy:.2f} = {Mheavy/a:.2f}a):")
    print(f"      ‖Y‖²_F = {np.sum(sv_2piece**2):.1f} = 784a²  ✓  AND rank-3 light spectrum ✓")
    assert np.isclose(np.sum(sv_2piece**2), 784*a**2, rtol=1e-9)
    print("      BUT: this is a HIERARCHICAL vacuum (3 small + 25 at ~5.6a) — NOT the")
    print("      democratic vacuum an so(8)-invariant potential selects (which gives (a)),")
    print("      and the 25 heavy directions have no symmetry/dynamical home (see (4)).")
    print("      ⇒ building such a Y reproduces both norms but is POSITED, not derived.")

    section("(4) Subalgebra obstruction: no clean '3 light + 25 eaten' from so(8)")
    # maximal subalgebra dims of so(8) (D4), dim 28
    dim_so = lambda n: n * (n - 1) // 2
    subalg = {
        "so(7) [=Spin(7), triality×3]": dim_so(7),
        "so(6)⊕so(2) = u(4)":          dim_so(6) + dim_so(2),
        "so(5)⊕so(3)":                 dim_so(5) + dim_so(3),
        "so(4)⊕so(4) = su(2)^4":       dim_so(4) + dim_so(4),
        "G2":                          14,
        "su(2) (a 3-dim subalgebra)":  3,
    }
    print(f"  dim so(8) = {dim_so(8)}")
    for name, d in subalg.items():
        print(f"    {name:32s} dim = {d}")
    maxproper = max(d for n, d in subalg.items() if d < 28)
    print(f"  max proper subalgebra dim = {maxproper} (= dim Spin(7))  <  25")
    print("  '25 eaten Goldstones' ⇒ unbroken H of dim 3 (e.g. su(2)); reachable only by a")
    print("  CHAIN, not a single step, and the 3 unbroken ≠ the 3 generations (which come")
    print("  from the Z₃ TRIALITY on the generation space, not an su(2) stabilizer of L).")
    assert maxproper == 21 and 25 not in subalg.values()

    section("(5) VERDICT")
    print("""  • The '6/784 gravity↔EW bonus' is NOT independent support — it is the EW bridge
    match v=784a² re-expressed (Σm=9Qa² is definitional).  [deflation]
  • A single Y cannot serve both roles: the so(8)-democratic vacuum gives 28 comparable
    singular values (wrong spectrum), while a rank-3+heavy Y reproduces both norms but
    is hierarchical (not the potential's vacuum) with 25 homeless heavy directions, and
    no single-step so(8)→H supplies them (max proper subalgebra = 21 < 25).
  • So the only consistent reading is TWO OBJECTS on DIFFERENT SPACES:
      – an EW condensate / operator on End(L) = M₂₈ (dim 784, sets v via ‖·‖²_F), and
      – the rank-3 Brannen kernel on the 3-dim GENERATION space (from Z₃ triality).
    They share only the scale a_ℓ and the conjecture 'physical scale = Frobenius² of the
    relevant mass bilinear' (R1).  This DISSOLVES the contradiction but does NOT derive
    v from a_ℓ: v_Higgs stays a conditional SECOND SCALE, tied to a_ℓ only via R1+R2,
    which still have no dynamical home (XiVacuum is 4-dim, not the 28-dim bilinear).
  • Net: G1's '1 input a_ℓ' headline remains CONDITIONAL on R1+R2; the rank tension is
    real for a single Y and is resolved only by the two-object reading.  The genuinely
    open task is a Lagrangian that produces BOTH objects from one scale — not built here.""")


if __name__ == "__main__":
    main()
