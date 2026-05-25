#!/usr/bin/env python3
"""
study.py — reproducible investigation of the lepton phase law  φ = Q/3.

The Brannen circulant mass matrix is   M = a (I + ξ S + ξ̄ S²),  ξ = t e^{iφ},
with S the 3-cycle (Z₃ generation shift).  Its eigenvalues are √m_k.  Koide's
Q = (Σ√m)² ... = (1+2t²)/3 fixes the amplitude (t²=1/2 ⇔ Q=2/3); the phase φ is
a SEPARATE degree of freedom (the family CP phase).  This script studies the
empirical relation  φ = Q/3  (lepton: φ = 2/9):  its forms, its sector-dependence,
its precision, and physical leads.

Run:  python3 study.py
"""
import numpy as np

PDG = {  # GeV, central PDG-2024
    'lepton'    : [0.51099895e-3, 105.6583755e-3, 1776.86e-3],
    'up-quark'  : [2.16e-3, 1.27, 172.69],
    'down-quark': [4.67e-3, 93.4e-3, 4.18],
}
DM_TAU = 0.12e-3  # GeV, dominant lepton uncertainty

def brannen_fit(masses):
    """Return Q, t², φ (principal Z₃ branch), cos3φ for a mass triple."""
    sm = np.sqrt(np.array(masses, float)); a = sm.sum()/3
    Q = (sm**2).sum()/sm.sum()**2
    t2 = (3*Q - 1)/2
    k = np.arange(3); ang = 2*np.pi*k/3
    C = np.sum((sm/a - 1)*np.cos(ang)); S = np.sum((sm/a - 1)*np.sin(ang))
    phi = (-np.arctan2(S, C)) % (2*np.pi/3)     # reduce to principal Z₃ branch
    return Q, t2, phi, np.cos(3*phi)

def main():
    print("="*72)
    print("1. The relation across sectors  (φ vs Q/3)")
    print("="*72)
    print(f"{'sector':12s} {'Q':>8s} {'t^2':>8s} {'phi':>9s} {'Q/3':>9s} "
          f"{'3phi':>8s} {'|phi-Q/3|':>10s}")
    for s, m in PDG.items():
        Q, t2, phi, _ = brannen_fit(m)
        print(f"{s:12s} {Q:8.5f} {t2:8.5f} {phi:9.5f} {Q/3:9.5f} {3*phi:8.5f} "
              f"{abs(phi-Q/3):10.5f}")
    print("  -> φ = Q/3 holds for LEPTONS only; quarks miss by O(0.1).  Lepton-specific.")

    print("\n" + "="*72)
    print("2. The lepton relation, three equivalent forms (all = 2 at the point)")
    print("="*72)
    Q, t2, phi, c3 = brannen_fit(PDG['lepton'])
    print(f"  3φ          = {3*phi:.6f}   (= Q  : phase = Koide ratio)")
    print(f"  3Q          = {3*Q:.6f}")
    print(f"  9φ          = {9*phi:.6f}")
    print(f"  1 + 2 t²    = {1+2*t2:.6f}   <- 9φ = 1+2|ξ|² : phase locked to modulus")
    print("  So the content is:  9·arg(ξ) = 1 + 2|ξ|²  (= 2 at the Koide point |ξ|²=1/2).")

    print("\n" + "="*72)
    print("3. Precision: φ=2/9 is Koide-tight (both m_τ-limited)")
    print("="*72)
    me, mmu, mtau = PDG['lepton']
    rows = [brannen_fit([me, mmu, mtau + s*DM_TAU]) for s in np.linspace(-1, 1, 5)]
    Qs, _, phis, _ = map(np.array, zip(*rows))
    dQ, dphi = (Qs.max()-Qs.min())/2, (phis.max()-phis.min())/2
    print(f"  Q  = {Q:.8f} ± {dQ:.1e}   |Q-2/3| = {abs(Q-2/3):.1e} = {abs(Q-2/3)/dQ:.2f}σ")
    print(f"  φ  = {phi:.8f} ± {dphi:.1e}   |φ-2/9| = {abs(phi-2/9):.1e} = {abs(phi-2/9)/dphi:.2f}σ")
    print("  -> φ=2/9 and Q=2/3 agree with data at the SAME ~0.9σ / 1e-5 level.")

    print("\n" + "="*72)
    print("4. Where φ puts the masses (physical: the lightest state)")
    print("="*72)
    for k in range(3):
        v = 1 + np.sqrt(2)*np.cos(phi + 2*np.pi*k/3)
        print(f"  k={k}: √m_k/a = 1+√2 cos(φ+2πk/3) = {v:+.4f}  -> m_k/a² = {v**2:.4f}")
    print(f"  electron (k=1) is small: cos(φ+2π/3)={np.cos(phi+2*np.pi/3):.4f} ≈ -1/√2={-1/np.sqrt(2):.4f}")
    print(f"  m_e -> 0 boundary at φ = π/12 = {np.pi/12:.5f}; actual φ = {phi:.5f}; gap {np.pi/12-phi:.5f}")
    print("  -> φ lies just inside the m_e→0 edge: the phase keeps the electron light but ≠0.")

    print("\n" + "="*72)
    print("5. 'Is it alone?'  The Q/3 = 2/9 chain and the ~0.22 cluster")
    print("="*72)
    print(f"  Koide:     Q = dimG₂/dimSpin7 = 14/21 = 2/3 = {2/3:.5f}")
    print(f"  phase:     φ = Q/3 = 2/9             = {2/9:.5f}   (lepton, 1e-5)")
    print(f"  gauge:     sin²θ_W(tree, Pati-Salam) = 2/9 = {2/9:.5f}   (structural; M_Z value 0.231)")
    print(f"  Cabibbo λ  ≈ 0.2250   (1% from 2/9, NOT 1e-5)")
    print( "  -> at 1e-5 ONLY the lepton phase hits 2/9; the gauge 2/9 is a tree/structural")
    print( "     claim (same Q/3), the CKM ~0.22 numbers are %-level near-misses, not siblings.")

if __name__ == '__main__':
    main()
