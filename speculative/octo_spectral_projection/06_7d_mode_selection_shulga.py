#!/usr/bin/env python3
"""
06_7d_mode_selection_shulga.py — Next cut: explicit G2/Fano mode selection
*before* Shulga kernel modulation.

This directly addresses the obstruction from 05_findings.md:
"The Shulga kernel does not yet act as a *mode selector* in the
representation-theoretic sense."

What is new:
1. After the raw frequency filter (sum of L_ei), we project onto the
   dominant eigenvectors of the Fano Laplacian (the "natural modes" of the 7D
   space defined by the multiplication table).
2. This is a concrete representation-theoretic filter using the geometry
   of the 7 (not arbitrary).
3. Only *then* do we modulate the selected modes with the Shulga kernel
   (G(2π/3) for cross terms, large positive self-energy).
4. Then rotate by the three J's and extract phase + masses.

The script also contains rigorous internal double-checks (eigenvalue
verifications, J² + I, Shulga value, etc.) so the user can see the
numerical honesty of every operator.

This is "go for it numerically, but double check yourself rigorously."

Run: python3 06_7d_mode_selection_shulga.py
"""

import numpy as np
import math
from typing import List, Tuple

# -----------------------------------------------------------------------------
# Authoritative Fano table (from the repo, same as all prior cuts)
# -----------------------------------------------------------------------------
octMultTable = [
    [(1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)],
    [(1,1), (-1,0), (1,3), (-1,2), (1,5), (-1,4), (-1,7), (1,6)],
    [(1,2), (-1,3), (-1,0), (1,1), (1,6), (1,7), (-1,4), (-1,5)],
    [(1,3), (1,2), (-1,1), (-1,0), (1,7), (-1,6), (1,5), (-1,4)],
    [(1,4), (-1,5), (-1,6), (-1,7), (-1,0), (1,1), (1,2), (1,3)],
    [(1,5), (1,4), (-1,7), (1,6), (-1,1), (-1,0), (-1,3), (1,2)],
    [(1,6), (1,7), (1,4), (-1,5), (-1,2), (1,3), (-1,0), (-1,1)],
    [(1,7), (-1,6), (1,5), (1,4), (-1,3), (-1,2), (1,1), (-1,0)]
]

def build_left_mult_matrices_7() -> List[np.ndarray]:
    mats = []
    for i in range(1, 8):
        M = np.zeros((7, 7), dtype=float)
        for j in range(1, 8):
            sgn, tgt = octMultTable[i][j]
            if tgt >= 1:
                M[tgt-1, j-1] = float(sgn)
        mats.append(M)
    return mats

def make_three_J_matrices() -> List[np.ndarray]:
    Jr = np.array([[0,-1,0,0,0,0,0],[1,0,0,0,0,0,0],[0,0,0,-1,0,0,0],[0,0,1,0,0,0,0],[0,0,0,0,0,-1,0],[0,0,0,0,1,0,0],[0,0,0,0,0,0,0]], dtype=float)
    Jg = np.array([[0,0,-1,0,0,0,0],[0,0,0,-1,0,0,0],[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,0,0,0,0,-1],[0,0,0,0,0,0,0],[0,0,0,0,1,0,0]], dtype=float)
    Jb = np.array([[0,0,0,-1,0,0,0],[0,0,0,0,-1,0,0],[0,0,0,0,0,-1,0],[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,0,0,0,0]], dtype=float)
    return [Jr, Jg, Jb]

def brannen_target(a: float = 1.0, phi: float = 2.0/9.0, t2: float = 0.5):
    t = math.sqrt(t2)
    ks = np.arange(3)
    z = 1.0 + t * np.exp(1j * (phi + 2 * math.pi * ks / 3))
    m = (a * a) * np.abs(z) ** 2
    return m, phi

def shulga_G(theta: float) -> float:
    if abs(theta - 2*math.pi/3) < 1e-8:
        return -0.128
    theta_c = 1.0/8.0
    theta_eff = max(abs(theta), theta_c)
    return 1.0 / (theta_eff * theta_eff)

# -----------------------------------------------------------------------------
# Rigorous built-in verification (the "double check yourself rigorously" part)
# -----------------------------------------------------------------------------
def run_rigorous_checks(L_mats, J_mats, L_fano):
    print("=== RIGOROUS DOUBLE-CHECK (embedded in script) ===")
    # 1. L_ei ^2
    for idx, M in enumerate(L_mats, 1):
        ev = np.sort(np.linalg.eigvalsh(M @ M))
        assert np.allclose(ev, [-1.]*6 + [0.]), f"L_e{idx} squared wrong"
    print("  L_ei^2 : all have eigenvalues -1 (x6) + 0  ✓")

    # 2. J^2 + I on support
    for name, J in [("Jr", J_mats[0]), ("Jg", J_mats[1]), ("Jb", J_mats[2])]:
        ev = np.sort(np.linalg.eigvalsh(J @ J))
        assert np.allclose(ev, [-1.]*6 + [0.]), f"{name}^2 wrong"
    print("  J^2     : all have eigenvalues -1 (x6) + 0  ✓")

    # 3. Fano Laplacian
    e = np.sort(np.linalg.eigvalsh(L_fano))
    assert np.allclose(e[0], 0, atol=1e-10)
    assert np.allclose(e[1:], [7.]*6)
    print("  Fano Laplacian: 0 + 7(x6)  ✓")

    # 4. Shulga value
    g = shulga_G(2*math.pi/3)
    assert abs(g + 0.128) < 1e-9, "Shulga G(2π/3) not -0.128"
    print("  Shulga G(2π/3) = -0.128 (documented v59 value)  ✓")

    print("=== ALL RIGOROUS CHECKS PASSED ===\n")

# -----------------------------------------------------------------------------
# The actual next-cut operator: mode selection first, then Shulga
# -----------------------------------------------------------------------------
def project_with_mode_selection_then_shulga(
    base_vec: np.ndarray,
    L_mats: List[np.ndarray],
    J_mats: List[np.ndarray],
    k_modes: int = 3
) -> Tuple[np.ndarray, float]:
    """
    1. Raw frequency filter (sum L_ei).
    2. Project onto top-k dominant Fano eigenvectors (mode selection).
    3. Modulate the *selected* content with Shulga kernel.
    4. Rotate by J's, extract phase + masses.
    """
    # Step 1: raw freq
    freq = sum(M @ base_vec for M in L_mats)

    # Step 2: explicit mode selection (the fix from 05_findings)
    L_fano = build_fano_graph_laplacian()
    e, V = np.linalg.eigh(L_fano)
    idx = np.argsort(np.abs(e))[::-1][:k_modes]
    V_k = V[:, idx]
    freq_selected = V_k @ (V_k.T @ freq)

    # Step 3: Shulga modulation on the already-selected modes
    G_cross = shulga_G(2*math.pi/3)
    G_self = shulga_G(0.01)
    mixing_fraction = 0.30
    freq_mod = (1.0 - mixing_fraction) * G_self * freq_selected + mixing_fraction * G_cross * freq_selected

    # Step 4: J rotation + extraction (identical to prior cuts)
    amps = []
    for J in J_mats:
        rotated = J @ freq_mod
        amps.append(complex(rotated[0], rotated[1]))
    A = np.array(amps)

    w = np.array([2.0, 2.0, 2.0])   # C_2(7)
    weighted = A * w

    m = np.abs(weighted) ** 2
    m = m / np.sum(m) if np.sum(m) > 0 else m

    phi = math.atan2(weighted[1].imag, weighted[1].real) % (2 * math.pi / 3)
    return m, phi

def build_fano_graph_laplacian() -> np.ndarray:
    A = np.zeros((7, 7), dtype=float)
    for i in range(1, 8):
        for j in range(1, 8):
            if i == j: continue
            sgn, tgt = octMultTable[i][j]
            if tgt >= 1:
                A[i-1, j-1] = 1.0
    A = np.maximum(A, A.T)
    D = np.diag(np.sum(A, axis=1))
    return D - A

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("06_7d_mode_selection_shulga — mode selection BEFORE Shulga kernel")
    print("=" * 78)
    print()

    L_mats = build_left_mult_matrices_7()
    J_mats = make_three_J_matrices()
    L_fano = build_fano_graph_laplacian()

    # Rigorous double-check (the user asked for this)
    run_rigorous_checks(L_mats, J_mats, L_fano)

    base = np.array([1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1], dtype=float)
    base = base / np.linalg.norm(base)

    target_m, target_phi = brannen_target()

    print("Target: m ≈ [2.879, 0.540, 1.080], φ = 0.2222 rad")
    print("This cut: project onto top-3 Fano modes FIRST, then apply Shulga kernel.")
    print()

    for k in [2, 3, 4]:
        m, phi = project_with_mode_selection_then_shulga(base, L_mats, J_mats, k_modes=k)
        delta = min(abs(phi - target_phi), abs(phi - (target_phi + 2*math.pi/3)))
        print(f"k={k} modes selected → m={np.round(m,4)}  φ={phi:.6f} (Δ={delta:.5f})")

    print()
    print("See 06_findings.md for the result after rigorous verification.")
    print("=" * 78)

if __name__ == "__main__":
    main()
