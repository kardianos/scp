#!/usr/bin/env python3
"""
V35 Parameter Space Mapping: Find the region where bound states match hydrogen.

Sweeps (hbar, m_eff, V_depth, n_power) and scores each point against
physical hydrogen observables. The 1D radial solver runs in ~0.1s per
evaluation, so 10,000 points takes ~15 minutes.

Targets:
  - Bohr ratio: <r>_0 / r_braid ~ 53,000 (within 0.5-2× = 26,000-106,000)
  - Energy ratio: E_1/E_0 ~ 0.25 (hydrogen 1/n² scaling)
  - Bound states: 2-10 for l=0
  - l=1 has at least 1 bound state

Output: data/param_sweep.tsv with all evaluations + scores
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
import sys, os, time

R_BRAID = 5.0       # braid radius in code units
R_CORE = 2.0        # repulsive core radius

# Target observables (hydrogen)
TARGET_BOHR_RATIO = 53000.0
TARGET_E_RATIO = 0.25    # E_1/E_0 for hydrogen

def solve_radial(hbar, m_eff, V_depth, n_power, l=0, r_min=0.5, N_r=5000):
    """Solve radial Schrodinger equation. Returns (energies, mean_r) for bound states."""
    # Adaptive r_max based on expected Bohr radius
    alpha = V_depth * 1.27  # effective coupling
    if hbar > 0 and m_eff > 0:
        a0_est = hbar**2 / (m_eff * alpha + 1e-30)
        r_max = max(100, min(20 * a0_est, 5e7))
    else:
        r_max = 1000

    r = np.linspace(r_min, r_max, N_r)
    dr = r[1] - r[0]

    # Potential: V(r) = -V_depth * 1.27 / (r/R_BRAID)^n_power  for r > R_CORE
    #            V(r) = repulsive for r < R_CORE
    V = np.zeros(N_r)
    for i in range(N_r):
        if r[i] < R_CORE:
            V[i] = V_depth * 10.0 * (R_CORE/r[i] - 1)  # repulsive core
        else:
            V[i] = -V_depth * 1.27 / (r[i]/R_BRAID)**n_power

    # Add centrifugal barrier
    if l > 0:
        V += hbar**2 * l*(l+1) / (2*m_eff * r**2 + 1e-30)

    # Hamiltonian: -hbar^2/(2m) d^2/dr^2 + V(r)
    kinetic = hbar**2 / (2 * m_eff * dr**2)
    diag = 2*kinetic + V
    offdiag = -kinetic * np.ones(N_r - 1)

    try:
        # Find lowest eigenvalues only
        n_eig = min(20, N_r - 2)
        energies, vectors = eigh_tridiagonal(diag, offdiag,
                                              select='i', select_range=(0, n_eig-1))
    except Exception:
        return np.array([]), np.array([])

    # Find bound states (E < 0)
    bound_mask = energies < 0
    bound_E = energies[bound_mask]
    bound_vecs = vectors[:, bound_mask]

    # Compute <r> for each bound state
    mean_r = np.zeros(len(bound_E))
    for j in range(len(bound_E)):
        psi2 = bound_vecs[:, j]**2
        norm = np.sum(psi2) * dr
        if norm > 1e-30:
            mean_r[j] = np.sum(r * psi2) * dr / norm

    return bound_E, mean_r


def score_point(hbar, m_eff, V_depth, n_power):
    """Score a parameter point against hydrogen observables. Lower = better."""
    # Solve l=0
    E0, r0 = solve_radial(hbar, m_eff, V_depth, n_power, l=0)
    n_bound = len(E0)

    if n_bound < 2:
        return 100.0, {'n_bound': n_bound, 'bohr_ratio': 0, 'e_ratio': 0, 'has_l1': False}

    bohr_ratio = r0[0] / R_BRAID
    e_ratio = E0[1] / E0[0] if abs(E0[0]) > 1e-30 else 0

    # Solve l=1
    E1, r1 = solve_radial(hbar, m_eff, V_depth, n_power, l=1)
    has_l1 = len(E1) > 0

    # Score components
    s_bohr = abs(np.log10(bohr_ratio / TARGET_BOHR_RATIO + 1e-30))  # log distance
    s_eratio = abs(e_ratio - TARGET_E_RATIO) / 0.25
    s_nbound = 0 if 2 <= n_bound <= 10 else abs(n_bound - 5) / 5.0
    s_l1 = 0 if has_l1 else 2.0

    total = s_bohr + s_eratio + s_nbound + s_l1

    info = {
        'n_bound': n_bound,
        'bohr_ratio': bohr_ratio,
        'e_ratio': e_ratio,
        'has_l1': has_l1,
        'E0': E0[0] if n_bound > 0 else 0,
        'r0': r0[0] if n_bound > 0 else 0,
        's_bohr': s_bohr,
        's_eratio': s_eratio,
    }
    return total, info


def main():
    os.makedirs('data', exist_ok=True)

    # Parameter space grid
    hbar_vals = np.logspace(-1, 4, 30)      # 0.1 to 10,000
    meff_vals = np.logspace(-4, 1, 25)       # 0.0001 to 10
    vdepth_vals = [0.5, 1.0, 2.0, 5.0]      # potential depth multiplier
    npower_vals = [1.0, 1.1, 1.189, 1.3]     # power law exponent

    total_evals = len(hbar_vals) * len(meff_vals) * len(vdepth_vals) * len(npower_vals)
    print(f"Parameter sweep: {total_evals} evaluations")
    print(f"  hbar:    {len(hbar_vals)} values [{hbar_vals[0]:.2f} - {hbar_vals[-1]:.0f}]")
    print(f"  m_eff:   {len(meff_vals)} values [{meff_vals[0]:.4f} - {meff_vals[-1]:.1f}]")
    print(f"  V_depth: {len(vdepth_vals)} values {vdepth_vals}")
    print(f"  n_power: {len(npower_vals)} values {npower_vals}")
    print()

    results = []
    best_score = 999
    best_params = None
    t0 = time.time()
    count = 0

    with open('data/param_sweep.tsv', 'w') as f:
        f.write("hbar\tm_eff\tV_depth\tn_power\tscore\tn_bound\tbohr_ratio\te_ratio\thas_l1\tE0\tr0\n")

        for V_depth in vdepth_vals:
            for n_power in npower_vals:
                for m_eff in meff_vals:
                    for hbar in hbar_vals:
                        score, info = score_point(hbar, m_eff, V_depth, n_power)
                        count += 1

                        f.write(f"{hbar:.6e}\t{m_eff:.6e}\t{V_depth:.2f}\t{n_power:.3f}\t"
                                f"{score:.4f}\t{info['n_bound']}\t{info['bohr_ratio']:.1f}\t"
                                f"{info['e_ratio']:.4f}\t{int(info['has_l1'])}\t"
                                f"{info.get('E0',0):.6e}\t{info.get('r0',0):.1f}\n")

                        if score < best_score:
                            best_score = score
                            best_params = (hbar, m_eff, V_depth, n_power)
                            best_info = info

                        if count % 1000 == 0:
                            elapsed = time.time() - t0
                            rate = count / elapsed
                            eta = (total_evals - count) / rate
                            print(f"  {count}/{total_evals} ({100*count/total_evals:.0f}%) "
                                  f"best={best_score:.3f} [{elapsed:.0f}s, eta {eta:.0f}s]")

        f.flush()

    elapsed = time.time() - t0
    print(f"\nComplete: {count} evaluations in {elapsed:.0f}s ({count/elapsed:.0f}/s)")

    # Report best results
    print(f"\n{'='*70}")
    print(f"BEST MATCH TO HYDROGEN")
    print(f"{'='*70}")
    print(f"Score: {best_score:.4f}")
    print(f"Parameters:")
    print(f"  hbar    = {best_params[0]:.4f}")
    print(f"  m_eff   = {best_params[1]:.6f}")
    print(f"  V_depth = {best_params[2]:.2f}")
    print(f"  n_power = {best_params[3]:.3f}")
    print(f"Observables:")
    print(f"  Bohr ratio:  {best_info['bohr_ratio']:.0f}  (target: {TARGET_BOHR_RATIO:.0f})")
    print(f"  E1/E0:       {best_info['e_ratio']:.4f}  (target: {TARGET_E_RATIO:.4f})")
    print(f"  N_bound:     {best_info['n_bound']}")
    print(f"  l=1 exists:  {best_info['has_l1']}")

    # Find top 20 matches
    print(f"\n--- Top 20 Parameter Combinations ---")
    print(f"{'Score':>8s}  {'hbar':>10s}  {'m_eff':>10s}  {'V_depth':>7s}  {'n':>5s}  "
          f"{'Bohr':>8s}  {'E1/E0':>6s}  {'N_b':>3s}  {'l1':>2s}")

    # Re-read and sort
    import csv
    rows = []
    with open('data/param_sweep.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    rows.sort(key=lambda r: float(r['score']))

    for row in rows[:20]:
        print(f"{float(row['score']):8.3f}  {float(row['hbar']):10.4f}  "
              f"{float(row['m_eff']):10.6f}  {float(row['V_depth']):7.2f}  "
              f"{float(row['n_power']):5.3f}  {float(row['bohr_ratio']):8.0f}  "
              f"{float(row['e_ratio']):6.3f}  {int(row['n_bound']):3d}  "
              f"{row['has_l1']:>2s}")

    # Identify the constraint surface: what combinations give Bohr ratio ~ 53,000?
    print(f"\n--- Constraint: Bohr ratio in [26000, 106000] ---")
    print(f"{'hbar':>10s}  {'m_eff':>10s}  {'V_depth':>7s}  {'n':>5s}  "
          f"{'Bohr':>8s}  {'E1/E0':>6s}  {'N_b':>3s}")
    viable = [r for r in rows
              if 26000 < float(r['bohr_ratio']) < 106000 and int(r['n_bound']) >= 2]
    for row in viable[:30]:
        print(f"{float(row['hbar']):10.4f}  {float(row['m_eff']):10.6f}  "
              f"{float(row['V_depth']):7.2f}  {float(row['n_power']):5.3f}  "
              f"{float(row['bohr_ratio']):8.0f}  {float(row['e_ratio']):6.3f}  "
              f"{int(row['n_bound']):3d}")

    # Key relationship: hbar^2 / m_eff determines the Bohr radius
    print(f"\n--- Key Relationship: hbar²/m_eff ---")
    print(f"For Bohr ratio ~ 53,000: hbar²/m_eff ≈ ?")
    for row in viable[:10]:
        h = float(row['hbar'])
        m = float(row['m_eff'])
        print(f"  hbar={h:.3f}, m_eff={m:.6f} → hbar²/m_eff = {h*h/m:.1f}")


if __name__ == '__main__':
    main()
