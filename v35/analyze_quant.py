#!/usr/bin/env python3
"""Analyze V35 field quantization scan results.

Reads timeseries TSV files from each epsilon run,
extracts key metrics, and outputs a summary table.
"""
import os
import sys
import numpy as np

DATADIR = "data"
EPSILONS = [0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]

def load_ts(eps):
    """Load timeseries for a given epsilon value."""
    if eps == 0:
        tag = "0"
    else:
        tag = str(eps)
    path = os.path.join(DATADIR, f"quant_eps_{tag}_timeseries.tsv")
    if not os.path.exists(path):
        return None
    data = np.loadtxt(path, skiprows=1)
    return data

def analyze():
    print("=" * 90)
    print("V35 Field Quantization Scan Results")
    print("=" * 90)
    print()

    # N=80, L=25 => dx = 50/79 = 0.6329
    dx = 50.0 / 79.0

    # Header
    fmt = "{:>10s} {:>12s} {:>10s} {:>10s} {:>10s} {:>10s} {:>12s} {:>8s}"
    print(fmt.format("eps", "hbar_sim", "E_drift%", "theta_rms", "E_pot(0)",
                      "E_pot(end)", "E_pot_retain", "survive"))
    print("-" * 90)

    results = []
    for eps in EPSILONS:
        data = load_ts(eps)
        if data is None:
            print(f"{eps:>10.4f}  -- NO DATA --")
            continue

        t = data[:, 0]
        E_total = data[:, 9]   # E_total column
        theta_rms_arr = data[:, 10]  # theta_rms column
        E_pot = data[:, 5]     # E_pot column

        E0 = E_total[0]
        E_end = E_total[-1]
        drift_pct = 100.0 * (E_end - E0) / (abs(E0) + 1e-30)

        trms_end = theta_rms_arr[-1]
        trms_max = np.max(theta_rms_arr)

        Ep0 = E_pot[0]
        Ep_end = E_pot[-1]
        Ep_retain = Ep_end / Ep0 if abs(Ep0) > 1e-30 else 0

        # Determine "survival": E_pot should stay negative (bound) and
        # theta_rms should not blow up. Simple heuristic: braid survives
        # if E_pot remains < 0.5 * E_pot(0) (i.e., retains at least 50%
        # of binding energy).
        # Also check if E_total doesn't drift more than 100%.
        survive = "YES" if (abs(drift_pct) < 100 and Ep_retain > 0.5) else "NO"

        hbar_sim = eps * dx if eps > 0 else 0

        results.append({
            'eps': eps, 'hbar_sim': hbar_sim, 'drift': drift_pct,
            'trms': trms_end, 'Ep0': Ep0, 'Ep_end': Ep_end,
            'Ep_retain': Ep_retain, 'survive': survive,
            'trms_max': trms_max
        })

        fmt2 = "{:>10.4f} {:>12.6e} {:>+10.3f} {:>10.4e} {:>10.1f} {:>10.1f} {:>12.3f} {:>8s}"
        print(fmt2.format(eps, hbar_sim, drift_pct, trms_end, Ep0, Ep_end,
                          Ep_retain, survive))

    print()
    print("=" * 90)
    print()

    # Find critical epsilon
    surviving = [r for r in results if r['survive'] == 'YES']
    if surviving:
        crit = max(surviving, key=lambda r: r['eps'])
        print(f"Critical epsilon (largest surviving): {crit['eps']}")
        print(f"  hbar_sim = {crit['hbar_sim']:.6e}")
        print(f"  E_drift = {crit['drift']:+.3f}%")
        print(f"  theta_rms = {crit['trms']:.4e}")
        print(f"  E_pot retention = {crit['Ep_retain']:.3f}")
    else:
        print("WARNING: No runs survived!")

    print()

    # Detailed time evolution for each
    print("TIME EVOLUTION (sampled at t=0, 10, 25, 50):")
    print("-" * 70)
    for eps in EPSILONS:
        data = load_ts(eps)
        if data is None:
            continue
        t = data[:, 0]
        E_total = data[:, 9]
        E_pot = data[:, 5]
        trms = data[:, 10]

        tag = f"eps={eps}"
        for t_sample in [0, 10, 25, 50]:
            idx = np.argmin(np.abs(t - t_sample))
            print(f"  {tag:>15s} t={t[idx]:5.1f}: E={E_total[idx]:+.2e} Ep={E_pot[idx]:+.1f} trms={trms[idx]:.3e}")
        print()

    return results

if __name__ == "__main__":
    results = analyze()
