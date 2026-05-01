#!/usr/bin/env python3
"""Windowed FFT and raw trace inspection for v54 points trajectory.

Diagnostic follow-up to freq_analysis.py, which revealed the default points.c
configuration is over-damped (KE decays 60x over T=400). Here we:
  1. Plot raw |v|(t) for a sample of particles to see actual motion
  2. Run separate FFTs over early / mid / late thirds
  3. Report per-window frequency distribution and KE
"""
import struct, sys
import numpy as np
from pathlib import Path
import csv
from freq_analysis import load_traj  # reuse the loader


def windowed_fft(vmag, dt_sample, frac_start, frac_end):
    """FFT over a sub-window of the time series (fraction of total length)."""
    n = vmag.shape[0]
    i0 = int(frac_start * n)
    i1 = int(frac_end * n)
    v = vmag[i0:i1]
    v = v - v.mean(axis=0, keepdims=True)
    L = v.shape[0]
    w = np.hanning(L)[:, None]
    F = np.fft.rfft(v * w, axis=0)
    freqs = np.fft.rfftfreq(L, d=dt_sample)
    P = np.abs(F)**2
    # Dominant freq (skip DC)
    if P.shape[0] > 1:
        idx = np.argmax(P[1:], axis=0) + 1
        f_dom = freqs[idx]
    else:
        f_dom = np.zeros(v.shape[1])
    ens = P.mean(axis=1)
    ke_mean = (vmag[i0:i1]**2).mean()  # proxy for kinetic energy
    return freqs, P, ens, f_dom, ke_mean, L


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'traj.bin'
    traj = load_traj(path)
    P = traj['samples']['p']
    dt_s = traj['dt'] * traj['sps']
    n_samples, n_tracked = P.shape

    # velocity magnitude
    vx = P['vx'].astype(np.float64); vy = P['vy'].astype(np.float64); vz = P['vz'].astype(np.float64)
    vmag = np.sqrt(vx*vx + vy*vy + vz*vz)
    t = np.arange(n_samples) * dt_s

    print(f"n_samples={n_samples}  n_tracked={n_tracked}  dt_sample={dt_s}  T={t[-1]:.1f}")
    print(f"\nGlobal |v| statistics over time (bin averages):")
    print(f"  {'t_window':>16} {'<|v|>':>10} {'<|v|²>':>10} {'max|v|':>10}")
    windows = [(0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.5), (0.5, 1.0)]
    for fa, fb in windows:
        i0, i1 = int(fa*n_samples), int(fb*n_samples)
        seg = vmag[i0:i1]
        print(f"  [{t[i0]:5.1f},{t[i1-1]:5.1f}]  {seg.mean():10.6f} {(seg**2).mean():10.6e} {seg.max():10.4f}")

    # --- Dump raw |v|(t) for 10 particles ---
    print(f"\nDumping v54/vmag_traces.tsv (10 sample particles, raw |v|(t))")
    sample_ids = list(range(0, n_tracked, max(1, n_tracked//10)))[:10]
    out = Path(path).parent
    with open(out/'vmag_traces.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['t'] + [f'vmag_p{i}' for i in sample_ids])
        # Decimate for readable file (every 20th sample = dt=0.4)
        for i in range(0, n_samples, 20):
            row = [f"{t[i]:.3f}"] + [f"{vmag[i, k]:.6f}" for k in sample_ids]
            w.writerow(row)

    # --- FFT per time window ---
    print("\nFFT per time window (dominant frequency population):")
    frac_windows = {
        'early (0-25%)': (0.00, 0.25),
        'mid   (25-50%)': (0.25, 0.50),
        'late  (50-100%)': (0.50, 1.00),
    }
    psd_rows = {}
    for name, (a, b) in frac_windows.items():
        freqs, Pspec, ens, f_dom, ke, L = windowed_fft(vmag, dt_s, a, b)
        print(f"\n  {name}  (L={L} samples, T={L*dt_s:.1f}, df={1/(L*dt_s):.5f})")
        print(f"    <|v|²> over window  : {ke:.6e}")
        print(f"    median f_dom        : {np.median(f_dom):.5f}  omega={2*np.pi*np.median(f_dom):.4f}")
        print(f"    mean   f_dom        : {np.mean(f_dom):.5f}  std={np.std(f_dom):.5f}")
        print(f"    range              : [{f_dom.min():.5f}, {f_dom.max():.5f}]")
        # Histogram
        n_bins = 20
        hc, he = np.histogram(f_dom, bins=n_bins, range=(0, max(f_dom.max()*1.2, 0.5)))
        mx = hc.max()
        for i in range(n_bins):
            if hc[i] > 0:
                bar = '#' * int(50 * hc[i] / max(mx, 1))
                print(f"    f=[{he[i]:5.3f},{he[i+1]:5.3f}) n={hc[i]:3d} {bar}")
        psd_rows[name] = (freqs, ens)

    # Save stratified PSD
    # Use the 'early' window for best-signal ensemble PSD
    freqs_early = psd_rows['early (0-25%)'][0]
    rows = []
    for i in range(len(freqs_early)):
        r = [f"{freqs_early[i]:.5f}", f"{2*np.pi*freqs_early[i]:.5f}"]
        for name in psd_rows:
            _, ens = psd_rows[name]
            # Window lengths may differ; pad with NaN if out of range
            if i < len(ens):
                r.append(f"{ens[i]:.5g}")
            else:
                r.append("nan")
        rows.append(r)
    with open(out/'psd_by_window.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['freq','omega'] + list(psd_rows.keys()))
        w.writerows(rows)
    print(f"\nWrote {out}/vmag_traces.tsv and {out}/psd_by_window.tsv")

    # Plot
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        ax = axes[0,0]
        for k in sample_ids[:5]:
            ax.plot(t, vmag[:, k], alpha=0.7, lw=0.6, label=f'pid{int(traj["ids"][k])}')
        ax.set_yscale('log')
        ax.set_xlabel('t'); ax.set_ylabel('|v|')
        ax.set_title('Raw |v|(t) for 5 tracked particles (log y)')
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[0,1]
        v2_t = (vmag**2).mean(axis=1)
        ax.semilogy(t, v2_t)
        ax.set_xlabel('t'); ax.set_ylabel('<|v|²>')
        ax.set_title(f'Ensemble <|v|²>(t) — decay from 0.025 to {v2_t[-1]:.1e}')
        ax.grid(alpha=0.3)

        ax = axes[1,0]
        colors = ['#1f77b4','#ff7f0e','#2ca02c']
        for (name, (freqs, ens)), c in zip(psd_rows.items(), colors):
            fmask = (freqs > 0) & (freqs < 5.0)
            ax.semilogy(freqs[fmask], ens[fmask], label=name, color=c, lw=0.8)
        ax.set_xlabel('frequency (1/t)'); ax.set_ylabel('PSD')
        ax.set_title('Ensemble PSD per time window')
        ax.legend(); ax.grid(alpha=0.3)

        ax = axes[1,1]
        # Dominant freq histogram for early window
        freqs, _, _, f_dom_early, _, _ = windowed_fft(vmag, dt_s, 0.0, 0.25)
        ax.hist(f_dom_early, bins=30, range=(0, f_dom_early.max()*1.2), color='#1f77b4', edgecolor='k')
        ax.set_xlabel('dominant frequency'); ax.set_ylabel('# particles')
        ax.set_title(f'Dominant freq: early window (T={0.25*400:.0f})')
        ax.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(out/'freq_windowed.png', dpi=120)
        print(f"Wrote {out}/freq_windowed.png")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == '__main__':
    main()
