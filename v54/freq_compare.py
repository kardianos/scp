#!/usr/bin/env python3
"""Compare FFT spectra between dissipative and conservative runs.

Reads traj.bin (dissipative baseline) and traj_cons.bin (conservative),
produces side-by-side spectra and tests the four ball-and-wire predictions:
  1. Clean peak in conservative, not in dissipative
  2. Peak near predicted natural frequency
  3. freq correlates with sqrt(neighbor_count)
  4. Spatially-close particles share frequencies
"""
import sys
import numpy as np
from pathlib import Path
import csv
from freq_analysis import load_traj


def per_particle_spectrum(traj, window_frac=(0.25, 1.0)):
    """FFT per particle over a time window."""
    P = traj['samples']['p']
    dt_s = traj['dt'] * traj['sps']
    n_samples, n_tracked = P.shape
    vmag = np.sqrt(P['vx']**2 + P['vy']**2 + P['vz']**2).astype(np.float64)

    i0 = int(window_frac[0] * n_samples)
    i1 = int(window_frac[1] * n_samples)
    seg = vmag[i0:i1]
    seg = seg - seg.mean(axis=0, keepdims=True)
    # Also detrend linearly (remove slow drift)
    L = seg.shape[0]
    x = np.arange(L)
    for k in range(n_tracked):
        slope, inter = np.polyfit(x, seg[:, k], 1)
        seg[:, k] -= slope * x + inter
    w = np.hanning(L)[:, None]
    F = np.fft.rfft(seg * w, axis=0)
    freqs = np.fft.rfftfreq(L, d=dt_s)
    Pspec = (np.abs(F) ** 2)

    dom_idx = np.argmax(Pspec[1:], axis=0) + 1
    freq_dom = freqs[dom_idx]

    mean_nn = P['nn'].astype(float).mean(axis=0)
    mean_x = P['x'].mean(axis=0); mean_y = P['y'].mean(axis=0); mean_z = P['z'].mean(axis=0)

    return {
        'freqs': freqs, 'Pspec': Pspec, 'ens_psd': Pspec.mean(axis=1),
        'freq_dom': freq_dom, 'mean_nn': mean_nn,
        'mean_x': mean_x, 'mean_y': mean_y, 'mean_z': mean_z,
        'n_samples_window': L, 'dt_s': dt_s,
        'ke_window': (vmag[i0:i1] ** 2).mean(),
    }


def ascii_spectrum(freqs, psd, f_max=2.0, n_bars=40, height=14):
    """Print an ASCII chart of the spectrum."""
    mask = (freqs > 0) & (freqs < f_max)
    f = freqs[mask]; p = psd[mask]
    # Bin into n_bars
    bin_edges = np.linspace(0, f_max, n_bars + 1)
    binned = np.zeros(n_bars)
    for i in range(n_bars):
        m = (f >= bin_edges[i]) & (f < bin_edges[i+1])
        if m.any():
            binned[i] = p[m].max()
    log_p = np.log10(binned + 1e-20)
    # Scale to [0, height]
    lo, hi = log_p.min(), log_p.max()
    if hi == lo:
        hi = lo + 1
    norm = (log_p - lo) / (hi - lo)
    rows = []
    for r in range(height, 0, -1):
        line = ''
        for i in range(n_bars):
            y = norm[i] * height
            line += '#' if y >= r else ' '
        rows.append(f'  {r:2d} | {line}')
    rows.append('     +' + '-' * n_bars)
    rows.append(f'       0{" " * (n_bars//2 - 3)}f={f_max/2:.2f}{" " * (n_bars//2 - 4)}{f_max:.2f}')
    print(f'  log10 power, max={hi:.2f}, min={lo:.2f}')
    print('\n'.join(rows))


def report(name, spec):
    print(f"\n===== {name} =====")
    print(f"  window samples : {spec['n_samples_window']}  T_window = {spec['n_samples_window']*spec['dt_s']:.1f}")
    print(f"  <|v|²> window  : {spec['ke_window']:.3e}")
    f_dom = spec['freq_dom']
    print(f"  dominant frequency population:")
    print(f"    median : {np.median(f_dom):.4f}   omega={2*np.pi*np.median(f_dom):.3f}")
    print(f"    mean   : {np.mean(f_dom):.4f}   std={np.std(f_dom):.4f}")
    print(f"    range  : [{np.min(f_dom):.4f}, {np.max(f_dom):.4f}]")
    # Histogram
    nn_bins = 30
    mx = f_dom.max()
    hc, he = np.histogram(f_dom, bins=nn_bins, range=(0, max(mx*1.2, 0.5)))
    mxc = hc.max()
    print(f"  histogram (bin width = {he[1]-he[0]:.4f}):")
    for i in range(nn_bins):
        if hc[i] > 0:
            bar = '#' * int(50 * hc[i] / max(mxc, 1))
            print(f"    f=[{he[i]:5.3f},{he[i+1]:5.3f}) n={hc[i]:4d} {bar}")

    # Correlations
    nn = spec['mean_nn']
    r_nn = np.corrcoef(nn, f_dom)[0, 1]
    r_sqrt = np.corrcoef(np.sqrt(nn + 1e-9), f_dom)[0, 1]
    print(f"  Pearson r(nn,       f_dom) = {r_nn:+.4f}")
    print(f"  Pearson r(sqrt(nn), f_dom) = {r_sqrt:+.4f}")

    # Stratify
    print(f"  stratified by mean neighbor count:")
    for lo, hi in [(0,6),(6,12),(12,18),(18,24),(24,32)]:
        m = (nn >= lo) & (nn < hi)
        if m.any():
            print(f"    nn=[{lo},{hi}):  n={int(m.sum()):3d}  median_f={np.median(f_dom[m]):.4f}  std_f={np.std(f_dom[m]):.4f}")

    # Spatial coherence
    mx_, my_, mz_ = spec['mean_x'], spec['mean_y'], spec['mean_z']
    iu, ju = np.triu_indices(len(f_dom), k=1)
    dist = np.sqrt((mx_[iu]-mx_[ju])**2 + (my_[iu]-my_[ju])**2 + (mz_[iu]-mz_[ju])**2)
    dfp = np.abs(f_dom[iu] - f_dom[ju])
    print(f"  frequency coherence: median |Δf| vs pair distance")
    for lo, hi in [(0, 1.2), (1.2, 2.4), (2.4, 4.0), (4.0, 6.0), (6.0, 9.0), (9.0, 25.0)]:
        m = (dist >= lo) & (dist < hi)
        if m.any():
            print(f"    dist=[{lo:4.1f},{hi:4.1f}):  n={int(m.sum()):6d}  median_|Δf|={np.median(dfp[m]):.4f}  mean={np.mean(dfp[m]):.4f}")

    # ASCII spectrum
    print(f"\n  Ensemble PSD (log power):")
    ascii_spectrum(spec['freqs'], spec['ens_psd'], f_max=2.0, n_bars=60, height=12)


def main():
    paths = {
        'DISSIPATIVE (baseline)': sys.argv[1] if len(sys.argv) > 1 else 'traj.bin',
        'CONSERVATIVE (new)':     sys.argv[2] if len(sys.argv) > 2 else 'traj_cons.bin',
    }
    specs = {}
    for name, p in paths.items():
        if not Path(p).exists():
            print(f"Skipping {name}: {p} not found")
            continue
        traj = load_traj(p)
        spec = per_particle_spectrum(traj, window_frac=(0.25, 1.0))
        specs[name] = spec
        report(name, spec)

    # Comparative TSV
    if len(specs) >= 2:
        # Interpolate onto common frequency grid if needed (they should match)
        names = list(specs.keys())
        out = Path(paths[names[0]]).parent
        with open(out / 'psd_comparison.tsv', 'w') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['freq', 'omega'] + [f'psd_{n}' for n in names])
            freqs = specs[names[0]]['freqs']
            for i in range(len(freqs)):
                row = [f"{freqs[i]:.5f}", f"{2*np.pi*freqs[i]:.5f}"]
                for n in names:
                    if i < len(specs[n]['ens_psd']):
                        row.append(f"{specs[n]['ens_psd'][i]:.5g}")
                    else:
                        row.append('nan')
                w.writerow(row)
        print(f"\nWrote {out}/psd_comparison.tsv")

        # Plot
        try:
            import matplotlib; matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(2, 2, figsize=(13, 10))
            colors = {'DISSIPATIVE (baseline)': '#c44', 'CONSERVATIVE (new)': '#2a7'}
            for ax, (title, lo, hi) in zip([axes[0,0]], [('Ensemble PSD (log y)', 0, 3.0)]):
                for n, s in specs.items():
                    m = (s['freqs'] > 0) & (s['freqs'] < hi)
                    ax.semilogy(s['freqs'][m], s['ens_psd'][m], label=n, color=colors.get(n, 'k'), lw=0.8)
                ax.set_xlabel('frequency (1/t)'); ax.set_ylabel('PSD')
                ax.set_title(title); ax.legend(); ax.grid(alpha=0.3)

            ax = axes[0,1]
            for n, s in specs.items():
                ax.hist(s['freq_dom'], bins=40, range=(0, max(s['freq_dom'].max()*1.2, 0.5)),
                        alpha=0.6, label=n, color=colors.get(n, 'k'))
            ax.set_xlabel('dominant freq (1/t)'); ax.set_ylabel('# particles')
            ax.set_title('Dominant frequency distribution'); ax.legend(); ax.grid(alpha=0.3)

            ax = axes[1,0]
            for n, s in specs.items():
                ax.scatter(s['mean_nn'], s['freq_dom'], s=10, alpha=0.6,
                           label=n, color=colors.get(n, 'k'))
            ax.set_xlabel('mean neighbor count'); ax.set_ylabel('dominant freq')
            ax.set_title('freq vs neighbor count'); ax.legend(); ax.grid(alpha=0.3)

            ax = axes[1,1]
            for n, s in specs.items():
                f_dom = s['freq_dom']
                iu, ju = np.triu_indices(len(f_dom), k=1)
                d = np.sqrt((s['mean_x'][iu]-s['mean_x'][ju])**2
                          + (s['mean_y'][iu]-s['mean_y'][ju])**2
                          + (s['mean_z'][iu]-s['mean_z'][ju])**2)
                df = np.abs(f_dom[iu] - f_dom[ju])
                # Bin by distance and plot median df
                edges = np.linspace(0, 20, 20)
                centers = 0.5*(edges[:-1]+edges[1:])
                meds = np.array([np.median(df[(d>=edges[i])&(d<edges[i+1])]) if ((d>=edges[i])&(d<edges[i+1])).any() else np.nan
                                 for i in range(len(edges)-1)])
                ax.plot(centers, meds, label=n, color=colors.get(n, 'k'), marker='o', lw=1.5)
            ax.set_xlabel('pair distance'); ax.set_ylabel('median |Δf|')
            ax.set_title('Frequency coherence vs distance'); ax.legend(); ax.grid(alpha=0.3)

            fig.tight_layout()
            out_png = out / 'freq_compare.png'
            fig.savefig(out_png, dpi=120)
            print(f"Wrote {out_png}")
        except Exception as e:
            print(f"(plot skipped: {e})")


if __name__ == '__main__':
    main()
