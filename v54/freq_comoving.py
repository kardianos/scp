#!/usr/bin/env python3
"""Co-moving-frame FFT analysis: subtract centroid + local neighbor drift
before FFT, so the residual is per-particle oscillation about its local
equilibrium, not bulk expansion of the whole system.
"""
import sys
import numpy as np
from pathlib import Path
import csv
from freq_analysis import load_traj


def comoving_velocity(P, k_local=12):
    """Return residual velocity = particle_v − k-nearest-neighbor mean_v
    computed per sample. Uses mean positions to identify neighbors once.
    """
    n_samples, n_tracked = P.shape
    mx = P['x'].mean(axis=0); my = P['y'].mean(axis=0); mz = P['z'].mean(axis=0)
    # k-nearest neighbors among tracked particles (stable set, by mean position)
    # Not necessarily physical contacts, but captures local bulk drift.
    pos = np.stack([mx, my, mz], axis=1)  # (n_tracked, 3)
    # pairwise distance
    diff = pos[:, None, :] - pos[None, :, :]
    dist = np.sqrt((diff**2).sum(axis=-1))
    np.fill_diagonal(dist, np.inf)
    knn = np.argsort(dist, axis=1)[:, :k_local]  # (n_tracked, k_local)

    vx = P['vx'].astype(np.float64)
    vy = P['vy'].astype(np.float64)
    vz = P['vz'].astype(np.float64)

    # For each particle, at each sample, subtract mean v of its k_local neighbors
    vx_res = np.empty_like(vx); vy_res = np.empty_like(vy); vz_res = np.empty_like(vz)
    for k in range(n_tracked):
        nbrs = knn[k]
        vx_res[:, k] = vx[:, k] - vx[:, nbrs].mean(axis=1)
        vy_res[:, k] = vy[:, k] - vy[:, nbrs].mean(axis=1)
        vz_res[:, k] = vz[:, k] - vz[:, nbrs].mean(axis=1)
    vmag_res = np.sqrt(vx_res**2 + vy_res**2 + vz_res**2)
    return vmag_res


def fft_per_particle(vmag, dt_s, window_frac=(0.25, 1.0)):
    n_samples, n_tracked = vmag.shape
    i0 = int(window_frac[0]*n_samples); i1 = int(window_frac[1]*n_samples)
    seg = vmag[i0:i1].copy()
    seg -= seg.mean(axis=0, keepdims=True)
    # Linear detrend per particle
    L = seg.shape[0]; x = np.arange(L)
    for k in range(n_tracked):
        s, i = np.polyfit(x, seg[:, k], 1)
        seg[:, k] -= s*x + i
    w = np.hanning(L)[:, None]
    F = np.fft.rfft(seg*w, axis=0)
    freqs = np.fft.rfftfreq(L, d=dt_s)
    Pspec = np.abs(F)**2
    dom = np.argmax(Pspec[1:], axis=0) + 1
    return freqs, Pspec, freqs[dom]


def ascii_spectrum(freqs, psd, f_max=2.0, n_bars=60, height=12):
    mask = (freqs > 0) & (freqs < f_max)
    f = freqs[mask]; p = psd[mask]
    bin_edges = np.linspace(0, f_max, n_bars + 1)
    binned = np.zeros(n_bars)
    for i in range(n_bars):
        m = (f >= bin_edges[i]) & (f < bin_edges[i+1])
        if m.any():
            binned[i] = p[m].max()
    log_p = np.log10(binned + 1e-20)
    lo, hi = log_p.min(), log_p.max()
    if hi <= lo: hi = lo + 1
    norm = (log_p - lo) / (hi - lo)
    for r in range(height, 0, -1):
        line = ''
        for i in range(n_bars):
            line += '#' if norm[i] * height >= r else ' '
        print(f'  {r:2d} | {line}')
    print(f'     +' + '-' * n_bars)
    print(f'     f=0{" " * (n_bars//2 - 5)}f={f_max/2:.2f}{" " * (n_bars//2 - 5)}f={f_max:.2f}')
    print(f'  log10 power: max={hi:.2f}, min={lo:.2f}')


def report(name, traj, spec):
    P = traj['samples']['p']
    mean_nn = P['nn'].astype(float).mean(axis=0)
    mean_x = P['x'].mean(axis=0); mean_y = P['y'].mean(axis=0); mean_z = P['z'].mean(axis=0)
    freqs, Pspec, f_dom = spec
    n_tracked = len(f_dom)

    print(f"\n===== {name} — CO-MOVING FRAME =====")
    print(f"  f_dom: median={np.median(f_dom):.4f}  mean={np.mean(f_dom):.4f}  std={np.std(f_dom):.4f}")
    print(f"         range=[{f_dom.min():.4f}, {f_dom.max():.4f}]   unique bins={len(np.unique(f_dom))}")

    r_nn = np.corrcoef(mean_nn, f_dom)[0,1]
    r_sqrt = np.corrcoef(np.sqrt(mean_nn+1e-9), f_dom)[0,1]
    print(f"  Pearson r(nn, f_dom)       = {r_nn:+.4f}")
    print(f"  Pearson r(sqrt(nn), f_dom) = {r_sqrt:+.4f}")

    print(f"  stratified by mean neighbor count:")
    for lo, hi in [(0,10),(10,14),(14,18),(18,22),(22,26),(26,32)]:
        m = (mean_nn >= lo) & (mean_nn < hi)
        if m.any():
            print(f"    nn=[{lo},{hi}):  n={int(m.sum()):3d}  median_f={np.median(f_dom[m]):.4f}  mean_f={np.mean(f_dom[m]):.4f}  std={np.std(f_dom[m]):.4f}")

    # Histogram
    nb = 30
    hc, he = np.histogram(f_dom, bins=nb, range=(0, max(f_dom.max()*1.2, 0.5)))
    mx = hc.max()
    print(f"  Histogram (bin {he[1]-he[0]:.4f}):")
    for i in range(nb):
        if hc[i] > 0:
            bar = '#' * int(40*hc[i]/max(mx,1))
            print(f"    f=[{he[i]:5.3f},{he[i+1]:5.3f})  n={hc[i]:3d}  {bar}")

    # Spatial coherence
    iu, ju = np.triu_indices(n_tracked, k=1)
    d = np.sqrt((mean_x[iu]-mean_x[ju])**2 + (mean_y[iu]-mean_y[ju])**2 + (mean_z[iu]-mean_z[ju])**2)
    df = np.abs(f_dom[iu] - f_dom[ju])
    print(f"  Frequency coherence (|Δf| vs pair distance):")
    for lo, hi in [(0, 1.2), (1.2, 2.4), (2.4, 4.0), (4.0, 6.0), (6.0, 9.0), (9.0, 25.0)]:
        m = (d >= lo) & (d < hi)
        if m.any():
            print(f"    dist=[{lo:4.1f},{hi:4.1f}):  n={int(m.sum()):6d}  median_|Δf|={np.median(df[m]):.4f}  mean={np.mean(df[m]):.4f}")

    print(f"\n  Ensemble PSD (log):")
    ascii_spectrum(freqs, Pspec.mean(axis=1), f_max=2.0, n_bars=60, height=12)


def main():
    paths = {
        'DISSIPATIVE': sys.argv[1] if len(sys.argv) > 1 else 'traj.bin',
        'CONSERVATIVE': sys.argv[2] if len(sys.argv) > 2 else 'traj_cons.bin',
    }
    specs = {}
    for name, p in paths.items():
        if not Path(p).exists():
            continue
        traj = load_traj(p)
        P = traj['samples']['p']
        dt_s = traj['dt'] * traj['sps']
        print(f"\nLoading {name}: {p}  (N={traj['n_total']}, n_tracked={traj['n_tracked']}, T={traj['n_samples']*dt_s:.1f})")
        vres = comoving_velocity(P, k_local=12)
        print(f"  co-moving |v_res|: mean={vres.mean():.4f}, max={vres.max():.4f}")
        spec = fft_per_particle(vres, dt_s, window_frac=(0.25, 1.0))
        specs[name] = (traj, spec)
        report(name, traj, spec)

    # Save comparison
    if specs:
        out = Path(list(paths.values())[0]).parent
        with open(out/'freq_comoving.tsv', 'w') as f:
            w = csv.writer(f, delimiter='\t')
            names = list(specs.keys())
            w.writerow(['pid','mean_nn'] + [f'f_dom_{n}' for n in names])
            for k in range(specs[names[0]][0]['n_tracked']):
                row = [int(specs[names[0]][0]['ids'][k]),
                       f"{specs[names[0]][0]['samples']['p']['nn'].astype(float).mean(axis=0)[k]:.3f}"]
                for n in names:
                    traj, spec = specs[n]
                    row.append(f"{spec[2][k]:.5f}")
                w.writerow(row)
        print(f"\nWrote {out}/freq_comoving.tsv")


if __name__ == '__main__':
    main()
