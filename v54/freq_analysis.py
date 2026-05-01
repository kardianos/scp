#!/usr/bin/env python3
"""FFT analysis of point-medium trajectories from v54/points.c.

Reads traj.bin, computes per-particle velocity-magnitude FFT, extracts
dominant frequencies, stratifies by neighbor count, and writes TSV
primary-source outputs plus PNG plots.
"""
import struct
import sys
import csv
from pathlib import Path
import numpy as np


def load_traj(path):
    raw = Path(path).read_bytes()
    off = 0
    magic = raw[off:off+4]; off += 4
    if magic != b'PNTS':
        raise ValueError(f"Bad magic: {magic!r}")
    ver, n_total, n_tracked, n_samples, sps = struct.unpack_from('<IIIII', raw, off); off += 20
    dt, r_touch, L = struct.unpack_from('<ddd', raw, off); off += 24
    ids = np.frombuffer(raw, dtype='<u4', count=n_tracked, offset=off).copy(); off += 4 * n_tracked

    particle_dt = np.dtype([
        ('x',  '<f4'), ('y',  '<f4'), ('z',  '<f4'),
        ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
        ('nn', 'u1'),  ('pad', 'V3'),
    ])
    sample_dt = np.dtype([
        ('step', '<u4'), ('pad', '<u4'), ('t', '<f8'),
        ('p', particle_dt, (n_tracked,)),
    ])
    assert sample_dt.itemsize == 16 + 28 * n_tracked, sample_dt.itemsize

    samples = np.frombuffer(raw, dtype=sample_dt, count=n_samples, offset=off)
    return {
        'ver': ver, 'n_total': n_total, 'n_tracked': n_tracked,
        'n_samples': n_samples, 'sps': sps, 'dt': dt,
        'r_touch': r_touch, 'L': L,
        'ids': ids, 'samples': samples,
    }


def analyze(traj):
    S = traj['samples']
    P = S['p']  # shape (n_samples, n_tracked)
    dt_sample = traj['dt'] * traj['sps']
    n_samples, n_tracked = P.shape

    # Velocity-magnitude per-particle time series (n_samples, n_tracked)
    vx = P['vx'].astype(np.float64)
    vy = P['vy'].astype(np.float64)
    vz = P['vz'].astype(np.float64)
    vmag = np.sqrt(vx*vx + vy*vy + vz*vz)

    # Remove DC, apply Hann window
    v = vmag - vmag.mean(axis=0, keepdims=True)
    window = np.hanning(n_samples)[:, None]
    vw = v * window

    # FFT per particle — frequencies are in 1/time (code units)
    F = np.fft.rfft(vw, axis=0)
    freqs = np.fft.rfftfreq(n_samples, d=dt_sample)
    Pspec = (np.abs(F) ** 2)  # power spectrum, shape (n_freqs, n_tracked)

    # Dominant frequency per particle (skip DC bin 0)
    dom_idx = np.argmax(Pspec[1:], axis=0) + 1
    freq_dom = freqs[dom_idx]
    omega_dom = 2 * np.pi * freq_dom

    # Mean neighbor count per particle (over all samples)
    mean_nn = P['nn'].astype(float).mean(axis=0)

    # Mean particle position over run
    mean_x = P['x'].mean(axis=0)
    mean_y = P['y'].mean(axis=0)
    mean_z = P['z'].mean(axis=0)
    mean_r = np.sqrt(mean_x**2 + mean_y**2 + mean_z**2)

    # Ensemble-averaged PSD
    ens_psd = Pspec.mean(axis=1)

    return {
        'freqs': freqs, 'ens_psd': ens_psd,
        'freq_dom': freq_dom, 'omega_dom': omega_dom,
        'mean_nn': mean_nn,
        'mean_x': mean_x, 'mean_y': mean_y, 'mean_z': mean_z, 'mean_r': mean_r,
        'Pspec': Pspec, 'dt_sample': dt_sample, 'T_total': dt_sample * n_samples,
    }


def write_tsv(path, rows, header):
    with open(path, 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'traj.bin'
    traj = load_traj(path)
    n_samp = traj['n_samples']
    n_track = traj['n_tracked']
    dt_s = traj['dt'] * traj['sps']
    T = dt_s * n_samp

    print(f"Loaded traj.bin:")
    print(f"  n_total       = {traj['n_total']}")
    print(f"  n_tracked     = {n_track}")
    print(f"  n_samples     = {n_samp}")
    print(f"  dt            = {traj['dt']}")
    print(f"  steps_per_samp= {traj['sps']}")
    print(f"  dt_sample     = {dt_s}")
    print(f"  T_total       = {T}")
    print(f"  Nyquist       = {1/(2*dt_s):.3f}  (freq),  {np.pi/dt_s:.3f} (omega)")
    print(f"  df            = {1/T:.5f}")
    print(f"  r_touch       = {traj['r_touch']}   L = {traj['L']}")

    R = analyze(traj)
    f_dom = R['freq_dom']
    om_dom = R['omega_dom']
    nn = R['mean_nn']

    print()
    print(f"Dominant-frequency statistics (per-particle; n={n_track}):")
    print(f"  median f    : {np.median(f_dom):.4f}   omega={2*np.pi*np.median(f_dom):.3f}")
    print(f"  mean   f    : {np.mean(f_dom):.4f}    std={np.std(f_dom):.4f}")
    print(f"  range  f    : [{np.min(f_dom):.4f}, {np.max(f_dom):.4f}]")
    print(f"  unique bins : {len(np.unique(f_dom))}")

    r_nn_f = np.corrcoef(nn, f_dom)[0, 1]
    print(f"\nPearson r(mean_neighbors, f_dom) = {r_nn_f:.4f}")
    r_nn_sqrt = np.corrcoef(np.sqrt(nn + 1e-9), f_dom)[0, 1]
    print(f"Pearson r(sqrt(mean_neighbors), f_dom) = {r_nn_sqrt:.4f}")

    print("\nStratified by mean neighbor count:")
    bins = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 28), (28, 40)]
    print(f"  {'nn_bin':>8}  {'count':>5}  {'med_f':>8}  {'std_f':>8}  {'med_om':>8}")
    for lo, hi in bins:
        mask = (nn >= lo) & (nn < hi)
        c = int(mask.sum())
        if c == 0:
            continue
        mf = np.median(f_dom[mask])
        sf = np.std(f_dom[mask])
        print(f"  [{lo:2d},{hi:2d}) {c:5d}  {mf:8.4f}  {sf:8.4f}  {2*np.pi*mf:8.3f}")

    # Histogram of f_dom
    n_bins = 40
    h_counts, h_edges = np.histogram(f_dom, bins=n_bins)
    print(f"\nHistogram of dominant frequency  (bin width = {(h_edges[1]-h_edges[0]):.4f}):")
    mx = h_counts.max()
    for i in range(n_bins):
        if h_counts[i] == 0 and i > 0 and h_counts[i-1] == 0:
            continue
        bar = '#' * int(60 * h_counts[i] / max(mx, 1))
        print(f"  f=[{h_edges[i]:5.3f},{h_edges[i+1]:5.3f})  n={h_counts[i]:4d} {bar}")

    # --- Frequency coherence: do spatially-close particles share frequencies? ---
    # For each pair of tracked particles, compute |Δf| and their mean separation.
    # If close particles lock in frequency, |Δf| should be small at small distance
    # and grow with distance (or saturate at the population std).
    print("\nFrequency coherence vs distance:")
    mx = R['mean_x']; my = R['mean_y']; mz = R['mean_z']
    nt = n_track
    # Distances: (nt, nt) — only upper triangle
    iu, ju = np.triu_indices(nt, k=1)
    dx = mx[iu] - mx[ju]
    dy = my[iu] - my[ju]
    dz = mz[iu] - mz[ju]
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    dfp = np.abs(f_dom[iu] - f_dom[ju])
    # Bin by distance and report median |Δf|
    dbin_edges = [0.0, 1.2, 2.4, 4.0, 6.0, 9.0, 14.0, 25.0]
    print(f"  {'dist_bin':>12}  {'pairs':>7}  {'med_|df|':>10}  {'mean_|df|':>10}")
    strat_coh = []
    for i in range(len(dbin_edges) - 1):
        lo, hi = dbin_edges[i], dbin_edges[i+1]
        m = (dist >= lo) & (dist < hi)
        c = int(m.sum())
        if c == 0:
            continue
        med = np.median(dfp[m])
        mean = np.mean(dfp[m])
        print(f"  [{lo:4.1f},{hi:4.1f})  {c:7d}  {med:10.5f}  {mean:10.5f}")
        strat_coh.append([f"{lo:.1f}", f"{hi:.1f}", c, f"{med:.5f}", f"{mean:.5f}"])
    # Save pair TSV (subsample if too large)
    out_dir = Path(path).parent
    write_tsv(out_dir/'freq_coherence_vs_distance.tsv', strat_coh,
              ['dist_lo', 'dist_hi', 'n_pairs', 'median_abs_df', 'mean_abs_df'])

    # --- Primary-source TSVs ---
    out_dir = Path(path).parent
    rows = []
    for k in range(n_track):
        rows.append([
            int(traj['ids'][k]),
            f"{R['mean_x'][k]:.4f}", f"{R['mean_y'][k]:.4f}", f"{R['mean_z'][k]:.4f}",
            f"{R['mean_r'][k]:.4f}", f"{R['mean_nn'][k]:.3f}",
            f"{R['freq_dom'][k]:.5f}", f"{R['omega_dom'][k]:.5f}",
        ])
    write_tsv(out_dir/'freq_per_particle.tsv', rows,
              ['pid','mean_x','mean_y','mean_z','mean_r','mean_nn','freq_dom','omega_dom'])

    rows_psd = [[f"{R['freqs'][i]:.5f}", f"{2*np.pi*R['freqs'][i]:.5f}", f"{R['ens_psd'][i]:.5g}"]
                for i in range(len(R['freqs']))]
    write_tsv(out_dir/'ensemble_psd.tsv', rows_psd, ['freq','omega','power'])

    # Stratified PSD for low/mid/high nn groups
    strat = {
        'low':  nn < 10,
        'mid':  (nn >= 10) & (nn < 18),
        'high': nn >= 18,
    }
    header = ['freq', 'omega'] + [f'psd_{k}' for k in strat]
    psd_strat = {k: R['Pspec'][:, m].mean(axis=1) if m.any() else np.zeros_like(R['freqs'])
                 for k, m in strat.items()}
    rows_s = []
    for i in range(len(R['freqs'])):
        rows_s.append([f"{R['freqs'][i]:.5f}", f"{2*np.pi*R['freqs'][i]:.5f}"] +
                      [f"{psd_strat[k][i]:.5g}" for k in strat])
    write_tsv(out_dir/'ensemble_psd_stratified.tsv', rows_s, header)

    print(f"\nWrote:")
    print(f"  {out_dir}/freq_per_particle.tsv")
    print(f"  {out_dir}/ensemble_psd.tsv")
    print(f"  {out_dir}/ensemble_psd_stratified.tsv")

    # Plots
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        ax = axes[0, 0]
        ax.hist(f_dom, bins=40, color='#2a7ab9', edgecolor='k')
        ax.set_xlabel('dominant frequency (1/t)')
        ax.set_ylabel('# particles')
        ax.set_title(f'Per-particle dominant freq (n={n_track}, T={T:.0f})')
        ax.grid(alpha=0.3)

        ax = axes[0, 1]
        fmask = R['freqs'] > 0
        ax.semilogy(R['freqs'][fmask], R['ens_psd'][fmask], color='k', label='ensemble')
        for k in strat:
            m = strat[k]
            if m.any():
                ax.semilogy(R['freqs'][fmask], psd_strat[k][fmask], label=f'{k} nn (n={int(m.sum())})', alpha=0.75)
        ax.set_xlabel('frequency (1/t)')
        ax.set_ylabel('power')
        ax.set_title('Ensemble-averaged PSD (stratified by mean_nn)')
        ax.legend()
        ax.grid(alpha=0.3)

        ax = axes[1, 0]
        sc = ax.scatter(nn, f_dom, c=R['mean_r'], s=15, cmap='viridis', alpha=0.8)
        ax.set_xlabel('mean neighbor count')
        ax.set_ylabel('dominant frequency (1/t)')
        ax.set_title(f'ω vs neighbor count  (Pearson r={r_nn_f:.3f})')
        plt.colorbar(sc, ax=ax, label='mean |r|')
        ax.grid(alpha=0.3)

        ax = axes[1, 1]
        sc = ax.scatter(R['mean_x'], R['mean_y'],
                        c=f_dom, s=15, cmap='plasma')
        ax.set_xlabel('mean x')
        ax.set_ylabel('mean y')
        ax.set_title('Spatial map colored by dominant frequency')
        plt.colorbar(sc, ax=ax, label='f_dom')
        ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig(out_dir / 'freq_analysis.png', dpi=120)
        print(f"  {out_dir}/freq_analysis.png")
    except Exception as e:
        print(f"(plotting skipped: {e})")


if __name__ == '__main__':
    main()
