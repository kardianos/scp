#!/usr/bin/env python3
"""
viz_hopfion.py — Visualize Faddeev-Skyrme hopfion snapshots.

Reads binary snapshot files from hopfion.c and produces:
  1. Midplane slice plots (n1, n2, n3, |pi|, energy density)
  2. 3D isosurface of n3 = 0 (the hopfion core tube)
  3. Time series from timeseries.dat

Usage:
  python viz_hopfion.py data/hopfion/          # all snapshots in dir
  python viz_hopfion.py data/hopfion/snap_00500.bin  # single snapshot
  python viz_hopfion.py --compare data/hopfion_T0 data/hopfion_T1 data/hopfion_T10
"""

import sys
import os
import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_snapshot(fname):
    """Read binary snapshot from hopfion.c."""
    with open(fname, 'rb') as f:
        N = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        step = struct.unpack('i', f.read(4))[0]
        time = struct.unpack('d', f.read(8))[0]
        nc = struct.unpack('i', f.read(4))[0]
        data = np.fromfile(f, dtype=np.float64, count=N**3 * nc)
        field = data.reshape(N, N, N, nc)
    h = L / N
    return N, L, h, step, time, field

def compute_energy_density(n, h):
    """Compute E2 and E4 energy density from field n(x,y,z,3)."""
    N = n.shape[0]
    # Derivatives via central differences
    dn = np.zeros((N, N, N, 3, 3))  # dn[x,y,z,component,direction]
    for d in range(3):
        dn[:,:,:,:,d] = (np.roll(n, -1, axis=d) - np.roll(n, 1, axis=d)) / (2*h)

    # E2 density: (1/2)|∇n|²
    grad_n2 = np.sum(dn**2, axis=(3,4))
    e2 = 0.5 * grad_n2

    # F_{ij} = n · (∂_in × ∂_jn)
    def F(i, j):
        cross = np.cross(dn[:,:,:,:,i], dn[:,:,:,:,j])
        return np.sum(n * cross, axis=-1)

    F01 = F(0, 1)
    F02 = F(0, 2)
    F12 = F(1, 2)
    e4 = 0.25 * (F01**2 + F02**2 + F12**2)

    return e2, e4

def plot_snapshot(fname, outdir=None):
    """Generate visualization panels for a single snapshot."""
    N, L, h, step, time, field = read_snapshot(fname)
    n = field  # (N,N,N,3) = (n1, n2, n3)
    mid = N // 2

    # Coordinate arrays
    x = np.linspace(-L/2 + h/2, L/2 - h/2, N)

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Faddeev-Skyrme hopfion  step={step}  t={time:.3f}  N={N}  L={L:.1f}',
                 fontsize=14)

    # Midplane z=0 (k=mid)
    n_mid = n[:, :, mid, :]  # (N, N, 3) in xz-plane... actually (x,y,3)

    # Panel 1: n3 (scalar, vacuum=-1, core=+1)
    ax = axes[0, 0]
    im = ax.imshow(n_mid[:,:,2].T, extent=[-L/2,L/2,-L/2,L/2],
                   origin='lower', cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_title('$n_3$ (z=0 slice)')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(im, ax=ax)

    # Panel 2: |pi| = sqrt(n1² + n2²) (pion magnitude)
    pi_mag = np.sqrt(n_mid[:,:,0]**2 + n_mid[:,:,1]**2)
    ax = axes[0, 1]
    im = ax.imshow(pi_mag.T, extent=[-L/2,L/2,-L/2,L/2],
                   origin='lower', cmap='inferno', vmin=0, vmax=1)
    ax.set_title('$|\\pi|$ = $\\sqrt{n_1^2+n_2^2}$ (z=0)')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(im, ax=ax)

    # Panel 3: n3 in xz-plane (j=mid)
    n_xz = n[:, mid, :, :]
    ax = axes[0, 2]
    im = ax.imshow(n_xz[:,:,2].T, extent=[-L/2,L/2,-L/2,L/2],
                   origin='lower', cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_title('$n_3$ (y=0 slice)')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax)

    # Panel 4: Energy density (z=0)
    e2, e4 = compute_energy_density(n, h)
    e_mid = (e2 + e4)[:, :, mid]
    ax = axes[1, 0]
    vmax = np.percentile(e_mid, 99)
    im = ax.imshow(e_mid.T, extent=[-L/2,L/2,-L/2,L/2],
                   origin='lower', cmap='hot', vmin=0, vmax=vmax)
    ax.set_title('Energy density (z=0)')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(im, ax=ax)

    # Panel 5: Energy density (y=0)
    e_xz = (e2 + e4)[:, mid, :]
    ax = axes[1, 1]
    im = ax.imshow(e_xz.T, extent=[-L/2,L/2,-L/2,L/2],
                   origin='lower', cmap='hot', vmin=0, vmax=vmax)
    ax.set_title('Energy density (y=0)')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax)

    # Panel 6: Pion direction arrows (z=0, subsampled)
    ax = axes[1, 2]
    step_q = max(1, N // 16)
    X, Y = np.meshgrid(x[::step_q], x[::step_q])
    U = n_mid[::step_q, ::step_q, 0]
    V = n_mid[::step_q, ::step_q, 1]
    C = n_mid[::step_q, ::step_q, 2]
    # Background: |pi|
    ax.imshow(pi_mag.T, extent=[-L/2,L/2,-L/2,L/2],
              origin='lower', cmap='Greys', alpha=0.3, vmin=0, vmax=1)
    ax.quiver(X, Y, U.T, V.T, C.T, cmap='RdBu_r', clim=(-1,1),
              scale=25, width=0.003)
    ax.set_title('Pion direction $(n_1, n_2)$ colored by $n_3$')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_xlim(-L/2, L/2); ax.set_ylim(-L/2, L/2)

    plt.tight_layout()
    if outdir is None:
        outdir = os.path.dirname(fname) or '.'
    base = os.path.splitext(os.path.basename(fname))[0]
    outpath = os.path.join(outdir, f'{base}.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')

def plot_3d_isosurface(fname, outdir=None, threshold=0.0):
    """Plot 3D isosurface of n3 = threshold (hopfion core tube)."""
    from skimage import measure

    N, L, h, step, time, field = read_snapshot(fname)
    n3 = field[:,:,:,2]  # n3 component

    try:
        verts, faces, _, _ = measure.marching_cubes(n3, level=threshold)
    except ValueError:
        print(f'  No isosurface at n3={threshold}')
        return

    # Convert verts from grid indices to physical coords
    verts = verts * h - L/2 + h/2

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
                    cmap='coolwarm', alpha=0.7, edgecolor='none')
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title(f'Hopfion core tube ($n_3$=0)  step={step}  t={time:.3f}')

    # Equal aspect
    maxr = L/2
    ax.set_xlim(-maxr, maxr); ax.set_ylim(-maxr, maxr); ax.set_zlim(-maxr, maxr)

    if outdir is None:
        outdir = os.path.dirname(fname) or '.'
    base = os.path.splitext(os.path.basename(fname))[0]
    outpath = os.path.join(outdir, f'{base}_3d.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')

def plot_timeseries(datadir, outdir=None):
    """Plot time series from timeseries.dat."""
    tsfile = os.path.join(datadir, 'timeseries.dat')
    if not os.path.exists(tsfile):
        print(f'No timeseries.dat in {datadir}')
        return

    data = np.loadtxt(tsfile, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)

    step = data[:, 0]
    t = data[:, 1]
    E2 = data[:, 2]
    E4 = data[:, 3]
    EV = data[:, 4]
    Ekin = data[:, 5]
    Q = data[:, 6]
    Teff = data[:, 7]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Time series: {os.path.basename(datadir)}', fontsize=14)

    ax = axes[0, 0]
    ax.plot(t, E2, label='$E_2$')
    ax.plot(t, E4, label='$E_4$')
    ax.plot(t, Ekin, label='$E_{kin}$')
    ax.plot(t, E2+E4+EV+Ekin, label='$E_{total}$', ls='--', color='k')
    ax.set_xlabel('time'); ax.set_ylabel('Energy')
    ax.legend(); ax.set_title('Energies')

    ax = axes[0, 1]
    ax.plot(t, Q, 'b-', linewidth=2)
    ax.axhline(1.0, color='r', ls='--', alpha=0.5, label='Q=1')
    ax.set_xlabel('time'); ax.set_ylabel('Hopf charge Q')
    ax.set_title('Topological charge')
    ax.legend()

    ax = axes[1, 0]
    ax.plot(t, Teff, 'r-')
    ax.set_xlabel('time'); ax.set_ylabel('$T_{eff}$')
    ax.set_title('Effective temperature')
    ax.set_yscale('log') if np.any(Teff > 0) else None

    ax = axes[1, 1]
    ax.plot(t, E2/np.maximum(E4, 1e-10), 'g-')
    ax.axhline(1.0, color='r', ls='--', alpha=0.5, label='virial $E_2/E_4=1$')
    ax.set_xlabel('time'); ax.set_ylabel('$E_2/E_4$')
    ax.set_title('Virial ratio')
    ax.legend()

    plt.tight_layout()
    if outdir is None:
        outdir = datadir
    outpath = os.path.join(outdir, 'timeseries.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')

def plot_comparison(dirs, outdir=None):
    """Compare time series across multiple runs (different temperatures)."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for datadir in dirs:
        tsfile = os.path.join(datadir, 'timeseries.dat')
        if not os.path.exists(tsfile):
            continue
        data = np.loadtxt(tsfile, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)

        t = data[:, 1]
        E2 = data[:, 2]; E4 = data[:, 3]; EV = data[:, 4]
        Ekin = data[:, 5]; Q = data[:, 6]; Teff = data[:, 7]
        label = os.path.basename(datadir)

        axes[0].plot(t, E2+E4+EV+Ekin, label=label)
        axes[1].plot(t, Q, label=label, linewidth=2)
        axes[2].plot(t, Teff, label=label)

    axes[0].set_xlabel('time'); axes[0].set_ylabel('$E_{total}$')
    axes[0].set_title('Total energy'); axes[0].legend()

    axes[1].set_xlabel('time'); axes[1].set_ylabel('Hopf charge Q')
    axes[1].set_title('Topological charge'); axes[1].legend()
    axes[1].axhline(1.0, color='r', ls='--', alpha=0.3)

    axes[2].set_xlabel('time'); axes[2].set_ylabel('$T_{eff}$')
    axes[2].set_title('Effective temperature'); axes[2].legend()

    plt.tight_layout()
    if outdir is None:
        outdir = os.path.dirname(dirs[0]) or '.'
    outpath = os.path.join(outdir, 'comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')


if __name__ == '__main__':
    import glob

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == '--compare':
        dirs = sys.argv[2:]
        plot_comparison(dirs)
        sys.exit(0)

    paths = sys.argv[1:]
    for path in paths:
        if os.path.isdir(path):
            # Plot timeseries
            plot_timeseries(path)
            # Plot all snapshots
            snaps = sorted(glob.glob(os.path.join(path, '*.bin')))
            for snap in snaps:
                plot_snapshot(snap)
                try:
                    plot_3d_isosurface(snap)
                except ImportError:
                    pass  # skimage not available
        elif path.endswith('.bin'):
            plot_snapshot(path)
            try:
                plot_3d_isosurface(path)
            except ImportError:
                pass
        else:
            print(f'Unknown path: {path}')
