#!/usr/bin/env python3
"""
viz_particle.py — Binary snapshot reader + visualization for proton.c

Reads binary snapshots (N, L, step, time, Q4 field data) and generates
2x2 panel figures showing xz-plane slices:
  1. Scalar field s(x,z)
  2. Pion magnitude |pi|(x,z)
  3. Energy density (from numerical derivatives)
  4. Pion arrows quiver plot over energy background

Usage:
  python viz_particle.py snapshot_0000.bin
  python viz_particle.py snapshot_*.bin          # batch mode
  python viz_particle.py -tseries timeseries.dat # time series plot
"""

import sys
import struct
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def read_snapshot(fname):
    """Read binary snapshot from proton.c."""
    with open(fname, 'rb') as f:
        N = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        step = struct.unpack('i', f.read(4))[0]
        time = struct.unpack('d', f.read(8))[0]
        data = np.fromfile(f, dtype=np.float64, count=N**3 * 4)
        q = data.reshape(N, N, N, 4)  # (s, f1, f2, f3)
    return N, L, step, time, q


def compute_derivatives_slice(q, N, L, j_slice):
    """Compute field derivatives on xz-plane at y=j_slice.
    Returns dq_dx[i,k,comp], dq_dy[i,k,comp], dq_dz[i,k,comp]."""
    h = 2.0 * L / N
    inv2h = 1.0 / (2.0 * h)
    dqdx = np.zeros((N, N, 4))
    dqdy = np.zeros((N, N, 4))
    dqdz = np.zeros((N, N, 4))

    for i in range(1, N-1):
        for k in range(1, N-1):
            for comp in range(4):
                dqdx[i, k, comp] = (q[i+1, j_slice, k, comp] - q[i-1, j_slice, k, comp]) * inv2h
                dqdy[i, k, comp] = (q[i, min(j_slice+1, N-1), k, comp] -
                                     q[i, max(j_slice-1, 0), k, comp]) * inv2h
                dqdz[i, k, comp] = (q[i, j_slice, k+1, comp] - q[i, j_slice, k-1, comp]) * inv2h
    return dqdx, dqdy, dqdz


def compute_energy_density_slice(q, N, L, j_slice):
    """Compute E2 energy density on xz-plane at y=j_slice."""
    dqdx, dqdy, dqdz = compute_derivatives_slice(q, N, L, j_slice)
    edens = 0.5 * (np.sum(dqdx**2, axis=-1) + np.sum(dqdy**2, axis=-1) + np.sum(dqdz**2, axis=-1))
    return edens


def q4_mul(a, b):
    """Quaternion multiply: a*b, where a,b have shape (..., 4)."""
    s = a[...,0]*b[...,0] - a[...,1]*b[...,1] - a[...,2]*b[...,2] - a[...,3]*b[...,3]
    f1 = a[...,0]*b[...,1] + a[...,1]*b[...,0] - a[...,2]*b[...,3] + a[...,3]*b[...,2]
    f2 = a[...,0]*b[...,2] + a[...,1]*b[...,3] + a[...,2]*b[...,0] - a[...,3]*b[...,1]
    f3 = a[...,0]*b[...,3] - a[...,1]*b[...,2] + a[...,2]*b[...,1] + a[...,3]*b[...,0]
    return np.stack([s, f1, f2, f3], axis=-1)


def compute_baryon_density_slice(q, N, L, j_slice):
    """Compute baryon density B^0 on xz-plane.
    B^0 = -(1/(2*pi^2)) * epsilon_{ijk} * <(q^-1 dq/dx_i)(q^-1 dq/dx_j)(q^-1 dq/dx_k)>_0
    For unit quaternion q, q^-1 = q_rev = (s, -f1, -f2, -f3)."""
    dqdx, dqdy, dqdz = compute_derivatives_slice(q, N, L, j_slice)
    bdens = np.zeros((N, N))

    for i in range(1, N-1):
        for k in range(1, N-1):
            qr = np.array([q[i, j_slice, k, 0], -q[i, j_slice, k, 1],
                           -q[i, j_slice, k, 2], -q[i, j_slice, k, 3]])
            n2 = np.dot(q[i, j_slice, k], q[i, j_slice, k])
            if n2 < 1e-20:
                continue
            qr /= n2  # q^{-1}

            # Left currents: A_d = q^{-1} * dq/dx_d
            def qmul(a, b):
                return np.array([
                    a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                    a[0]*b[1]+a[1]*b[0]-a[2]*b[3]+a[3]*b[2],
                    a[0]*b[2]+a[1]*b[3]+a[2]*b[0]-a[3]*b[1],
                    a[0]*b[3]-a[1]*b[2]+a[2]*b[1]+a[3]*b[0]])

            Ax = qmul(qr, dqdx[i, k])
            Ay = qmul(qr, dqdy[i, k])
            Az = qmul(qr, dqdz[i, k])

            # B^0 = -(1/(2*pi^2)) * epsilon_{ijk} * (Ai x Aj) . Ak (quaternion triple product)
            # = -(1/(2*pi^2)) * [Ax, Ay, Az] where [A,B,C] = <A*B*C>_0 (scalar part of triple product)
            # Using epsilon: xyz - xzy + yzx - yxz + zxy - zyx = 2*(xyz - xzy + yzx)
            # Simpler: B^0 = -(1/(2*pi^2)) * sum over cyclic perms of epsilon
            AxAy = qmul(Ax, Ay)
            AxAz = qmul(Ax, Az)
            AyAz = qmul(Ay, Az)

            # <Ax*Ay*Az>_0 = dot(AxAy, Az_rev) but for pure imaginaries...
            # Actually: <ABC>_0 = scalar part of A*B*C
            AxAyAz = qmul(AxAy, Az)
            AyAxAz = qmul(qmul(Ay, Ax), Az)

            # B^0 = -(1/(2*pi^2)) * (AxAyAz[0] - AyAxAz[0])  ... but need full epsilon sum
            # Full: eps_123*A1*A2*A3 = A_x*A_y*A_z - A_x*A_z*A_y + A_y*A_z*A_x - A_y*A_x*A_z + A_z*A_x*A_y - A_z*A_y*A_x
            # = [Ax,Ay]*Az + [Ay,Az]*Ax + [Az,Ax]*Ay  (where [A,B] = AB - BA)
            # Scalar part = 2*(AxAyAz[0] - AyAxAz[0])... no, let me be careful.

            # epsilon_{ijk} A_i A_j A_k = 6 * A_x A_y A_z (antisymmetric sum over 3!)
            # Actually epsilon sums to: A1A2A3 - A1A3A2 + A2A3A1 - A2A1A3 + A3A1A2 - A3A2A1
            AxAzAy = qmul(AxAz, Ay)
            # By quaternion identities for pure imaginary (trace-free) left currents:
            # The scalar part of [A,B]*C = <(AB-BA)C>_0 = <ABC>_0 - <BAC>_0
            # Total = <AxAyAz - AxAzAy + AyAzAx - AyAxAz + AzAxAy - AzAyAx>_0
            # = 2*(<AxAyAz>_0 - <AxAzAy>_0 + <AyAzAx>_0)  (cyclic property)

            # Simplest: just compute scalar part of epsilon sum
            t1 = AxAyAz[0]      # xyz
            t2 = AxAzAy[0]      # xzy
            AyAzAx = qmul(AyAz, Ax)
            t3 = AyAzAx[0]      # yzx
            AyAxAz_s = qmul(qmul(Ay, Ax), Az)[0]  # yxz
            AzAxAy = qmul(qmul(Az, Ax), Ay)
            t5 = AzAxAy[0]      # zxy
            AzAyAx = qmul(qmul(Az, Ay), Ax)
            t6 = AzAyAx[0]      # zyx

            eps_sum = t1 - t2 + t3 - AyAxAz_s + t5 - t6
            bdens[i, k] = -eps_sum / (2.0 * np.pi**2)

    return bdens


def plot_snapshot(fname, outdir=None):
    """Generate 3x2 panel figure for a single snapshot."""
    N, L, step, time, q = read_snapshot(fname)
    h = 2.0 * L / N
    j_mid = N // 2  # y=0 slice

    # Coordinate arrays for xz-plane
    x = np.linspace(-L + 0.5*h, L - 0.5*h, N)
    extent = [-L, L, -L, L]

    # Extract xz-slice at y=0
    s_slice = q[:, j_mid, :, 0]      # scalar field
    f1_slice = q[:, j_mid, :, 1]
    f2_slice = q[:, j_mid, :, 2]
    f3_slice = q[:, j_mid, :, 3]
    pi_mag = np.sqrt(f1_slice**2 + f2_slice**2 + f3_slice**2)

    # Energy density (2nd order, fast)
    edens = compute_energy_density_slice(q, N, L, j_mid)

    # Baryon density on xz-plane
    bdens_xz = compute_baryon_density_slice(q, N, L, j_mid)

    # Also compute baryon density on xy-plane (z=0)
    # We need to swap axes: for xy-plane, "j_slice" is the z-index
    k_mid = N // 2
    # Transpose q to get (x, y, z, comp) -> for xy-plane at z=k_mid
    # We can reuse by permuting: q_xy[i,j] = q[i,j,k_mid]
    q_xy = np.zeros((N, N, N, 4))
    q_xy[:] = q  # same data, but we need z-slice
    # For xy-plane: indices are (x, y) at fixed z=k_mid
    # Reinterpret: treat the k_mid slice of the z-axis
    # compute_baryon_density_slice uses j_slice for the second index (y)
    # For xy-plane, we want derivatives d/dx, d/dy, d/dz at z=k_mid
    # This requires a different slicing approach. Let's just do it inline.
    bdens_xy = np.zeros((N, N))
    h_inv2 = 1.0 / (2.0 * h)
    for i in range(1, N-1):
        for j in range(1, N-1):
            qr = np.array([q[i,j,k_mid,0], -q[i,j,k_mid,1],
                           -q[i,j,k_mid,2], -q[i,j,k_mid,3]])
            n2 = np.dot(q[i,j,k_mid], q[i,j,k_mid])
            if n2 < 1e-20:
                continue
            qr /= n2

            def qmul(a, b):
                return np.array([
                    a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                    a[0]*b[1]+a[1]*b[0]-a[2]*b[3]+a[3]*b[2],
                    a[0]*b[2]+a[1]*b[3]+a[2]*b[0]-a[3]*b[1],
                    a[0]*b[3]-a[1]*b[2]+a[2]*b[1]+a[3]*b[0]])

            dx = (q[i+1,j,k_mid] - q[i-1,j,k_mid]) * h_inv2
            dy = (q[i,j+1,k_mid] - q[i,j-1,k_mid]) * h_inv2
            dz = (q[i,j,min(k_mid+1,N-1)] - q[i,j,max(k_mid-1,0)]) * h_inv2

            Ax = qmul(qr, dx)
            Ay = qmul(qr, dy)
            Az = qmul(qr, dz)

            t1 = qmul(qmul(Ax,Ay),Az)[0]
            t2 = qmul(qmul(Ax,Az),Ay)[0]
            t3 = qmul(qmul(Ay,Az),Ax)[0]
            t4 = qmul(qmul(Ay,Ax),Az)[0]
            t5 = qmul(qmul(Az,Ax),Ay)[0]
            t6 = qmul(qmul(Az,Ay),Ax)[0]
            bdens_xy[i,j] = -(t1 - t2 + t3 - t4 + t5 - t6) / (2.0*np.pi**2)

    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    fig.suptitle(f'Step {step}, t = {time:.3f}  (N={N}, L={L:.1f})',
                 fontsize=14, fontweight='bold')

    # Panel 1: Scalar field s (xz-plane)
    ax = axes[0, 0]
    vmax = max(abs(s_slice.min()), abs(s_slice.max()), 0.1)
    im = ax.imshow(s_slice.T, origin='lower', extent=extent,
                   cmap='RdBu', vmin=-vmax, vmax=vmax, aspect='equal')
    ax.set_title('Scalar field s(x,z) [y=0]')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Panel 2: Pion magnitude (xz-plane)
    ax = axes[0, 1]
    im = ax.imshow(pi_mag.T, origin='lower', extent=extent,
                   cmap='inferno', vmin=0, aspect='equal')
    ax.set_title(r'$|\vec{\pi}|(x,z)$ [y=0]')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Panel 3: Energy density (xz-plane)
    ax = axes[1, 0]
    edens_clip = np.clip(edens, 1e-6, None)
    im = ax.imshow(edens_clip.T, origin='lower', extent=extent,
                   cmap='hot', norm=LogNorm(vmin=max(edens_clip.min(), 1e-4),
                                            vmax=edens_clip.max()),
                   aspect='equal')
    ax.set_title('Energy density [y=0, log]')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Panel 4: Baryon density B^0 (xz-plane)
    ax = axes[1, 1]
    bmax = max(abs(bdens_xz.min()), abs(bdens_xz.max()), 1e-6)
    im = ax.imshow(bdens_xz.T, origin='lower', extent=extent,
                   cmap='RdBu_r', vmin=-bmax, vmax=bmax, aspect='equal')
    ax.set_title(r'Baryon density $B^0$(x,z) [y=0]')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    plt.colorbar(im, ax=ax, shrink=0.8, label=r'$B^0$')

    # Panel 5: Baryon density B^0 (xy-plane)
    ax = axes[2, 0]
    bmax_xy = max(abs(bdens_xy.min()), abs(bdens_xy.max()), 1e-6)
    im = ax.imshow(bdens_xy.T, origin='lower', extent=extent,
                   cmap='RdBu_r', vmin=-bmax_xy, vmax=bmax_xy, aspect='equal')
    ax.set_title(r'Baryon density $B^0$(x,y) [z=0]')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, shrink=0.8, label=r'$B^0$')

    # Panel 6: Pion arrows over energy background (xz-plane)
    ax = axes[2, 1]
    im = ax.imshow(edens_clip.T, origin='lower', extent=extent,
                   cmap='gray_r', norm=LogNorm(vmin=max(edens_clip.min(), 1e-4),
                                               vmax=edens_clip.max()),
                   aspect='equal', alpha=0.5)
    skip = max(1, N // 32)
    xs = x[::skip]
    f1_sub = f1_slice[::skip, ::skip]
    f3_sub = f3_slice[::skip, ::skip]
    X, Z = np.meshgrid(xs, xs)
    ax.quiver(X, Z, f1_sub.T, f3_sub.T, color='blue', alpha=0.7,
              scale=None, headwidth=3)
    ax.set_title(r'$\pi_1, \pi_3$ arrows [y=0]')
    ax.set_xlabel('x'); ax.set_ylabel('z')
    ax.set_xlim(-L, L); ax.set_ylim(-L, L)

    plt.tight_layout()

    if outdir is None:
        outdir = os.path.dirname(fname) or '.'
    base = os.path.splitext(os.path.basename(fname))[0]
    outpath = os.path.join(outdir, f'{base}.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')
    return outpath


def plot_timeseries(fname, outdir=None):
    """Plot time series data from proton.c ASCII output."""
    data = np.loadtxt(fname, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Columns: step, time, E2, E4, E_kin, Q, T_eff
    step = data[:, 0]
    t = data[:, 1]
    E2 = data[:, 2]
    E4 = data[:, 3]
    E_kin = data[:, 4]
    Q = data[:, 5]
    T_eff = data[:, 6]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f'Time series: {os.path.basename(fname)}',
                 fontsize=13, fontweight='bold')

    # Energy
    ax = axes[0, 0]
    ax.plot(t, E2, label='E2', color='blue')
    ax.plot(t, E4, label='E4', color='red')
    ax.plot(t, E_kin, label='E_kin', color='green')
    ax.plot(t, E2+E4+E_kin, label='E_total', color='black', ls='--')
    ax.set_xlabel('time'); ax.set_ylabel('Energy')
    ax.set_title('Energy components')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Topological charge
    ax = axes[0, 1]
    ax.plot(t, Q, color='purple', linewidth=2)
    ax.set_xlabel('time'); ax.set_ylabel('Q')
    ax.set_title('Topological charge')
    ax.axhline(y=round(Q[0]), color='gray', ls='--', alpha=0.5)
    ax.grid(True, alpha=0.3)

    # Effective temperature
    ax = axes[1, 0]
    ax.plot(t, T_eff, color='orange', linewidth=2)
    ax.set_xlabel('time'); ax.set_ylabel('T_eff')
    ax.set_title('Effective temperature')
    ax.grid(True, alpha=0.3)

    # E2/E4 ratio (virial)
    ax = axes[1, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(E4 > 0, E2/E4, np.nan)
    ax.plot(t, ratio, color='brown', linewidth=1.5)
    ax.axhline(y=1.0, color='gray', ls='--', alpha=0.5, label='virial (E2=E4)')
    ax.set_xlabel('time'); ax.set_ylabel('E2/E4')
    ax.set_title('Virial ratio')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if outdir is None:
        outdir = os.path.dirname(fname) or '.'
    base = os.path.splitext(os.path.basename(fname))[0]
    outpath = os.path.join(outdir, f'{base}.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {outpath}')


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == '-tseries':
        if len(sys.argv) < 3:
            print('Usage: viz_particle.py -tseries timeseries.dat')
            sys.exit(1)
        plot_timeseries(sys.argv[2])
    else:
        files = []
        for arg in sys.argv[1:]:
            files.extend(glob.glob(arg))
        files.sort()
        if not files:
            print(f'No files found matching {sys.argv[1:]}')
            sys.exit(1)
        for f in files:
            print(f'Processing {f}...')
            plot_snapshot(f)


if __name__ == '__main__':
    main()
