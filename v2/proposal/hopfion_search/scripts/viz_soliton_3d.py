#!/usr/bin/env python3
"""
viz_soliton_3d.py — Visualization of B=1 Skyrmion field structure

Creates publication-quality figures showing:
  1. Single soliton anatomy (energy, baryon density, scalar field, hedgehog arrows)
  2. Two-soliton product ansatz at several separations
  3. 3D isosurface of baryon density with hedgehog arrows

Outputs: results/viz/soliton_*.png
"""

import os
import sys
import numpy as np
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

# ============================================================
# Load 1D profile
# ============================================================

PROFILE = os.path.join(os.path.dirname(__file__),
                       '..', 'data', 'profiles', 'profile_sigma_e1.dat')
data = np.loadtxt(PROFILE)
r_data, f_data, fp_data = data[:, 0], data[:, 1], data[:, 2]

f_of_r  = interp1d(r_data, f_data,  kind='cubic', fill_value=0.0, bounds_error=False)
fp_of_r = interp1d(r_data, fp_data, kind='cubic', fill_value=0.0, bounds_error=False)

OUTDIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'viz')
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Hedgehog field constructors
# ============================================================

def hedgehog_q(x, y, z, cx=0, cy=0, cz=0):
    """Hedgehog quaternion q = (cos f, sin f * r_hat) centered at (cx,cy,cz)."""
    dx, dy, dz = x - cx, y - cy, z - cz
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    rs = np.maximum(r, 1e-12)
    f = f_of_r(r)
    sf = np.sin(f)
    return np.cos(f), sf * dx/rs, sf * dy/rs, sf * dz/rs


def qmul(a, b):
    """Hamilton product of two quaternion 4-tuples (arrays)."""
    return (a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
            a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
            a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
            a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0])


def product_ansatz(x, y, z, D):
    """Product ansatz for two B=1 solitons separated by D along x-axis."""
    q1 = hedgehog_q(x, y, z, cx=-D/2)
    q2 = hedgehog_q(x, y, z, cx=+D/2)
    return qmul(q1, q2)

# ============================================================
# Physics quantities (analytic on hedgehog background)
# ============================================================

def energy_density_hedgehog(R):
    """E2 + E4 density for a single hedgehog at origin (e=1, rho0=1)."""
    rs = np.maximum(R, 1e-12)
    f = f_of_r(R)
    fp = fp_of_r(R)
    sf2 = np.sin(f)**2
    e2 = 0.5 * (fp**2 + 2*sf2/rs**2)
    e4 = fp**2 * sf2/rs**2 + sf2**2 / (2*rs**4)
    return e2 + e4


def baryon_density_hedgehog(R):
    """B^0 for a single hedgehog."""
    rs = np.maximum(R, 1e-12)
    f = f_of_r(R)
    fp = fp_of_r(R)
    return -fp * np.sin(f)**2 / (2 * np.pi**2 * rs**2)

# ============================================================
# Numerical derivatives for general (non-hedgehog) fields
# ============================================================

def compute_baryon_3d(q, h):
    """
    Compute B^0 on a 3D grid using the triple-product formula:
      B0 = -(1/(2pi^2)) <(q_rev * dq_x)(q_rev * dq_y)(q_rev * dq_z)>_0 / |q|^6
    Uses 4th-order central differences.
    """
    inv2pi2 = 1.0 / (2 * np.pi**2)
    s, f1, f2, f3 = q

    # 4th-order central differences: (-1, 8, -8, 1) / (12h)
    def deriv(arr, axis):
        return (-np.roll(arr, -2, axis=axis) + 8*np.roll(arr, -1, axis=axis)
                - 8*np.roll(arr, 1, axis=axis) + np.roll(arr, 2, axis=axis)) / (12*h)

    ds  = [deriv(s, ax) for ax in range(3)]
    df1 = [deriv(f1, ax) for ax in range(3)]
    df2 = [deriv(f2, ax) for ax in range(3)]
    df3 = [deriv(f3, ax) for ax in range(3)]

    # Left currents A_d = q_rev * dq_d  (quaternion product, q_rev = (s, -f1, -f2, -f3))
    def left_current(d):
        a0 = s*ds[d] + f1*df1[d] + f2*df2[d] + f3*df3[d]
        a1 = s*df1[d] - f1*ds[d] - f2*df3[d] + f3*df2[d]
        a2 = s*df2[d] - f2*ds[d] - f3*df1[d] + f1*df3[d]
        a3 = s*df3[d] - f3*ds[d] - f1*df2[d] + f2*df1[d]
        return a0, a1, a2, a3

    A0 = left_current(0)
    A1 = left_current(1)
    A2 = left_current(2)

    # Triple product A0*A1*A2, scalar part
    # First: A01 = A0 * A1
    a01_0 = A0[0]*A1[0] - A0[1]*A1[1] - A0[2]*A1[2] - A0[3]*A1[3]
    a01_1 = A0[0]*A1[1] + A0[1]*A1[0] + A0[2]*A1[3] - A0[3]*A1[2]
    a01_2 = A0[0]*A1[2] - A0[1]*A1[3] + A0[2]*A1[0] + A0[3]*A1[1]
    a01_3 = A0[0]*A1[3] + A0[1]*A1[2] - A0[2]*A1[1] + A0[3]*A1[0]
    # Scalar part of A01 * A2
    scal = a01_0*A2[0] - a01_1*A2[1] - a01_2*A2[2] - a01_3*A2[3]

    norm2 = s**2 + f1**2 + f2**2 + f3**2
    norm6 = norm2**3
    B0 = np.where(norm6 > 1e-30, -inv2pi2 * scal / norm6, 0.0)
    return B0


def compute_energy_3d(q, h):
    """Compute E2 + E4 density on a 3D grid."""
    s, f1, f2, f3 = q

    def deriv(arr, axis):
        return (-np.roll(arr, -2, axis=axis) + 8*np.roll(arr, -1, axis=axis)
                - 8*np.roll(arr, 1, axis=axis) + np.roll(arr, 2, axis=axis)) / (12*h)

    # E2 = (1/2) sum_d |dq/dx_d|^2
    e2 = np.zeros_like(s)
    dq = [[None]*3 for _ in range(4)]
    fields = [s, f1, f2, f3]
    for d in range(3):
        for a in range(4):
            dq[a][d] = deriv(fields[a], d)
        e2 += 0.5 * sum(dq[a][d]**2 for a in range(4))

    # E4: need left currents and commutators
    # A_d = q_rev * dq_d
    def left_current(d):
        a0 = s*dq[0][d] + f1*dq[1][d] + f2*dq[2][d] + f3*dq[3][d]
        a1 = s*dq[1][d] - f1*dq[0][d] - f2*dq[3][d] + f3*dq[2][d]
        a2 = s*dq[2][d] - f2*dq[0][d] - f3*dq[1][d] + f1*dq[3][d]
        a3 = s*dq[3][d] - f3*dq[0][d] - f1*dq[2][d] + f2*dq[1][d]
        return a0, a1, a2, a3

    A = [left_current(d) for d in range(3)]
    norm2 = s**2 + f1**2 + f2**2 + f3**2
    norm4 = norm2**2

    # E4 = -(1/(4e^2)) sum_{d1<d2} [A_{d1}, A_{d2}]^2 / |q|^4
    # [A,B] = A*B - B*A, and [A,B]^2 = scalar part of [A,B]*[A,B]
    e4 = np.zeros_like(s)
    for d1 in range(3):
        for d2 in range(d1+1, 3):
            # comm = A[d1]*A[d2] - A[d2]*A[d1]
            # Since A[d].0 = scalar part, commutator scalar part = 0
            # Only vector parts contribute to commutator
            def qm(a, b):
                return (a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
                        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
                        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
                        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0])

            ab = qm(A[d1], A[d2])
            ba = qm(A[d2], A[d1])
            comm = tuple(ab[i] - ba[i] for i in range(4))
            # comm^2 = comm * comm, scalar part
            comm2_s = comm[0]**2 + comm[1]**2 + comm[2]**2 + comm[3]**2
            # Note: for pure quaternion commutator, [A,B]^2.s = -|[A,B]|^2
            # But comm includes scalar part too. For E4: -(1/4e^2) * comm * comm -> scalar
            # Actually: comm*comm scalar = comm.s^2 - |comm.vec|^2
            comm2_s = (comm[0]**2 - comm[1]**2 - comm[2]**2 - comm[3]**2)
            e4 += np.where(norm4 > 1e-30,
                           -0.25 * comm2_s / norm4, 0.0)  # e=1

    return e2 + e4

# ============================================================
# Figure 1: Single soliton anatomy
# ============================================================

def fig_single_soliton():
    print("Figure 1: Single soliton anatomy...")
    N = 600
    L = 5.0
    x = np.linspace(-L, L, N)
    X, Z = np.meshgrid(x, x)
    R = np.sqrt(X**2 + Z**2)

    e_tot = energy_density_hedgehog(R)
    B0 = baryon_density_hedgehog(R)
    s_field = np.cos(f_of_r(R))
    sinf = np.sin(f_of_r(R))
    R_safe = np.maximum(R, 1e-12)
    pion_x = sinf * X / R_safe
    pion_z = sinf * Z / R_safe

    fig, axes = plt.subplots(2, 2, figsize=(12, 10.5))
    fig.suptitle(r'B = 1 Skyrmion  ($\sigma$-model, e = 1, $\rho_0$ = 1)',
                 fontsize=15, fontweight='bold')

    # (a) Energy density
    ax = axes[0, 0]
    vmax = np.percentile(e_tot, 99.8)
    im = ax.pcolormesh(X, Z, e_tot, cmap='inferno', shading='gouraud',
                       vmin=0, vmax=vmax, rasterized=True)
    ax.set_title(r'(a)  Energy density $\varepsilon_2 + \varepsilon_4$', fontsize=12)
    ax.set_xlabel('x  [code units]')
    ax.set_ylabel('z  [code units]')
    ax.set_aspect('equal')
    cb = plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label(r'$\varepsilon$', fontsize=11)

    # (b) Baryon density
    ax = axes[0, 1]
    vmax_b = np.percentile(B0, 99.8)
    im = ax.pcolormesh(X, Z, B0, cmap='YlOrRd', shading='gouraud',
                       vmin=0, vmax=vmax_b, rasterized=True)
    ax.set_title(r'(b)  Baryon density $B^0$', fontsize=12)
    ax.set_xlabel('x  [code units]')
    ax.set_ylabel('z  [code units]')
    ax.set_aspect('equal')
    cb = plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label(r'$B^0$', fontsize=11)

    # (c) Scalar field
    ax = axes[1, 0]
    im = ax.pcolormesh(X, Z, s_field, cmap='RdBu', shading='gouraud',
                       norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
                       rasterized=True)
    ax.set_title(r'(c)  Scalar field $s = \cos f(r)$', fontsize=12)
    ax.set_xlabel('x  [code units]')
    ax.set_ylabel('z  [code units]')
    ax.set_aspect('equal')
    cb = plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label('s', fontsize=11)
    ax.contour(X, Z, s_field, levels=[0], colors='k', linewidths=0.7, linestyles='--')

    # (d) Hedgehog pion field
    ax = axes[1, 1]
    ax.pcolormesh(X, Z, e_tot, cmap='Greys', shading='gouraud', alpha=0.25,
                  vmin=0, vmax=vmax, rasterized=True)
    step = N // 22
    sl = slice(step//2, None, step)
    pion_mag = np.sqrt(pion_x**2 + pion_z**2)
    Q = ax.quiver(X[sl, sl], Z[sl, sl], pion_x[sl, sl], pion_z[sl, sl],
                  pion_mag[sl, sl], cmap='plasma', scale=20, width=0.005,
                  headwidth=3.5, headlength=4, clim=[0, 1])
    ax.set_title(r'(d)  Pion field $(\pi_1, \pi_3)$ in xz-plane', fontsize=12)
    ax.set_xlabel('x  [code units]')
    ax.set_ylabel('z  [code units]')
    ax.set_aspect('equal')
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    cb = plt.colorbar(Q, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label(r'$|\vec{\pi}|$', fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(OUTDIR, 'soliton_single.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close()

# ============================================================
# Figure 2: Two-soliton product ansatz
# ============================================================

def fig_two_soliton():
    print("Figure 2: Two-soliton product ansatz...")
    N3d = 128
    L = 7.0
    h = 2*L / N3d
    x = np.linspace(-L + h/2, L - h/2, N3d)
    X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')

    separations = [3.0, 5.0, 8.0]
    fig, axes = plt.subplots(2, len(separations), figsize=(5*len(separations), 9.5))
    fig.suptitle(r'Two-soliton product ansatz  $q = q_1 \cdot q_2$  (B = 2)',
                 fontsize=15, fontweight='bold')

    for col, D in enumerate(separations):
        print(f"  D = {D}...")

        # Compute product ansatz on 2D xz-slice only (y=0) — no 3D grid needed
        Xs, Zs = np.meshgrid(x, x, indexing='ij')
        Ys = np.zeros_like(Xs)
        q = product_ansatz(Xs, Ys, Zs, D)

        s_slice = q[0]
        # Pion magnitude: shows "where the soliton is" without derivatives
        pion_mag = np.sqrt(q[1]**2 + q[2]**2 + q[3]**2)

        # Top row: scalar field
        ax = axes[0, col]
        im = ax.pcolormesh(Xs.T, Zs.T, s_slice.T, cmap='RdBu', shading='gouraud',
                           norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
                           rasterized=True)
        ax.contour(Xs.T, Zs.T, s_slice.T, levels=[0], colors='k',
                   linewidths=0.7, linestyles='--')
        ax.set_title(f'D = {D:.0f}:  scalar s(x, z)', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_aspect('equal')
        if col == len(separations) - 1:
            plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label='s')

        # Bottom row: pion field magnitude |F| = sin(f) for hedgehog, varies for product
        ax = axes[1, col]
        im = ax.pcolormesh(Xs.T, Zs.T, pion_mag.T, cmap='inferno', shading='gouraud',
                           vmin=0, vmax=1, rasterized=True)
        ax.contour(Xs.T, Zs.T, pion_mag.T, levels=[0.5, 0.9], colors='w',
                   linewidths=0.6, linestyles='--')
        ax.set_title(f'D = {D:.0f}:  pion magnitude ' + r'$|\vec{\pi}|$', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_aspect('equal')
        if col == len(separations) - 1:
            plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label=r'$|\vec{\pi}|$')

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(OUTDIR, 'soliton_two_body.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close()

# ============================================================
# Figure 3: 3D isosurface of baryon density + hedgehog arrows
# ============================================================

def fig_3d_isosurface():
    print("Figure 3: 3D isosurface...")
    N3d = 100
    L = 4.0
    h = 2*L / N3d
    x = np.linspace(-L + h/2, L - h/2, N3d)
    X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')
    R3 = np.sqrt(X3**2 + Y3**2 + Z3**2)

    # Topological distance: 1-cos(f) -> 0 in vacuum, 2 at anti-vacuum core
    topo = 1.0 - np.cos(f_of_r(R3))

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Outer isosurface at 1-s = 0.3 (outer halo)
    verts, faces, _, _ = marching_cubes(topo, level=0.3, spacing=(h, h, h))
    verts = verts - np.array([L, L, L])
    mesh = Poly3DCollection(verts[faces], alpha=0.12, linewidths=0,
                            facecolor='#91bfdb', edgecolor='none')
    ax.add_collection3d(mesh)

    # Middle isosurface at 1-s = 1.0 (S3 equator, f=pi/2)
    verts2, faces2, _, _ = marching_cubes(topo, level=1.0, spacing=(h, h, h))
    verts2 = verts2 - np.array([L, L, L])
    mesh2 = Poly3DCollection(verts2[faces2], alpha=0.3, linewidths=0,
                             facecolor='#fc8d59', edgecolor='none')
    ax.add_collection3d(mesh2)

    # Inner isosurface at 1-s = 1.8 (deep core)
    verts3, faces3, _, _ = marching_cubes(topo, level=1.8, spacing=(h, h, h))
    verts3 = verts3 - np.array([L, L, L])
    mesh3 = Poly3DCollection(verts3[faces3], alpha=0.65, linewidths=0,
                             facecolor='#d73027', edgecolor='none')
    ax.add_collection3d(mesh3)

    # Hedgehog arrows on two spherical shells
    for r_shell, col, lw in [(2.2, '#2166ac', 1.2), (1.0, '#1a9850', 1.4)]:
        n_ph = 12
        n_th = 7
        phi_arr = np.linspace(0, 2*np.pi, n_ph, endpoint=False)
        theta_arr = np.linspace(0.12*np.pi, 0.88*np.pi, n_th)
        for th in theta_arr:
            for ph in phi_arr:
                px = r_shell * np.sin(th) * np.cos(ph)
                py = r_shell * np.sin(th) * np.sin(ph)
                pz = r_shell * np.cos(th)
                f_val = float(f_of_r(r_shell))
                sf = np.sin(f_val)
                scale = 0.45
                ax.quiver(px, py, pz,
                          scale * sf * np.sin(th)*np.cos(ph),
                          scale * sf * np.sin(th)*np.sin(ph),
                          scale * sf * np.cos(th),
                          color=col, arrow_length_ratio=0.3, linewidth=lw)

    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.set_zlim(-3.2, 3.2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(r'B = 1 Skyrmion:  topological structure + pion hedgehog',
                 fontsize=13, fontweight='bold', pad=15)
    ax.view_init(elev=25, azim=135)

    ax.text2D(0.02, 0.96,
              'Isosurfaces of 1 - cos f(r):\n'
              '  Blue: outer halo (f ~ 33$\\degree$)\n'
              '  Orange: S$^3$ equator (f = 90$\\degree$)\n'
              '  Red: anti-vacuum core (f ~ 154$\\degree$)',
              transform=ax.transAxes, fontsize=9, va='top',
              bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

    path = os.path.join(OUTDIR, 'soliton_3d.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close()

# ============================================================
# Figure 4: 3D two-soliton isosurface
# ============================================================

def fig_3d_two_soliton():
    print("Figure 4: 3D two-soliton isosurface...")
    N3d = 128
    L = 6.0
    h = 2*L / N3d
    x = np.linspace(-L + h/2, L - h/2, N3d)
    X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')

    D = 4.0
    q = product_ansatz(X3, Y3, Z3, D)
    # Topological distance from vacuum: 1-s = 1-cos(f) -> 0 in vacuum, 2 at core
    topo_dist = 1.0 - q[0]

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Outer isosurface at 1-s = 0.5 (f ~ 60 deg)
    try:
        verts, faces, _, _ = marching_cubes(topo_dist, level=0.5, spacing=(h, h, h))
        verts = verts - np.array([L, L, L])
        mesh = Poly3DCollection(verts[faces], alpha=0.15, linewidths=0,
                                facecolor='#91bfdb', edgecolor='none')
        ax.add_collection3d(mesh)
    except Exception as e:
        print(f"  Outer isosurface failed: {e}")

    # Middle isosurface at 1-s = 1.0 (f = pi/2, equator of S3)
    try:
        verts2, faces2, _, _ = marching_cubes(topo_dist, level=1.0, spacing=(h, h, h))
        verts2 = verts2 - np.array([L, L, L])
        mesh2 = Poly3DCollection(verts2[faces2], alpha=0.35, linewidths=0,
                                 facecolor='#fc8d59', edgecolor='none')
        ax.add_collection3d(mesh2)
    except Exception as e:
        print(f"  Middle isosurface failed: {e}")

    # Inner isosurface at 1-s = 1.8 (deep core, near anti-vacuum)
    try:
        verts3, faces3, _, _ = marching_cubes(topo_dist, level=1.8, spacing=(h, h, h))
        verts3 = verts3 - np.array([L, L, L])
        mesh3 = Poly3DCollection(verts3[faces3], alpha=0.7, linewidths=0,
                                 facecolor='#d73027', edgecolor='none')
        ax.add_collection3d(mesh3)
    except Exception as e:
        print(f"  Inner isosurface failed: {e}")

    lim = 4.5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(f'Two-soliton product ansatz  (D = {D})',
                 fontsize=13, fontweight='bold', pad=15)
    ax.view_init(elev=25, azim=135)

    ax.text2D(0.02, 0.96,
              'Isosurfaces of 1 - s(x):\n'
              '  Outer (blue): s = 0.5\n'
              '  Middle (orange): s = 0 (S$^3$ equator)\n'
              '  Inner (red): s = -0.8 (near anti-vacuum)',
              transform=ax.transAxes, fontsize=9, va='top', family='monospace',
              bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

    path = os.path.join(OUTDIR, 'soliton_3d_two_body.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close()

# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("Skyrmion 3D Visualization")
    print("=" * 50)
    fig_single_soliton()
    fig_two_soliton()
    fig_3d_isosurface()
    fig_3d_two_soliton()
    print("\nAll figures saved to:", OUTDIR)
