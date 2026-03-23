#!/usr/bin/env python3
"""
Generate initial conditions for a toroidal braid.

Takes the linear braid3 (three phase-shifted fields on a helical tube along z)
and bends it into a torus so the wave travels continuously around the ring.
All three fields are co-located on ONE ring with phase offsets.

Single torus lies in the xy-plane at z=0.
Multi-torus: torus 1 in xy-plane, torus 2 in xz-plane, torus 3 in yz-plane.

Output: .npz with keys phi (3,N,N,N), vel (3,N,N,N), N, L.
"""

import numpy as np
import argparse
import sys


# --- Default parameters ---
R_MAJOR  = 15.0                    # torus major radius
R_MINOR  = 3.0                     # tube Gaussian envelope width
R_HELIX  = 1.0                     # strand helix radius within tube
N_OSC    = 3                       # oscillations around the torus (integer for closure)
N_TWIST  = 3                       # strand twists around the tube
A_AMP    = 0.8                     # amplitude
DELTA    = [0.0, 3.0005, 4.4325]  # phase offsets (V28 optimized)
M2       = 2.25                    # mass^2 for omega
N_GRID   = 256                     # grid resolution
L_BOX    = 30.0                    # half-box size (grid spans [-L, L])


def build_torus_braid(N, L, R_major, R_minor, R_helix, n_osc, n_twist,
                      A, delta, m2, chirality_sign):
    """
    Build phi (3,N,N,N) and vel (3,N,N,N) for a single torus braid
    lying in the xy-plane centered at the origin.

    Returns (phi, vel) as float64 arrays.
    """
    dx = 2.0 * L / N
    coords = np.linspace(-L + dx / 2, L - dx / 2, N)  # cell centers

    # 1D coordinate arrays
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')

    # Torus geometry
    r_cyl = np.sqrt(X**2 + Y**2)             # cylindrical radius
    theta_tor = np.arctan2(Y, X)              # toroidal angle [0, 2pi)
    d_tube = np.sqrt((r_cyl - R_major)**2 + Z**2)  # distance to torus center-circle

    # Poloidal angle in tube cross-section: angle from the outward-radial direction
    # in the (r_cyl - R_major, z) plane
    psi = np.arctan2(Z, r_cyl - R_major)     # poloidal angle

    inv2r2 = 1.0 / (2.0 * R_minor**2)

    # Wave parameters
    k = n_osc / R_major                       # effective wavenumber
    omega = np.sqrt(k**2 + m2)
    omega_signed = chirality_sign * omega      # +omega for U, -omega for D

    phi = np.zeros((3, N, N, N), dtype=np.float64)
    vel = np.zeros((3, N, N, N), dtype=np.float64)

    for a in range(3):
        field_val = np.zeros((N, N, N), dtype=np.float64)
        field_vel = np.zeros((N, N, N), dtype=np.float64)

        # Phase of the traveling wave at this field component
        phase = n_osc * theta_tor + delta[a]

        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)

        # Sum over 3 strands
        for s in range(3):
            # Strand s center in tube cross-section at poloidal angle alpha_s
            alpha_s = n_twist * theta_tor + 2.0 * np.pi * s / 3.0

            # Strand center position relative to tube center-circle:
            #   radial offset (in r_cyl direction): R_helix * cos(alpha_s)
            #   z offset: R_helix * sin(alpha_s)
            # Distance from grid point to strand center:
            dr = (r_cyl - R_major) - R_helix * np.cos(alpha_s)
            dz = Z - R_helix * np.sin(alpha_s)
            d_s2 = dr**2 + dz**2

            envelope = np.exp(-d_s2 * inv2r2)

            field_val += A * envelope * cos_phase
            field_vel += omega_signed * A * envelope * sin_phase

        phi[a] = field_val
        vel[a] = field_vel

    return phi, vel


def rotate_xz(phi, vel):
    """
    Rotate fields from xy-plane torus to xz-plane torus.
    Permute axes: (x, y, z) -> (x, z, y), i.e. swap y <-> z.
    """
    phi_rot = np.zeros_like(phi)
    vel_rot = np.zeros_like(vel)
    for a in range(3):
        phi_rot[a] = np.swapaxes(phi[a], 1, 2)
        vel_rot[a] = np.swapaxes(vel[a], 1, 2)
    return phi_rot, vel_rot


def rotate_yz(phi, vel):
    """
    Rotate fields from xy-plane torus to yz-plane torus.
    Permute axes: (x, y, z) -> (z, y, x), i.e. swap x <-> z.
    """
    phi_rot = np.zeros_like(phi)
    vel_rot = np.zeros_like(vel)
    for a in range(3):
        phi_rot[a] = np.swapaxes(phi[a], 0, 2)
        vel_rot[a] = np.swapaxes(vel[a], 0, 2)
    return phi_rot, vel_rot


def print_stats(phi, vel, L, N, label=""):
    """Print diagnostic statistics for the fields."""
    prefix = f"[{label}] " if label else ""
    dx = 2.0 * L / N
    dV = dx**3

    print(f"{prefix}Grid: N={N}, L={L}, dx={dx:.4f}")
    for a in range(3):
        rms_phi = np.sqrt(np.mean(phi[a]**2))
        rms_vel = np.sqrt(np.mean(vel[a]**2))
        print(f"{prefix}  phi_{a}: [{phi[a].min():.4f}, {phi[a].max():.4f}] rms={rms_phi:.4f}")
        print(f"{prefix}  vel_{a}: [{vel[a].min():.4f}, {vel[a].max():.4f}] rms={rms_vel:.4f}")

    P = phi[0] * phi[1] * phi[2]
    max_abs_P = np.max(np.abs(P))
    int_abs_P = np.sum(np.abs(P)) * dV
    P_pos = np.sum(P[P > 0]) * dV
    P_neg = np.sum(P[P < 0]) * dV
    print(f"{prefix}  max|P| = {max_abs_P:.6f}")
    print(f"{prefix}  int|P| = {int_abs_P:.4f}")
    print(f"{prefix}  P>0 integral = {P_pos:.4f},  P<0 integral = {P_neg:.4f}")

    # Simple energy estimate: E_pot ~ (1/2) * m^2 * sum(phi^2) * dV
    E_pot = 0.0
    for a in range(3):
        E_pot += 0.5 * M2 * np.sum(phi[a]**2) * dV
    E_kin = 0.0
    for a in range(3):
        E_kin += 0.5 * np.sum(vel[a]**2) * dV
    print(f"{prefix}  E_pot(mass) = {E_pot:.2f},  E_kin = {E_kin:.2f}")

    # Check field is nonzero in the tube
    nonzero_frac = np.sum(np.abs(phi[0]) > 1e-6) / phi[0].size
    print(f"{prefix}  Nonzero fraction (|phi_0|>1e-6): {nonzero_frac:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate toroidal braid initial conditions')
    parser.add_argument('-o', '--output', default='torus_braid.npz',
                        help='Output .npz file path')
    parser.add_argument('--R_major', type=float, default=R_MAJOR,
                        help=f'Torus major radius (default: {R_MAJOR})')
    parser.add_argument('--R_minor', type=float, default=R_MINOR,
                        help=f'Tube Gaussian width (default: {R_MINOR})')
    parser.add_argument('--R_helix', type=float, default=R_HELIX,
                        help=f'Strand helix radius (default: {R_HELIX})')
    parser.add_argument('--n_osc', type=int, default=N_OSC,
                        help=f'Oscillations around torus (default: {N_OSC})')
    parser.add_argument('--n_twist', type=int, default=N_TWIST,
                        help=f'Strand twists around tube (default: {N_TWIST})')
    parser.add_argument('--A', type=float, default=A_AMP,
                        help=f'Amplitude (default: {A_AMP})')
    parser.add_argument('--m2', type=float, default=M2,
                        help=f'Mass squared (default: {M2})')
    parser.add_argument('--N', type=int, default=N_GRID,
                        help=f'Grid resolution (default: {N_GRID})')
    parser.add_argument('--L', type=float, default=L_BOX,
                        help=f'Half-box size (default: {L_BOX})')
    parser.add_argument('--chirality', type=str, default='U',
                        help='Chirality: U (right-handed) or D (left-handed). '
                             'For multi-torus, one char per torus, e.g. UUD (default: U)')
    parser.add_argument('--n_tori', type=int, default=1,
                        help='Number of tori (1, 2, or 3)')
    parser.add_argument('--delta', type=float, nargs=3, default=DELTA,
                        help=f'Phase offsets (default: {DELTA})')
    args = parser.parse_args()

    N = args.N
    L = args.L
    n_tori = args.n_tori

    if n_tori < 1 or n_tori > 3:
        print("ERROR: n_tori must be 1, 2, or 3")
        sys.exit(1)

    # Parse chirality string
    chir_str = args.chirality.upper()
    if len(chir_str) == 1:
        chir_str = chir_str * n_tori  # replicate single char for all tori
    if len(chir_str) != n_tori:
        print(f"ERROR: chirality string '{args.chirality}' length {len(chir_str)} "
              f"!= n_tori {n_tori}")
        sys.exit(1)
    for ch in chir_str:
        if ch not in ('U', 'D'):
            print(f"ERROR: chirality must be U or D, got '{ch}'")
            sys.exit(1)

    chirality_signs = [+1.0 if ch == 'U' else -1.0 for ch in chir_str]

    print(f"Toroidal braid generator")
    print(f"  R_major={args.R_major}, R_minor={args.R_minor}, R_helix={args.R_helix}")
    print(f"  n_osc={args.n_osc}, n_twist={args.n_twist}, A={args.A}, m2={args.m2}")
    print(f"  N={N}, L={L}, dx={2*L/N:.4f}")
    print(f"  n_tori={n_tori}, chirality={chir_str}")
    print(f"  delta={args.delta}")
    print()

    phi_total = np.zeros((3, N, N, N), dtype=np.float64)
    vel_total = np.zeros((3, N, N, N), dtype=np.float64)

    for t_idx in range(n_tori):
        chir_sign = chirality_signs[t_idx]
        chir_label = chir_str[t_idx]

        print(f"--- Torus {t_idx+1}/{n_tori} (chirality={chir_label}) ---")

        phi_t, vel_t = build_torus_braid(
            N, L,
            R_major=args.R_major,
            R_minor=args.R_minor,
            R_helix=args.R_helix,
            n_osc=args.n_osc,
            n_twist=args.n_twist,
            A=args.A,
            delta=args.delta,
            m2=args.m2,
            chirality_sign=chir_sign,
        )

        # Rotate torus into correct plane
        if t_idx == 0:
            # Torus 1: xy-plane (no rotation)
            pass
        elif t_idx == 1:
            # Torus 2: xz-plane (swap y <-> z)
            phi_t, vel_t = rotate_xz(phi_t, vel_t)
        elif t_idx == 2:
            # Torus 3: yz-plane (swap x <-> z)
            phi_t, vel_t = rotate_yz(phi_t, vel_t)

        print_stats(phi_t, vel_t, L, N, label=f"torus{t_idx+1}")
        print()

        phi_total += phi_t
        vel_total += vel_t

    if n_tori > 1:
        print("--- Combined ---")
        print_stats(phi_total, vel_total, L, N, label="total")
        print()

    # Save
    np.savez(args.output,
             phi=phi_total,
             vel=vel_total,
             N=np.int32(N),
             L=np.float64(L))

    size_mb = phi_total.nbytes * 2 / 1e6  # phi + vel
    print(f"Saved {args.output}: phi {phi_total.shape}, vel {vel_total.shape}, "
          f"{size_mb:.1f} MB")


if __name__ == '__main__':
    main()
