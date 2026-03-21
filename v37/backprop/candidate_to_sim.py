"""
Convert a backprop candidate (.npy) to an initial condition for the C simulator.

Upscales from training resolution (N=32/48) to simulation resolution (N=128)
using trilinear interpolation, then writes a raw binary file that v37_knot.c
can load.

Usage:
    python candidate_to_sim.py best_phi.npy best_vel.npy --N_out 128 --L 15 -o candidate_init.bin
"""

import numpy as np
from scipy.ndimage import zoom
import argparse
import struct
import sys


def upscale_field(phi, N_out):
    """Trilinear interpolation from phi's resolution to N_out."""
    N_in = phi.shape[1]
    if N_in == N_out:
        return phi.copy()
    factor = N_out / N_in
    result = np.zeros((3, N_out, N_out, N_out), dtype=np.float64)
    for a in range(3):
        result[a] = zoom(phi[a].astype(np.float64), factor, order=1)
    return result


def write_init_bin(phi, vel, N, L, path):
    """Write initial condition in the format v37 simulations can read.

    Format: int N, double L, double t=0, int nfields=6,
            then 3 phi arrays + 3 theta arrays (all zeros for theta),
            then 3 phi_vel arrays + 3 theta_vel arrays.
    """
    N3 = N * N * N
    with open(path, 'wb') as f:
        f.write(struct.pack('i', N))
        f.write(struct.pack('d', L))
        f.write(struct.pack('d', 0.0))  # t=0
        f.write(struct.pack('i', 6))    # 6 fields (3 phi + 3 theta)

        # phi fields
        for a in range(3):
            f.write(phi[a].astype(np.float64).tobytes())
        # theta fields (zeros)
        zeros = np.zeros(N3, dtype=np.float64)
        for a in range(3):
            f.write(zeros.tobytes())

        # phi velocities
        for a in range(3):
            f.write(vel[a].astype(np.float64).tobytes())
        # theta velocities (zeros)
        for a in range(3):
            f.write(zeros.tobytes())

    size_gb = (12 * N3 * 8 + 20) / 1e9
    print(f"Wrote {path}: N={N}, L={L}, {size_gb:.2f} GB")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phi_path', help='Path to best_phi.npy')
    parser.add_argument('vel_path', help='Path to best_vel.npy')
    parser.add_argument('--N_out', type=int, default=128)
    parser.add_argument('--L', type=float, default=15.0)
    parser.add_argument('-o', '--output', default='candidate_init.bin')
    args = parser.parse_args()

    phi = np.load(args.phi_path)
    vel = np.load(args.vel_path)
    N_in = phi.shape[1]

    print(f"Input: {phi.shape}, N_in={N_in}")
    print(f"Target: N_out={args.N_out}")

    if N_in != args.N_out:
        print(f"Upscaling {N_in} -> {args.N_out}...")
        phi = upscale_field(phi, args.N_out)
        vel = upscale_field(vel, args.N_out)

    # Stats
    for a in range(3):
        P = phi[0] * phi[1] * phi[2]
        print(f"phi_{a}: [{phi[a].min():.4f}, {phi[a].max():.4f}] std={phi[a].std():.4f}")
    print(f"max|P| = {np.max(np.abs(P)):.4f}")

    write_init_bin(phi, vel, args.N_out, args.L, args.output)


if __name__ == '__main__':
    main()
