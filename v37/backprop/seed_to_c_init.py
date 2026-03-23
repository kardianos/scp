"""
Convert a backprop candidate (.npz) to an initial condition binary for v37_seedrun.c.

Loads phi (3x32x32x32) and vel (3x32x32x32) from an .npz file, upscales to N=128
using trilinear interpolation, adds zero theta fields, and writes a raw binary file.

Binary format:
    int N, double L, double t=0, int nfields=12
    Then 12 arrays of float64 in row-major order (each N^3 doubles):
      phi0, phi1, phi2, theta0, theta1, theta2,
      phi0_vel, phi1_vel, phi2_vel, theta0_vel, theta1_vel, theta2_vel

Usage:
    python seed_to_c_init.py candidate.npz --N_out 128 --L 15 -o seed_init.bin
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
        return phi.copy().astype(np.float64)
    factor = N_out / N_in
    result = np.zeros((3, N_out, N_out, N_out), dtype=np.float64)
    for a in range(3):
        result[a] = zoom(phi[a].astype(np.float64), factor, order=1)
    return result


def write_seed_bin(phi, vel, N, L, path):
    """Write seed binary in the format v37_seedrun.c expects.

    Header: int N, double L, double t=0, int nfields=12
    Data (12 arrays, each N^3 float64, row-major):
      phi0, phi1, phi2,           (position fields)
      theta0, theta1, theta2,     (angle fields, zeros)
      phi0_vel, phi1_vel, phi2_vel,   (position velocities)
      theta0_vel, theta1_vel, theta2_vel  (angle velocities, zeros)
    """
    N3 = N * N * N
    zeros = np.zeros(N3, dtype=np.float64)

    with open(path, 'wb') as f:
        # Header
        f.write(struct.pack('i', N))
        f.write(struct.pack('d', L))
        f.write(struct.pack('d', 0.0))  # t=0
        f.write(struct.pack('i', 12))   # 12 arrays total

        # phi fields (3)
        for a in range(3):
            f.write(phi[a].astype(np.float64).tobytes())

        # theta fields (3, all zeros)
        for a in range(3):
            f.write(zeros.tobytes())

        # phi velocities (3)
        for a in range(3):
            f.write(vel[a].astype(np.float64).tobytes())

        # theta velocities (3, all zeros)
        for a in range(3):
            f.write(zeros.tobytes())

    header_bytes = 4 + 8 + 8 + 4  # N + L + t + nfields = 24 bytes
    data_bytes = 12 * N3 * 8
    total_bytes = header_bytes + data_bytes
    size_mb = total_bytes / 1e6
    print(f"Wrote {path}: N={N}, L={L}, nfields=12, {size_mb:.1f} MB ({total_bytes} bytes)")


def main():
    parser = argparse.ArgumentParser(
        description='Convert backprop candidate .npz to v37_seedrun binary')
    parser.add_argument('npz_path', help='Path to candidate .npz file (phi: 3xNxNxN, vel: 3xNxNxN)')
    parser.add_argument('--N_out', type=int, default=128,
                        help='Output grid resolution (default: 128)')
    parser.add_argument('--L', type=float, default=15.0,
                        help='Half-box size (default: 15.0)')
    parser.add_argument('-o', '--output', default='seed_init.bin',
                        help='Output binary file path')
    parser.add_argument('--phi_key', default='phi',
                        help='Key for phi array in .npz (default: phi)')
    parser.add_argument('--vel_key', default='vel',
                        help='Key for velocity array in .npz (default: vel)')
    args = parser.parse_args()

    # Load .npz
    data = np.load(args.npz_path)
    available_keys = list(data.keys())
    print(f"Loaded {args.npz_path}, keys: {available_keys}")

    if args.phi_key not in data:
        print(f"ERROR: key '{args.phi_key}' not found. Available: {available_keys}")
        sys.exit(1)
    if args.vel_key not in data:
        print(f"ERROR: key '{args.vel_key}' not found. Available: {available_keys}")
        sys.exit(1)

    phi = data[args.phi_key]
    vel = data[args.vel_key]

    # Validate shapes
    if phi.ndim != 4 or phi.shape[0] != 3:
        print(f"ERROR: phi shape {phi.shape}, expected (3, N, N, N)")
        sys.exit(1)
    if vel.ndim != 4 or vel.shape[0] != 3:
        print(f"ERROR: vel shape {vel.shape}, expected (3, N, N, N)")
        sys.exit(1)

    N_in = phi.shape[1]
    if phi.shape != (3, N_in, N_in, N_in):
        print(f"ERROR: phi not cubic: {phi.shape}")
        sys.exit(1)
    if vel.shape != (3, N_in, N_in, N_in):
        print(f"ERROR: vel not cubic: {vel.shape}")
        sys.exit(1)

    print(f"Input: phi {phi.shape}, vel {vel.shape}, N_in={N_in}")
    print(f"Target: N_out={args.N_out}, L={args.L}")

    # Upscale if needed
    if N_in != args.N_out:
        print(f"Upscaling {N_in} -> {args.N_out} (trilinear)...")
        phi = upscale_field(phi, args.N_out)
        vel = upscale_field(vel, args.N_out)
        print(f"Upscaled: phi {phi.shape}, vel {vel.shape}")
    else:
        phi = phi.astype(np.float64)
        vel = vel.astype(np.float64)

    # Statistics
    for a in range(3):
        print(f"  phi_{a}: [{phi[a].min():.4f}, {phi[a].max():.4f}] "
              f"rms={np.sqrt(np.mean(phi[a]**2)):.4f}")
    for a in range(3):
        print(f"  vel_{a}: [{vel[a].min():.4f}, {vel[a].max():.4f}] "
              f"rms={np.sqrt(np.mean(vel[a]**2)):.4f}")
    P = phi[0] * phi[1] * phi[2]
    print(f"  max|P| = {np.max(np.abs(P)):.6f}")
    print(f"  int|P| = {np.sum(np.abs(P)) * (2*args.L/args.N_out)**3:.4f}")

    write_seed_bin(phi, vel, args.N_out, args.L, args.output)


if __name__ == '__main__':
    main()
