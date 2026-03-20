#!/usr/bin/env python3
"""
Visualize a 6-field Cosserat snapshot: position (φ) + angle (θ) fields.

Three channels mapped to RGB:
  RED:   bound field — |P| = |φ₀φ₁φ₂| (triple product, braid core)
  GREEN: unbound field — Σφ² in regions where |P| is small (background fabric)
  BLUE:  angular field — Σθ² (rotation/torsion excitation)

Usage: python3 visualize.py <snapshot.bin> [slice_axis] [output.png]
  slice_axis: 'x', 'y', or 'z' (default: 'z' = xy-plane through braid center)
"""

import sys
import struct
import numpy as np

def read_snapshot(path):
    with open(path, 'rb') as f:
        N = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        t = struct.unpack('d', f.read(8))[0]
        nf = struct.unpack('i', f.read(4))[0]

        fields = []
        for _ in range(nf):
            data = np.frombuffer(f.read(N**3 * 8), dtype=np.float64)
            fields.append(data.reshape(N, N, N))

    return N, L, t, nf, fields

def find_braid_center(phi, N):
    """Find braid center as energy-weighted centroid above 5× average."""
    phi2 = sum(p**2 for p in phi[:3])
    avg = phi2.mean()
    thresh = 5.0 * avg
    mask = phi2 > thresh

    if mask.sum() == 0:
        return N//2, N//2, N//2

    coords = np.mgrid[0:N, 0:N, 0:N]
    weights = phi2 * mask
    total = weights.sum()
    ci = int(round((coords[0] * weights).sum() / total))
    cj = int(round((coords[1] * weights).sum() / total))
    ck = int(round((coords[2] * weights).sum() / total))
    return ci, cj, ck

def make_slice_image(fields, N, L, center, axis='z'):
    """Create RGB image from a 2D slice through the braid center."""
    phi = fields[:3]   # position fields
    theta = fields[3:] if len(fields) >= 6 else [np.zeros_like(fields[0])]*3

    ci, cj, ck = center

    # Extract 2D slices
    if axis == 'z':
        sl = np.s_[:, :, ck]
        xlabel, ylabel = 'x', 'y'
    elif axis == 'y':
        sl = np.s_[:, cj, :]
        xlabel, ylabel = 'x', 'z'
    else:  # x
        sl = np.s_[ci, :, :]
        xlabel, ylabel = 'y', 'z'

    # Compute quantities on the slice
    phi2 = sum(p[sl]**2 for p in phi)             # total field energy
    P = phi[0][sl] * phi[1][sl] * phi[2][sl]       # triple product
    absP = np.abs(P)                                # binding strength
    theta2 = sum(t[sl]**2 for t in theta)           # angular excitation

    # Normalize each channel
    def normalize(arr, percentile=99.5):
        vmax = np.percentile(arr, percentile)
        if vmax < 1e-30:
            return np.zeros_like(arr)
        return np.clip(arr / vmax, 0, 1)

    # RED: bound field (|P| = braid core)
    red = normalize(absP)

    # GREEN: unbound field (Σφ² where |P| is small = background fabric)
    # Use Σφ² weighted by (1 - |P|/|P|_max) to suppress the core
    P_norm = absP / (absP.max() + 1e-30)
    unbound = phi2 * (1.0 - P_norm**0.5)  # suppress core
    green = normalize(unbound)

    # BLUE: angular field (Σθ²)
    blue = normalize(theta2)

    # Compose RGB
    img = np.stack([red, green, blue], axis=-1)

    # Apply gamma for visibility
    img = np.power(img, 0.4)

    return img, xlabel, ylabel

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 visualize.py <snapshot.bin> [axis] [output.png]")
        sys.exit(1)

    path = sys.argv[1]
    axis = sys.argv[2] if len(sys.argv) > 2 else 'z'
    outpath = sys.argv[3] if len(sys.argv) > 3 else None

    print(f"Reading {path}...")
    N, L, t, nf, fields = read_snapshot(path)
    print(f"  N={N}, L={L}, t={t}, nf={nf}")

    center = find_braid_center(fields, N)
    print(f"  Braid center: i={center[0]}, j={center[1]}, k={center[2]}")

    # Compute field statistics
    phi2 = sum(f**2 for f in fields[:3])
    absP = np.abs(fields[0] * fields[1] * fields[2])
    if nf >= 6:
        theta2 = sum(f**2 for f in fields[3:6])
        print(f"  Σφ²: mean={phi2.mean():.4e}, max={phi2.max():.4e}")
        print(f"  |P|:  mean={absP.mean():.4e}, max={absP.max():.4e}")
        print(f"  Σθ²: mean={theta2.mean():.4e}, max={theta2.max():.4e}")
    else:
        print(f"  Σφ²: mean={phi2.mean():.4e}, max={phi2.max():.4e}")
        print(f"  |P|:  mean={absP.mean():.4e}, max={absP.max():.4e}")
        print("  No θ fields in this snapshot")

    # Make images for all three axes
    for ax in ['x', 'y', 'z']:
        img, xl, yl = make_slice_image(fields, N, L, center, ax)

        # Save as raw PPM (no matplotlib needed)
        h, w, _ = img.shape
        img_bytes = (img * 255).astype(np.uint8)

        fname = outpath if (outpath and ax == axis) else f"slice_{ax}_t{int(t)}.ppm"
        with open(fname, 'wb') as f:
            f.write(f"P6\n{w} {h}\n255\n".encode())
            f.write(img_bytes.tobytes())
        print(f"  Saved {fname} ({w}×{h}, {ax}-slice, {xl}-{yl} plane)")

    # Also save a combined profile (1D radial)
    print("\n  Radial profiles (from braid center):")
    ci, cj, ck = center
    dx = 2*L/(N-1)
    max_r = min(L, 30)
    n_bins = int(max_r / 0.5)
    bins_phi2 = np.zeros(n_bins)
    bins_absP = np.zeros(n_bins)
    bins_theta2 = np.zeros(n_bins)
    bins_count = np.zeros(n_bins)

    for i in range(N):
        for j in range(N):
            for k in range(N):
                r = np.sqrt(((i-ci)*dx)**2 + ((j-cj)*dx)**2 + ((k-ck)*dx)**2)
                b = int(r / 0.5)
                if b >= n_bins:
                    continue
                idx = i*N*N + j*N + k
                p2 = sum(fields[a].flat[idx]**2 for a in range(3))
                ap = abs(fields[0].flat[idx] * fields[1].flat[idx] * fields[2].flat[idx])
                t2 = sum(fields[a].flat[idx]**2 for a in range(3, min(6, nf))) if nf >= 6 else 0
                bins_phi2[b] += p2
                bins_absP[b] += ap
                bins_theta2[b] += t2
                bins_count[b] += 1

    # Save profile
    prof_path = f"radial_profile_t{int(t)}.tsv"
    with open(prof_path, 'w') as f:
        f.write("r\tphi2\tabsP\ttheta2\tcounts\n")
        for b in range(n_bins):
            if bins_count[b] > 0:
                r = (b + 0.5) * 0.5
                f.write(f"{r:.2f}\t{bins_phi2[b]/bins_count[b]:.6e}\t"
                        f"{bins_absP[b]/bins_count[b]:.6e}\t"
                        f"{bins_theta2[b]/bins_count[b]:.6e}\t"
                        f"{int(bins_count[b])}\n")
    print(f"  Saved {prof_path}")

    print("\n  Color key:")
    print("    RED   = bound field (|P| = |φ₀φ₁φ₂|, braid core)")
    print("    GREEN = unbound field (Σφ² where |P| small, fabric)")
    print("    BLUE  = angular field (Σθ², rotation/torsion)")

if __name__ == '__main__':
    main()
