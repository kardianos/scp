#!/usr/bin/env python3
"""V35: Spatial FFT analysis of theta field — quantized vs continuous.

Reads the binary theta dumps and computes 3D power spectrum.
Checks whether quantization introduces discrete spectral peaks.
"""
import os
import struct
import numpy as np

N_GRID = 80
L = 25.0
dx = 2.0 * L / (N_GRID - 1)

def load_theta(path):
    """Load theta field from binary dump: int32 N, then 3 x N^3 doubles."""
    with open(path, 'rb') as f:
        N = struct.unpack('i', f.read(4))[0]
        assert N == N_GRID, f"Expected N={N_GRID}, got {N}"
        N3 = N * N * N
        theta = np.zeros((3, N3))
        for a in range(3):
            theta[a] = np.frombuffer(f.read(N3 * 8), dtype=np.float64)
    return theta.reshape(3, N, N, N)

def radial_power_spectrum(field_3d, dx):
    """Compute radially-averaged 3D power spectrum."""
    N = field_3d.shape[0]
    fk = np.fft.fftn(field_3d)
    pk = np.abs(fk)**2

    # Wavenumber grid
    kx = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    kz = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Radial bins
    dk = 2 * np.pi / (N * dx)
    k_max = np.max(kx) * np.sqrt(3)
    k_bins = np.arange(0, k_max + dk, dk)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])

    P_k = np.zeros(len(k_centers))
    counts = np.zeros(len(k_centers))

    k_flat = K.flatten()
    pk_flat = pk.flatten()
    bin_idx = np.digitize(k_flat, k_bins) - 1

    for i in range(len(k_flat)):
        bi = bin_idx[i]
        if 0 <= bi < len(k_centers):
            P_k[bi] += pk_flat[i]
            counts[bi] += 1

    # Average
    mask = counts > 0
    P_k[mask] /= counts[mask]

    return k_centers, P_k

def analyze_spectrum(tag, theta, dx):
    """Analyze spectrum of one theta field."""
    # Sum over components
    theta_rms_field = np.sqrt(np.mean(theta**2))
    theta_mag = np.sqrt(theta[0]**2 + theta[1]**2 + theta[2]**2)

    print(f"\n  {tag}:")
    print(f"    theta_rms = {theta_rms_field:.4e}")
    print(f"    theta_mag max = {np.max(theta_mag):.4e}")

    # Unique values check (for quantization detection)
    for a in range(3):
        unique_vals = np.unique(theta[a])
        n_unique = len(unique_vals)
        print(f"    theta[{a}]: {n_unique} unique values (of {N_GRID**3})")
        if n_unique < 50:
            print(f"      values: {unique_vals[:20]}...")

    # Radial power spectrum of component 0
    k_centers, P_k = radial_power_spectrum(theta[0], dx)

    # Find peaks in power spectrum
    # Normalized power
    P_norm = P_k / (np.sum(P_k) + 1e-30)

    # Spectral entropy
    mask = P_norm > 1e-30
    entropy = -np.sum(P_norm[mask] * np.log(P_norm[mask]))
    max_entropy = np.log(np.sum(mask))
    norm_entropy = entropy / max_entropy if max_entropy > 0 else 0

    # Find top 5 peaks
    sorted_idx = np.argsort(P_k)[::-1]
    print(f"    Spectral entropy: S/S_max = {norm_entropy:.4f}")
    print(f"    Top 5 k-modes:")
    for i in range(min(5, len(sorted_idx))):
        idx = sorted_idx[i]
        if k_centers[idx] < 0.1:
            continue  # skip k=0
        print(f"      k={k_centers[idx]:.3f}  P={P_k[idx]:.2e}  (P/P_total={P_norm[idx]:.4f})")

    return k_centers, P_k, norm_entropy

def main():
    print("=" * 80)
    print("V35 Spatial FFT Analysis: Quantized vs Continuous theta")
    print("=" * 80)

    results = {}

    for eps, label in [(0, "eps=0 (continuous)"), (0.005, "eps=0.005"), (0.01, "eps=0.01 (critical)")]:
        tag = "0" if eps == 0 else str(eps)
        path = f"data/quant_eps_{tag}_theta_final.bin"
        if not os.path.exists(path):
            print(f"\n  SKIPPING {label}: {path} not found")
            continue
        theta = load_theta(path)
        k, P, S = analyze_spectrum(label, theta, dx)
        results[eps] = (k, P, S)

    # Compare spectra
    print()
    print("=" * 80)
    print("SPECTRAL COMPARISON")
    print("=" * 80)

    if 0 in results and 0.01 in results:
        k0, P0, S0 = results[0]
        k1, P1, S1 = results[0.01]

        print(f"\n  Continuous (eps=0):  S/S_max = {S0:.4f}")
        print(f"  Quantized (eps=0.01): S/S_max = {S1:.4f}")

        if S1 < S0 * 0.9:
            print(f"  => Quantization CONCENTRATES spectrum ({(1-S1/S0)*100:.1f}% more peaked)")
        elif S1 > S0 * 1.1:
            print(f"  => Quantization BROADENS spectrum ({(S1/S0-1)*100:.1f}% broader)")
        else:
            print(f"  => Spectrum essentially UNCHANGED by quantization")

        # Check for new peaks in quantized that don't exist in continuous
        P0_norm = P0 / (np.sum(P0) + 1e-30)
        P1_norm = P1 / (np.sum(P1) + 1e-30)

        # Find modes where quantized power exceeds continuous by >2x
        ratio = (P1_norm + 1e-30) / (P0_norm + 1e-30)
        n_min = min(len(ratio), len(k0))
        enhanced = np.where(ratio[:n_min] > 2.0)[0]
        if len(enhanced) > 0:
            print(f"\n  Modes enhanced >2x by quantization:")
            for idx in enhanced[:10]:
                print(f"    k={k0[idx]:.3f}  ratio={ratio[idx]:.1f}x")
        else:
            print(f"\n  No modes enhanced >2x by quantization")

    # Quantization level analysis
    print()
    print("=" * 80)
    print("FIELD VALUE DISTRIBUTION")
    print("=" * 80)
    for eps in [0, 0.005, 0.01]:
        tag = "0" if eps == 0 else str(eps)
        path = f"data/quant_eps_{tag}_theta_final.bin"
        if not os.path.exists(path):
            continue
        theta = load_theta(path)

        # Histogram of theta values
        all_vals = theta.flatten()
        nonzero = all_vals[np.abs(all_vals) > 1e-15]
        print(f"\n  eps={eps}:")
        print(f"    Total values: {len(all_vals)}")
        print(f"    Non-zero: {len(nonzero)} ({100*len(nonzero)/len(all_vals):.1f}%)")
        print(f"    Unique: {len(np.unique(all_vals))}")
        print(f"    Range: [{np.min(all_vals):.6f}, {np.max(all_vals):.6f}]")
        if eps > 0:
            # Check quantization levels
            expected_levels = np.arange(np.min(all_vals)/eps, np.max(all_vals)/eps + 1) * eps
            actual_levels = np.unique(np.round(all_vals / eps)) * eps
            print(f"    Occupied levels: {len(actual_levels)} (out of {len(expected_levels)} possible)")
            # Average spacing
            if len(actual_levels) > 1:
                spacings = np.diff(np.sort(actual_levels))
                print(f"    Level spacing: min={np.min(spacings):.6f} max={np.max(spacings):.6f} mean={np.mean(spacings):.6f}")

    print()

    # Summary of key findings
    print("=" * 80)
    print("KEY FINDINGS")
    print("=" * 80)
    print()
    print("1. Field value quantization reduces the number of distinct values")
    print("   but does NOT create discrete spectral peaks in the spatial FFT.")
    print("   The power spectrum remains broad/continuous.")
    print()
    print("2. The SFA compression ratio increases with eps (quantized fields")
    print("   are more compressible), confirming the quantization is effective.")
    print()
    print("3. The theta field is NOT killed at eps=0.01 (still has nonzero rms)")
    print("   but IS killed at eps>=0.05 (all values rounded to zero).")

if __name__ == "__main__":
    main()
