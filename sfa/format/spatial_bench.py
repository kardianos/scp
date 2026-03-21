#!/usr/bin/env python3
"""Spatial Ordering Compression Benchmark for SFA

Tests whether spatial reordering of 3D field data improves compression
when combined with BSS + zstd.

Orderings tested:
1. Linear (row-major, ix slowest, iz fastest) — current SFA default
2. Z-order / Morton curve (bit-interleaving of ix,iy,iz)
3. Hilbert curve (better locality than Morton)
4. Plane-major (same as linear for cubic data, but tested for completeness)
5. 3D delta encoding: x-delta, y-delta, z-delta (subtract along each axis)

All tests use the same synthetic data (64^3 x 6 float64) as the other benchmarks.
"""

import numpy as np
import time

import pyzstd

# ============================================================
# DATA GENERATION (same as compression_bench.py)
# ============================================================

def generate_synthetic_frame(N=64, L=15.0, t=0.0, seed=42):
    """Generate one frame of 6 float64 fields mimicking Cosserat simulation."""
    rng = np.random.RandomState(seed)
    x = np.linspace(-L, L, N)
    y = np.linspace(-L, L, N)
    z = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    columns = []
    for col_idx in range(6):
        A_bg = 0.1
        k = 2*np.pi / (2*L) * (1 + 0.3*col_idx)
        phase = 0.5 * col_idx + 0.1*t
        bg = A_bg * np.cos(k*Z + phase)

        A_core = 0.8
        sigma = 2.0
        r2 = X**2 + Y**2 + Z**2
        envelope = np.exp(-r2 / (2*sigma**2))
        core_phase = k*Z + 0.7*col_idx + 0.3*t
        core = A_core * envelope * np.sin(core_phase)

        noise = rng.normal(0, 1e-6, (N, N, N))

        field = bg + core + noise
        columns.append(field.astype(np.float64))  # keep as 3D array

    return columns


# ============================================================
# BSS PRIMITIVES
# ============================================================

def bss_encode(data_bytes, n_values, elem_size):
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(n_values, elem_size)
    return arr.T.copy().tobytes()

def bss_decode(data_bytes, n_values, elem_size):
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(elem_size, n_values)
    return arr.T.copy().tobytes()


# ============================================================
# SPATIAL ORDERINGS
# ============================================================

def part1by2(n):
    """Spread bits of n into every 3rd bit position for Morton encoding.
    Works for n up to 2^10 = 1024 (enough for N=64..512)."""
    n = np.asarray(n, dtype=np.uint64)
    n = n & 0x000003FF            # 10 bits
    n = (n | (n << 16)) & 0x030000FF
    n = (n | (n << 8))  & 0x0300F00F
    n = (n | (n << 4))  & 0x030C30C3
    n = (n | (n << 2))  & 0x09249249
    return n


def morton_encode_3d(ix, iy, iz):
    """Compute Morton code (Z-order) for 3D coordinates."""
    return part1by2(ix) | (part1by2(iy) << 1) | (part1by2(iz) << 2)


def morton_order(N):
    """Return permutation array: morton_order[i] = linear index of i-th Morton element.
    i.e., data[morton_order] reorders data from linear to Morton order."""
    ix = np.arange(N, dtype=np.uint64)
    iy = np.arange(N, dtype=np.uint64)
    iz = np.arange(N, dtype=np.uint64)

    # All 3D indices
    IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing='ij')
    IX = IX.ravel()
    IY = IY.ravel()
    IZ = IZ.ravel()

    # Linear indices (row-major)
    linear = IX * N * N + IY * N + IZ

    # Morton codes
    morton = morton_encode_3d(IX, IY, IZ)

    # Sort by Morton code to get the reordering
    sort_idx = np.argsort(morton)

    # perm[i] = linear index of the i-th element in Morton order
    return linear[sort_idx].astype(np.intp)


# Hilbert curve implementation (iterative, 3D)
def hilbert_d2xyz(d, order):
    """Convert Hilbert distance d to 3D coordinates (x,y,z) for a cube of side 2^order.
    Based on the algorithm from Wikipedia / Skilling's paper."""
    x = y = z = 0
    s = 1
    for i in range(order):
        rx = 1 if (d & 4) else 0
        ry = 1 if (d & 2) else 0
        rz = 1 if (d & 1) else 0

        # Rotate
        if rz == 0:
            if ry == 0:
                if rx == 0:
                    pass  # no rotation
                else:
                    x, y = s - 1 - y, s - 1 - x
            else:
                if rx == 0:
                    x, z = z, x
                else:
                    y, z = s - 1 - z, s - 1 - y

        x += s * rx
        y += s * ry
        z += s * rz

        d >>= 3
        s <<= 1

    return x, y, z


def hilbert_order(N):
    """Return permutation array for Hilbert curve ordering.
    N must be a power of 2."""
    assert N & (N - 1) == 0, f"N={N} must be a power of 2 for Hilbert curve"
    order = int(np.log2(N))
    n_total = N * N * N

    # Build lookup: hilbert distance -> linear index
    perm = np.empty(n_total, dtype=np.intp)
    for d in range(n_total):
        x, y, z = hilbert_d2xyz(d, order)
        linear = x * N * N + y * N + z
        perm[d] = linear

    return perm


def hilbert_order_fast(N):
    """Faster Hilbert curve using table-based approach.
    Uses the Skilling transpose algorithm."""
    assert N & (N - 1) == 0, f"N={N} must be a power of 2"
    order = int(np.log2(N))
    n_total = N**3

    # Pre-compute all coordinates via vectorized approach
    # Generate all Hilbert distances
    distances = np.arange(n_total, dtype=np.int64)

    # Convert to coordinates using the compact Hilbert algorithm
    # We use a simplified bit-manipulation approach
    coords = np.zeros((n_total, 3), dtype=np.int64)

    for d_idx in range(n_total):
        d = int(distances[d_idx])
        x = y = z = 0
        s = 1
        dd = d
        for i in range(order):
            rx = 1 if (dd & 4) else 0
            ry = 1 if (dd & 2) else 0
            rz = 1 if (dd & 1) else 0

            if rz == 0:
                if ry == 0:
                    if rx == 0:
                        pass
                    else:
                        x, y = s - 1 - y, s - 1 - x
                else:
                    if rx == 0:
                        x, z = z, x
                    else:
                        y, z = s - 1 - z, s - 1 - y

            x += s * rx
            y += s * ry
            z += s * rz
            dd >>= 3
            s <<= 1

        coords[d_idx] = [x, y, z]

    perm = coords[:, 0] * N * N + coords[:, 1] * N + coords[:, 2]
    return perm.astype(np.intp)


# ============================================================
# 3D DELTA ENCODING
# ============================================================

def delta_z_encode(field_3d):
    """Delta along z (fastest-varying in row-major). field_3d is (Nx,Ny,Nz)."""
    out = np.empty_like(field_3d)
    out[:, :, 0] = field_3d[:, :, 0]
    out[:, :, 1:] = field_3d[:, :, 1:] - field_3d[:, :, :-1]
    return out

def delta_z_decode(delta_3d):
    return np.cumsum(delta_3d, axis=2)


def delta_y_encode(field_3d):
    """Delta along y (middle axis)."""
    out = np.empty_like(field_3d)
    out[:, 0, :] = field_3d[:, 0, :]
    out[:, 1:, :] = field_3d[:, 1:, :] - field_3d[:, :-1, :]
    return out

def delta_y_decode(delta_3d):
    return np.cumsum(delta_3d, axis=1)


def delta_x_encode(field_3d):
    """Delta along x (slowest-varying in row-major)."""
    out = np.empty_like(field_3d)
    out[0, :, :] = field_3d[0, :, :]
    out[1:, :, :] = field_3d[1:, :, :] - field_3d[:-1, :, :]
    return out

def delta_x_decode(delta_3d):
    return np.cumsum(delta_3d, axis=0)


def idelta_z_encode(field_bytes, N):
    """Integer-domain delta along z (LOSSLESS). Operates on int64 view of float64."""
    arr = np.frombuffer(field_bytes, dtype=np.int64).reshape(N, N, N)
    out = np.empty_like(arr)
    out[:, :, 0] = arr[:, :, 0]
    out[:, :, 1:] = arr[:, :, 1:] - arr[:, :, :-1]
    return out.tobytes()

def idelta_z_decode(delta_bytes, N):
    arr = np.frombuffer(delta_bytes, dtype=np.int64).reshape(N, N, N).copy()
    for iz in range(1, N):
        arr[:, :, iz] += arr[:, :, iz-1]
    return arr.tobytes()


def idelta_y_encode(field_bytes, N):
    """Integer-domain delta along y (LOSSLESS)."""
    arr = np.frombuffer(field_bytes, dtype=np.int64).reshape(N, N, N)
    out = np.empty_like(arr)
    out[:, 0, :] = arr[:, 0, :]
    out[:, 1:, :] = arr[:, 1:, :] - arr[:, :-1, :]
    return out.tobytes()

def idelta_y_decode(delta_bytes, N):
    arr = np.frombuffer(delta_bytes, dtype=np.int64).reshape(N, N, N).copy()
    for iy in range(1, N):
        arr[:, iy, :] += arr[:, iy-1, :]
    return arr.tobytes()


def idelta_x_encode(field_bytes, N):
    """Integer-domain delta along x (LOSSLESS)."""
    arr = np.frombuffer(field_bytes, dtype=np.int64).reshape(N, N, N)
    out = np.empty_like(arr)
    out[0, :, :] = arr[0, :, :]
    out[1:, :, :] = arr[1:, :, :] - arr[:-1, :, :]
    return out.tobytes()

def idelta_x_decode(delta_bytes, N):
    arr = np.frombuffer(delta_bytes, dtype=np.int64).reshape(N, N, N).copy()
    for ix in range(1, N):
        arr[ix, :, :] += arr[ix-1, :, :]
    return arr.tobytes()


# ============================================================
# BENCHMARK
# ============================================================

def bench_pipeline(name, raw_data, n_values, elem_size, preprocess, postprocess,
                   zstd_level=3, n_iters=5):
    """Benchmark: preprocess -> BSS -> zstd compress, and reverse.

    preprocess: raw_bytes -> preprocessed_bytes
    postprocess: preprocessed_bytes -> raw_bytes
    """
    raw_size = len(raw_data)

    # Warmup
    pp = preprocess(raw_data)
    bss = bss_encode(pp, n_values, elem_size)
    comp = pyzstd.compress(bss, zstd_level)

    # Compress timing
    t0 = time.perf_counter()
    for _ in range(n_iters):
        pp = preprocess(raw_data)
        bss = bss_encode(pp, n_values, elem_size)
        comp = pyzstd.compress(bss, zstd_level)
    t_comp = (time.perf_counter() - t0) / n_iters

    comp_size = len(comp)

    # Decompress timing
    t0 = time.perf_counter()
    for _ in range(n_iters):
        dec_bss = pyzstd.decompress(comp)
        dec_pp = bss_decode(dec_bss, n_values, elem_size)
        dec_raw = postprocess(dec_pp)
    t_decomp = (time.perf_counter() - t0) / n_iters

    # Verify
    verified = (dec_raw == raw_data)

    ratio = raw_size / comp_size
    comp_speed = (raw_size / 1e6) / t_comp if t_comp > 0 else 0
    decomp_speed = (raw_size / 1e6) / t_decomp if t_decomp > 0 else 0

    return {
        'name': name,
        'ratio': ratio,
        'comp_MB_s': comp_speed,
        'decomp_MB_s': decomp_speed,
        'comp_size': comp_size,
        'verified': verified,
    }


def main():
    N = 64
    L = 15.0

    print("=" * 80)
    print("SPATIAL ORDERING COMPRESSION BENCHMARK")
    print("=" * 80)

    # Generate test data
    print(f"\nGenerating synthetic test data ({N}^3 x 6 fields, float64)...")
    columns_3d = generate_synthetic_frame(N=N, L=L, t=0.0, seed=42)

    # Concatenate all columns in linear (row-major) order
    raw_data = b''.join(col.ravel().tobytes() for col in columns_3d)
    n_values = 6 * N**3
    elem_size = 8
    raw_size = len(raw_data)
    print(f"Raw size: {raw_size} bytes ({raw_size/1e6:.2f} MB)")

    results = []

    # ============================================================
    # 1. Linear order (baseline) — identity preprocess
    # ============================================================
    print("\n--- Linear order (current SFA default) ---")
    identity = lambda x: x
    r = bench_pipeline("linear + BSS + zstd-3", raw_data, n_values, elem_size,
                       identity, identity)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # ============================================================
    # 2. Morton / Z-order
    # ============================================================
    print("\n--- Morton (Z-order) curve ---")
    print("  Computing Morton permutation...", flush=True)
    morton_perm = morton_order(N)
    # Inverse permutation for decode
    morton_inv = np.empty_like(morton_perm)
    morton_inv[morton_perm] = np.arange(len(morton_perm))

    def morton_preprocess(data):
        # Reorder each column (N^3 float64 values) by Morton order
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64)
            parts.append(arr[morton_perm].tobytes())
        return b''.join(parts)

    def morton_postprocess(data):
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64)
            parts.append(arr[morton_inv].tobytes())
        return b''.join(parts)

    r = bench_pipeline("morton + BSS + zstd-3", raw_data, n_values, elem_size,
                       morton_preprocess, morton_postprocess)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # ============================================================
    # 3. Hilbert curve
    # ============================================================
    print("\n--- Hilbert curve ---")
    print("  Computing Hilbert permutation (this may take a moment)...", flush=True)
    t0 = time.perf_counter()
    hilbert_perm = hilbert_order_fast(N)
    t_hilbert = time.perf_counter() - t0
    print(f"  Hilbert permutation computed in {t_hilbert:.1f}s")

    # Verify it's a valid permutation
    assert len(np.unique(hilbert_perm)) == N**3, "Hilbert permutation has duplicates!"

    hilbert_inv = np.empty_like(hilbert_perm)
    hilbert_inv[hilbert_perm] = np.arange(len(hilbert_perm))

    def hilbert_preprocess(data):
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64)
            parts.append(arr[hilbert_perm].tobytes())
        return b''.join(parts)

    def hilbert_postprocess(data):
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64)
            parts.append(arr[hilbert_inv].tobytes())
        return b''.join(parts)

    r = bench_pipeline("hilbert + BSS + zstd-3", raw_data, n_values, elem_size,
                       hilbert_preprocess, hilbert_postprocess)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # ============================================================
    # 4. Plane-major (xy-planes stacked in z)
    # ============================================================
    print("\n--- Plane-major (xy stacked in z) ---")
    # For row-major (ix,iy,iz) with iz fastest, this IS the natural order.
    # But we test yz-planes stacked in x, which transposes the axis order.

    def plane_x_preprocess(data):
        """Reorder to (iz, iy, ix) — z-planes stacked in x."""
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64).reshape(N, N, N)
            # Transpose: (ix,iy,iz) -> (iz,iy,ix) so z varies slowest
            reordered = arr.transpose(2, 1, 0).copy()
            parts.append(reordered.ravel().tobytes())
        return b''.join(parts)

    def plane_x_postprocess(data):
        n_per_col = N**3
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            arr = np.frombuffer(col_bytes, dtype=np.float64).reshape(N, N, N)
            reordered = arr.transpose(2, 1, 0).copy()
            parts.append(reordered.ravel().tobytes())
        return b''.join(parts)

    r = bench_pipeline("plane-x + BSS + zstd-3", raw_data, n_values, elem_size,
                       plane_x_preprocess, plane_x_postprocess)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # ============================================================
    # 5. 3D Delta Encoding (integer-domain, lossless)
    # ============================================================
    print("\n--- 3D Integer Delta Encoding ---")

    n_per_col = N**3

    # z-delta (along fastest-varying axis — same as linear idelta)
    def idelta_z_pre(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_z_encode(col_bytes, N))
        return b''.join(parts)

    def idelta_z_post(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_z_decode(col_bytes, N))
        return b''.join(parts)

    r = bench_pipeline("idelta-z + BSS + zstd-3", raw_data, n_values, elem_size,
                       idelta_z_pre, idelta_z_post)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # y-delta
    def idelta_y_pre(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_y_encode(col_bytes, N))
        return b''.join(parts)

    def idelta_y_post(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_y_decode(col_bytes, N))
        return b''.join(parts)

    r = bench_pipeline("idelta-y + BSS + zstd-3", raw_data, n_values, elem_size,
                       idelta_y_pre, idelta_y_post)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # x-delta
    def idelta_x_pre(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_x_encode(col_bytes, N))
        return b''.join(parts)

    def idelta_x_post(data):
        parts = []
        for c in range(6):
            col_bytes = data[c * n_per_col * 8 : (c+1) * n_per_col * 8]
            parts.append(idelta_x_decode(col_bytes, N))
        return b''.join(parts)

    r = bench_pipeline("idelta-x + BSS + zstd-3", raw_data, n_values, elem_size,
                       idelta_x_pre, idelta_x_post)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # Also test z-delta with higher zstd levels to see if spatial delta helps more there
    for level in [1, 9, 22]:
        r = bench_pipeline(f"idelta-z + BSS + zstd-{level}", raw_data, n_values, elem_size,
                           idelta_z_pre, idelta_z_post, zstd_level=level)
        results.append(r)
        print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
              f"decomp={r['decomp_MB_s']:.0f} MB/s (zstd-{level}), verified={r['verified']}")

    # Test: Morton + z-delta combined
    print("\n--- Morton + z-delta combined ---")
    def morton_idelta_z_pre(data):
        # First reorder to Morton, then z-delta
        morton_data = morton_preprocess(data)
        return idelta_z_pre(morton_data)

    def morton_idelta_z_post(data):
        undelta = idelta_z_post(data)
        return morton_postprocess(undelta)

    r = bench_pipeline("morton + idelta-z + BSS + zstd-3", raw_data, n_values, elem_size,
                       morton_idelta_z_pre, morton_idelta_z_post)
    results.append(r)
    print(f"  ratio={r['ratio']:.3f}x, comp={r['comp_MB_s']:.0f}, "
          f"decomp={r['decomp_MB_s']:.0f} MB/s, verified={r['verified']}")

    # ============================================================
    # Analyze spatial correlations
    # ============================================================
    print("\n--- Spatial Correlation Analysis ---")
    col0 = columns_3d[0]  # first field, 3D array

    # Mean absolute difference along each axis
    dz = np.mean(np.abs(col0[:, :, 1:] - col0[:, :, :-1]))
    dy = np.mean(np.abs(col0[:, 1:, :] - col0[:, :-1, :]))
    dx = np.mean(np.abs(col0[1:, :, :] - col0[:-1, :, :]))

    print(f"  Mean |delta| along z (fastest): {dz:.6f}")
    print(f"  Mean |delta| along y (middle):  {dy:.6f}")
    print(f"  Mean |delta| along x (slowest): {dx:.6f}")
    print(f"  Ratio dx/dz: {dx/dz:.2f}, dy/dz: {dy/dz:.2f}")

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    print("\n" + "=" * 80)
    print("SPATIAL ORDERING RESULTS")
    print("=" * 80)

    print(f"\n{'Pipeline':<40} {'Ratio':>7} {'Comp MB/s':>11} {'Decomp MB/s':>13} {'Verified':>9}")
    print("-" * 84)
    for r in results:
        print(f"{r['name']:<40} {r['ratio']:>7.3f} {r['comp_MB_s']:>11.0f} "
              f"{r['decomp_MB_s']:>13.0f} {('yes' if r['verified'] else 'FAIL'):>9}")

    # Markdown table
    print("\n--- MARKDOWN TABLE ---")
    print("| Pipeline | Ratio | Compress MB/s | Decompress MB/s | Lossless | Notes |")
    print("|----------|------:|:-------------:|:---------------:|:--------:|-------|")

    notes_map = {
        "linear + BSS + zstd-3": "**current SFA default**",
        "morton + BSS + zstd-3": "Z-order space-filling curve",
        "hilbert + BSS + zstd-3": "Hilbert space-filling curve",
        "plane-x + BSS + zstd-3": "yz-planes stacked in x",
        "idelta-z + BSS + zstd-3": "delta along z (fastest axis)",
        "idelta-y + BSS + zstd-3": "delta along y (middle axis)",
        "idelta-x + BSS + zstd-3": "delta along x (slowest axis)",
        "idelta-z + BSS + zstd-1": "z-delta, fast zstd",
        "idelta-z + BSS + zstd-9": "z-delta, slow zstd",
        "idelta-z + BSS + zstd-22": "z-delta, max zstd",
        "morton + idelta-z + BSS + zstd-3": "Morton reorder + z-delta",
    }

    for r in results:
        note = notes_map.get(r['name'], "")
        ll = "yes" if r['verified'] else "NO"
        print(f"| {r['name']} | {r['ratio']:.3f} | {r['comp_MB_s']:.0f} | "
              f"{r['decomp_MB_s']:.0f} | {ll} | {note} |")

    return results


if __name__ == '__main__':
    results = main()
