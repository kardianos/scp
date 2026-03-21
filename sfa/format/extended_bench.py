#!/usr/bin/env python3
"""Extended SFA Compression Benchmark

Tests:
1. BSS-only (no compression) — entropy analysis per byte stream
2. Alternative codecs: lz4, brotli, snappy, gzip, lzma (with and without BSS)

Uses same synthetic data generation as compression_bench.py for fair comparison.
"""

import numpy as np
import time
import sys
import math

import pyzstd
import lz4.frame
import brotli
import snappy
import gzip
import lzma

# ============================================================
# DATA GENERATION (same as compression_bench.py)
# ============================================================

def generate_synthetic_frame(N=64, L=15.0, t=0.0, seed=42):
    """Generate one frame of 6 float64 fields mimicking Cosserat simulation."""
    rng = np.random.RandomState(seed)
    dx = 2*L / (N-1)
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
        columns.append(field.astype(np.float64).ravel())

    return columns


def frame_to_bytes(columns):
    """Concatenate column arrays into single byte buffer."""
    parts = []
    for col in columns:
        parts.append(col.tobytes())
    return b''.join(parts)


# ============================================================
# BSS PRIMITIVES
# ============================================================

def bss_encode(data_bytes, n_values, elem_size):
    """Byte Stream Split: transpose bytes for float compression."""
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(n_values, elem_size)
    transposed = arr.T.copy()  # (elem_size, n_values)
    return transposed.tobytes()


def bss_decode(data_bytes, n_values, elem_size):
    """Reverse Byte Stream Split."""
    arr = np.frombuffer(data_bytes, dtype=np.uint8).reshape(elem_size, n_values)
    return arr.T.copy().tobytes()


# ============================================================
# ENTROPY ANALYSIS
# ============================================================

def byte_entropy(data):
    """Compute Shannon entropy of byte stream in bits/byte."""
    arr = np.frombuffer(data, dtype=np.uint8)
    counts = np.bincount(arr, minlength=256)
    probs = counts / len(arr)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))


def analyze_bss_entropy(raw_data, n_values, elem_size=8):
    """Analyze entropy of each BSS byte stream separately."""
    bss = bss_encode(raw_data, n_values, elem_size)
    stream_size = n_values  # bytes per stream

    results = []
    for b in range(elem_size):
        stream = bss[b * stream_size : (b + 1) * stream_size]
        h = byte_entropy(stream)
        # Theoretical compressed size for this stream (bits -> bytes)
        theo_bytes = h * stream_size / 8.0
        results.append({
            'stream': b,
            'entropy': h,
            'raw_bytes': stream_size,
            'theoretical_bytes': theo_bytes,
        })

    return results


# ============================================================
# CODEC WRAPPERS (compress/decompress pairs)
# ============================================================

def make_zstd(level):
    def compress(data):
        return pyzstd.compress(data, level)
    def decompress(data):
        return pyzstd.decompress(data)
    return compress, decompress, f"zstd-{level}"


def make_lz4(level=0):
    """lz4 frame compression. level=0 is default (fast)."""
    def compress(data):
        return lz4.frame.compress(data, compression_level=level)
    def decompress(data):
        return lz4.frame.decompress(data)
    label = f"lz4" if level == 0 else f"lz4-{level}"
    return compress, decompress, label


def make_lz4_hc(level=9):
    """lz4 high-compression mode."""
    def compress(data):
        return lz4.frame.compress(data, compression_level=level)
    def decompress(data):
        return lz4.frame.decompress(data)
    return compress, decompress, f"lz4hc-{level}"


def make_brotli(quality):
    def compress(data):
        return brotli.compress(data, quality=quality)
    def decompress(data):
        return brotli.decompress(data)
    return compress, decompress, f"brotli-{quality}"


def make_snappy():
    def compress(data):
        return snappy.compress(data)
    def decompress(data):
        return snappy.decompress(data)
    return compress, decompress, "snappy"


def make_gzip(level):
    def compress(data):
        return gzip.compress(data, compresslevel=level)
    def decompress(data):
        return gzip.decompress(data)
    return compress, decompress, f"gzip-{level}"


def make_lzma(preset):
    def compress(data):
        return lzma.compress(data, preset=preset)
    def decompress(data):
        return lzma.decompress(data)
    return compress, decompress, f"lzma-{preset}"


# ============================================================
# BENCHMARK HARNESS
# ============================================================

def bench_pipeline(name, raw_data, n_values, elem_size, use_bss,
                   compress_fn, decompress_fn, n_iters=5):
    """Benchmark a single pipeline configuration.

    Returns dict with ratio, compress/decompress speeds, and verification status.
    """
    raw_size = len(raw_data)

    # --- Compress ---
    # Warmup
    if use_bss:
        bss_data = bss_encode(raw_data, n_values, elem_size)
        compressed = compress_fn(bss_data)
    else:
        compressed = compress_fn(raw_data)

    # Timed runs
    t0 = time.perf_counter()
    for _ in range(n_iters):
        if use_bss:
            bss_data = bss_encode(raw_data, n_values, elem_size)
            compressed = compress_fn(bss_data)
        else:
            compressed = compress_fn(raw_data)
    t_comp = (time.perf_counter() - t0) / n_iters

    comp_size = len(compressed)

    # --- Decompress ---
    # Warmup
    if use_bss:
        dec_bss = decompress_fn(compressed)
        decompressed = bss_decode(dec_bss, n_values, elem_size)
    else:
        decompressed = decompress_fn(compressed)

    # Timed runs
    t0 = time.perf_counter()
    for _ in range(n_iters):
        if use_bss:
            dec_bss = decompress_fn(compressed)
            decompressed = bss_decode(dec_bss, n_values, elem_size)
        else:
            decompressed = decompress_fn(compressed)
    t_decomp = (time.perf_counter() - t0) / n_iters

    # Verify roundtrip
    verified = (decompressed == raw_data)

    ratio = raw_size / comp_size if comp_size > 0 else float('inf')
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


def bench_bss_only(raw_data, n_values, elem_size=8, n_iters=10):
    """Benchmark BSS alone (no compression). Measures transpose throughput."""
    raw_size = len(raw_data)

    # Warmup
    bss_data = bss_encode(raw_data, n_values, elem_size)
    _ = bss_decode(bss_data, n_values, elem_size)

    # Encode timing
    t0 = time.perf_counter()
    for _ in range(n_iters):
        bss_data = bss_encode(raw_data, n_values, elem_size)
    t_enc = (time.perf_counter() - t0) / n_iters

    # Decode timing
    t0 = time.perf_counter()
    for _ in range(n_iters):
        decoded = bss_decode(bss_data, n_values, elem_size)
    t_dec = (time.perf_counter() - t0) / n_iters

    # Verify roundtrip
    verified = (decoded == raw_data)

    enc_speed = (raw_size / 1e6) / t_enc if t_enc > 0 else 0
    dec_speed = (raw_size / 1e6) / t_dec if t_dec > 0 else 0

    # BSS doesn't change size
    assert len(bss_data) == raw_size, "BSS should preserve data size"

    return {
        'name': 'BSS only',
        'ratio': 1.0,
        'comp_MB_s': enc_speed,
        'decomp_MB_s': dec_speed,
        'comp_size': raw_size,
        'verified': verified,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 80)
    print("EXTENDED SFA COMPRESSION BENCHMARK")
    print("=" * 80)

    # Generate test data (same as compression_bench.py)
    print("\nGenerating synthetic test data (64^3 x 6 fields, float64)...")
    columns = generate_synthetic_frame(N=64, L=15.0, t=0.0, seed=42)
    raw_data = frame_to_bytes(columns)
    n_values = sum(len(c) for c in columns)
    elem_size = 8
    raw_size = len(raw_data)

    print(f"Raw size: {raw_size} bytes ({raw_size/1e6:.2f} MB)")
    print(f"Values: {n_values} float64")

    # ============================================================
    # PART 1: BSS-only analysis — entropy per byte stream
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 1: BSS ENTROPY ANALYSIS")
    print("=" * 80)

    # Raw (interleaved) entropy
    raw_entropy = byte_entropy(raw_data)
    print(f"\nRaw data entropy: {raw_entropy:.3f} bits/byte")
    print(f"Theoretical minimum size (raw): {raw_entropy * raw_size / 8:.0f} bytes "
          f"(ratio {8.0/raw_entropy:.2f}x)")

    # BSS per-stream entropy
    print(f"\nBSS per-stream entropy (byte 0 = most significant / exponent, byte 7 = least significant):")
    print(f"{'Stream':>8} {'Entropy':>10} {'Theo size':>12} {'Theo ratio':>12} {'Description':>20}")
    print("-" * 70)

    stream_results = analyze_bss_entropy(raw_data, n_values, elem_size)
    total_theo = 0
    for sr in stream_results:
        b = sr['stream']
        h = sr['entropy']
        theo = sr['theoretical_bytes']
        total_theo += theo
        ratio_str = f"{8.0/h:.2f}x" if h > 0 else "inf"

        if b == 0:
            desc = "sign+exponent MSB"
        elif b == 1:
            desc = "exponent LSB+mantissa"
        elif b < 4:
            desc = f"mantissa (high)"
        else:
            desc = f"mantissa (low)"

        print(f"{b:>8} {h:>10.3f} {theo:>12.0f} {ratio_str:>12} {desc:>20}")

    total_raw_bits = raw_size * 8
    total_theo_bits = sum(sr['entropy'] * sr['raw_bytes'] for sr in stream_results)
    overall_theo_ratio = raw_size / (total_theo_bits / 8.0)

    print(f"\n--- Summary ---")
    print(f"Total theoretical compressed size (per-stream optimal): {total_theo/1e6:.2f} MB "
          f"(ratio {overall_theo_ratio:.3f}x)")
    print(f"Best achieved (BSS+zstd-22 from previous study): ratio 1.55x")
    print(f"Headroom: theoretical {overall_theo_ratio:.3f}x vs achieved 1.55x "
          f"=> zstd captures {1.55/overall_theo_ratio*100:.0f}% of theoretical gain"
          if overall_theo_ratio > 1 else "")

    # BSS-only throughput
    bss_result = bench_bss_only(raw_data, n_values, elem_size, n_iters=10)
    print(f"\nBSS transpose throughput: "
          f"encode {bss_result['comp_MB_s']:.0f} MB/s, "
          f"decode {bss_result['decomp_MB_s']:.0f} MB/s")
    print(f"Roundtrip verified: {bss_result['verified']}")

    # ============================================================
    # PART 2: Alternative codecs
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 2: ALTERNATIVE CODECS")
    print("=" * 80)

    # Define all codec configurations
    codecs = [
        # (compress_fn, decompress_fn, label)
        make_zstd(1),
        make_zstd(3),
        make_zstd(6),
        make_zstd(9),
        make_zstd(15),
        make_zstd(22),
        make_lz4(0),       # default (fast) mode
        make_lz4_hc(9),    # high compression
        make_brotli(1),
        make_brotli(4),
        make_brotli(6),
        make_brotli(9),
        make_brotli(11),   # max quality
        make_snappy(),
        make_gzip(1),
        make_gzip(6),
        make_gzip(9),
        make_lzma(1),
        make_lzma(3),
        make_lzma(6),
    ]

    all_results = []

    # Add raw baseline
    all_results.append({
        'name': 'raw (no compress)',
        'ratio': 1.0,
        'comp_MB_s': float('inf'),
        'decomp_MB_s': float('inf'),
        'comp_size': raw_size,
        'verified': True,
    })

    # Add BSS-only baseline
    all_results.append(bss_result)

    n_iters = 5  # iterations per measurement

    for comp_fn, decomp_fn, label in codecs:
        # Without BSS
        print(f"  Testing {label} (raw)...", end='', flush=True)
        try:
            r = bench_pipeline(
                label, raw_data, n_values, elem_size,
                use_bss=False, compress_fn=comp_fn, decompress_fn=decomp_fn,
                n_iters=n_iters
            )
            all_results.append(r)
            print(f" ratio={r['ratio']:.2f}x, "
                  f"comp={r['comp_MB_s']:.0f} MB/s, "
                  f"decomp={r['decomp_MB_s']:.0f} MB/s"
                  f"{'' if r['verified'] else ' VERIFY FAILED!'}")
        except Exception as e:
            print(f" FAILED: {e}")
            all_results.append({
                'name': label,
                'ratio': 0, 'comp_MB_s': 0, 'decomp_MB_s': 0,
                'comp_size': 0, 'verified': False,
            })

        # With BSS
        print(f"  Testing BSS + {label}...", end='', flush=True)
        try:
            r = bench_pipeline(
                f"BSS + {label}", raw_data, n_values, elem_size,
                use_bss=True, compress_fn=comp_fn, decompress_fn=decomp_fn,
                n_iters=n_iters
            )
            all_results.append(r)
            print(f" ratio={r['ratio']:.2f}x, "
                  f"comp={r['comp_MB_s']:.0f} MB/s, "
                  f"decomp={r['decomp_MB_s']:.0f} MB/s"
                  f"{'' if r['verified'] else ' VERIFY FAILED!'}")
        except Exception as e:
            print(f" FAILED: {e}")
            all_results.append({
                'name': f"BSS + {label}",
                'ratio': 0, 'comp_MB_s': 0, 'decomp_MB_s': 0,
                'comp_size': 0, 'verified': False,
            })

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    print("\n" + "=" * 80)
    print("FULL RESULTS TABLE")
    print("=" * 80)

    print(f"\n{'Pipeline':<28} {'Ratio':>7} {'Comp MB/s':>11} {'Decomp MB/s':>13} {'Verified':>9}")
    print("-" * 72)
    for r in all_results:
        ratio_s = f"{r['ratio']:.2f}" if r['ratio'] < 1000 else "1.00"
        comp_s = f"{r['comp_MB_s']:.0f}" if r['comp_MB_s'] < 1e9 else "inf"
        decomp_s = f"{r['decomp_MB_s']:.0f}" if r['decomp_MB_s'] < 1e9 else "inf"
        ver_s = "yes" if r['verified'] else "FAIL"
        print(f"{r['name']:<28} {ratio_s:>7} {comp_s:>11} {decomp_s:>13} {ver_s:>9}")

    # ============================================================
    # COMPARISON: Codec-only vs BSS+Codec
    # ============================================================
    print("\n" + "=" * 80)
    print("BSS IMPROVEMENT FACTOR (ratio with BSS / ratio without BSS)")
    print("=" * 80)

    # Group by base codec name
    codec_names = set()
    by_name = {}
    for r in all_results:
        by_name[r['name']] = r

    print(f"\n{'Codec':<20} {'Raw ratio':>10} {'BSS ratio':>10} {'BSS factor':>11} {'BSS comp':>10} {'BSS decomp':>12}")
    print("-" * 78)

    for comp_fn, decomp_fn, label in codecs:
        raw_r = by_name.get(label)
        bss_r = by_name.get(f"BSS + {label}")
        if raw_r and bss_r and raw_r['ratio'] > 0 and bss_r['ratio'] > 0:
            factor = bss_r['ratio'] / raw_r['ratio']
            print(f"{label:<20} {raw_r['ratio']:>10.2f} {bss_r['ratio']:>10.2f} "
                  f"{factor:>10.2f}x {bss_r['comp_MB_s']:>10.0f} {bss_r['decomp_MB_s']:>12.0f}")

    # ============================================================
    # PARETO FRONT: ratio vs decompress speed
    # ============================================================
    print("\n" + "=" * 80)
    print("PARETO ANALYSIS: Best ratio-vs-speed tradeoffs")
    print("=" * 80)

    # Filter valid BSS results
    bss_results = [r for r in all_results if r['name'].startswith('BSS + ') and r['ratio'] > 1.0]
    bss_results.sort(key=lambda r: r['ratio'], reverse=True)

    print(f"\nBSS + codec results, sorted by ratio (best to worst):")
    print(f"{'Pipeline':<28} {'Ratio':>7} {'Comp MB/s':>11} {'Decomp MB/s':>13}")
    print("-" * 62)
    for r in bss_results:
        print(f"{r['name']:<28} {r['ratio']:>7.3f} {r['comp_MB_s']:>11.0f} {r['decomp_MB_s']:>13.0f}")

    # Identify Pareto front (non-dominated in ratio AND decompress speed)
    print(f"\nPareto-optimal pipelines (not dominated in both ratio and decompress speed):")
    pareto = []
    for r in bss_results:
        dominated = False
        for other in bss_results:
            if other is r:
                continue
            if other['ratio'] >= r['ratio'] and other['decomp_MB_s'] >= r['decomp_MB_s']:
                if other['ratio'] > r['ratio'] or other['decomp_MB_s'] > r['decomp_MB_s']:
                    dominated = True
                    break
        if not dominated:
            pareto.append(r)

    print(f"{'Pipeline':<28} {'Ratio':>7} {'Comp MB/s':>11} {'Decomp MB/s':>13}")
    print("-" * 62)
    for r in pareto:
        print(f"{r['name']:<28} {r['ratio']:>7.3f} {r['comp_MB_s']:>11.0f} {r['decomp_MB_s']:>13.0f}")

    # ============================================================
    # GENERATE MARKDOWN TABLE
    # ============================================================

    # Build the table in the requested format
    print("\n" + "=" * 80)
    print("MARKDOWN TABLE (for COMPRESSION_STUDY.md)")
    print("=" * 80)

    # Curated list of interesting pipelines
    curated_names = [
        'raw (no compress)',
        'BSS only',
        # zstd
        'zstd-1', 'zstd-3', 'zstd-6', 'zstd-9', 'zstd-22',
        'BSS + zstd-1', 'BSS + zstd-3', 'BSS + zstd-6', 'BSS + zstd-9', 'BSS + zstd-22',
        # lz4
        'lz4', 'BSS + lz4',
        'lz4hc-9', 'BSS + lz4hc-9',
        # brotli
        'brotli-1', 'BSS + brotli-1',
        'brotli-4', 'BSS + brotli-4',
        'brotli-6', 'BSS + brotli-6',
        'brotli-9', 'BSS + brotli-9',
        'brotli-11', 'BSS + brotli-11',
        # snappy
        'snappy', 'BSS + snappy',
        # gzip
        'gzip-1', 'BSS + gzip-1',
        'gzip-6', 'BSS + gzip-6',
        'gzip-9', 'BSS + gzip-9',
        # lzma
        'lzma-1', 'BSS + lzma-1',
        'lzma-3', 'BSS + lzma-3',
        'lzma-6', 'BSS + lzma-6',
    ]

    # Entropy analysis table
    print("\nEntropy per BSS byte stream:")
    print("| Stream | Description | Entropy (bits/byte) | Theoretical ratio | Notes |")
    print("|--------|-------------|--------------------:|------------------:|-------|")
    for sr in stream_results:
        b = sr['stream']
        h = sr['entropy']
        if b == 0:
            desc = "sign + exponent MSB"
            note = "very low entropy, highly compressible"
        elif b == 1:
            desc = "exponent LSB + mantissa[51:48]"
            note = "low entropy"
        elif b == 2:
            desc = "mantissa[47:40]"
            note = "moderate entropy"
        elif b < 5:
            desc = f"mantissa[{47-8*(b-2)}:{40-8*(b-2)}]"
            note = "high entropy"
        else:
            desc = f"mantissa (low bits)"
            note = "near-random, incompressible"

        theo_r = 8.0 / h if h > 0 else float('inf')
        print(f"| {b} | {desc} | {h:.3f} | {theo_r:.2f}x | {note} |")

    print(f"\nTotal raw entropy: {raw_entropy:.3f} bits/byte (theoretical ratio: {8.0/raw_entropy:.2f}x)")
    print(f"Per-stream optimal total: {total_theo/1e6:.3f} MB from {raw_size/1e6:.2f} MB "
          f"(ratio {overall_theo_ratio:.3f}x)")

    # Main comparison table
    print("\nFull codec comparison:")
    print("| Pipeline | Ratio | Compress MB/s | Decompress MB/s | Notes |")
    print("|----------|------:|:-------------:|:---------------:|-------|")

    for name in curated_names:
        r = by_name.get(name)
        if r is None:
            continue

        if r['comp_MB_s'] > 1e9:
            comp_s = "---"
        else:
            comp_s = f"{r['comp_MB_s']:.0f}"

        if r['decomp_MB_s'] > 1e9:
            decomp_s = "---"
        else:
            decomp_s = f"{r['decomp_MB_s']:.0f}"

        ratio_s = f"{r['ratio']:.3f}" if r['ratio'] < 100 else "1.000"

        # Notes
        note = ""
        if name == 'raw (no compress)':
            note = "baseline"
        elif name == 'BSS only':
            note = "transpose only, no size change"
        elif name == 'BSS + zstd-3':
            note = "**current SFA default**"
        elif name == 'BSS + zstd-1':
            note = "fast zstd"
        elif 'lz4' in name and 'hc' not in name and 'BSS' in name:
            note = "fastest codec"
        elif 'snappy' in name and 'BSS' in name:
            note = "Google, speed-optimized"
        elif 'lzma' in name and 'BSS' in name:
            note = "slow but high ratio"
        elif 'brotli-11' in name and 'BSS' in name:
            note = "max brotli quality"

        if not r['verified'] and r['ratio'] > 0:
            note = "ROUNDTRIP FAILED"

        print(f"| {name} | {ratio_s} | {comp_s} | {decomp_s} | {note} |")

    return all_results, stream_results, overall_theo_ratio


if __name__ == '__main__':
    all_results, stream_results, theo_ratio = main()
